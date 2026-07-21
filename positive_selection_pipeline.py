import os
import sys
import argparse
import subprocess
import logging
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from scipy.stats import chi2
import re
import shutil
import shlex

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def run_command(cmd, log_error=True):
    try:
        subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except subprocess.CalledProcessError as e:
        if log_error:
            logging.error(f"Command failed: {cmd}")
            logging.error(f"Error message: {e.stderr.decode('utf-8')}")
        raise

def validate_and_translate(fasta_in, prefix, out_dir):
    prot_fasta = os.path.join(out_dir, "alignments", f"{prefix}.prot.fasta")
    clean_nuc_fasta = os.path.join(out_dir, "alignments", f"{prefix}.clean_nuc.fasta")
    id_mapping_file = os.path.join(out_dir, "reports", f"{prefix}_id_mapping.txt")

    valid_records = []
    prot_records = []

    with open(fasta_in, "r") as f:
        records = list(SeqIO.parse(f, "fasta"))

    id_map = {}

    for idx, rec in enumerate(records):
        seq_str = str(rec.seq).upper()

        # PAML has a strict 30-character limit and fails unpredictably with certain special characters.
        # To guarantee 100% compatibility across MAFFT, FastTree, and PAML, we use generic, safe internal IDs.
        safe_id = f"Seq{idx+1}"
        id_map[safe_id] = str(rec.id)

        # Ensure length is multiple of 3
        if len(seq_str) % 3 != 0:
            logging.warning(f"Sequence {rec.id} length is not a multiple of 3. Truncating.")
            seq_str = seq_str[:-(len(seq_str)%3)]

        # Translate to protein
        my_seq = Seq(seq_str)
        try:
            prot_seq = my_seq.translate(table=11) # Bacterial by default
        except Exception as e:
            logging.warning(f"Translation failed for {rec.id}: {e}")
            continue

        # Handle stop codons
        if "*" in prot_seq[:-1]:
            logging.warning(f"Internal stop codon found in {rec.id}. Skipping.")
            continue

        # Remove terminal stop codon if present for codeml
        if prot_seq.endswith("*"):
            prot_seq = prot_seq[:-1]
            seq_str = seq_str[:-3]

        clean_rec = SeqRecord(Seq(seq_str), id=safe_id, description="")
        prot_rec = SeqRecord(prot_seq, id=safe_id, description="")

        valid_records.append(clean_rec)
        prot_records.append(prot_rec)

    if len(valid_records) < 3:
        raise ValueError(f"Not enough valid sequences in {fasta_in} (minimum 3 required). Found {len(valid_records)}.")

    SeqIO.write(valid_records, clean_nuc_fasta, "fasta")
    SeqIO.write(prot_records, prot_fasta, "fasta")

    # Save the mapping file so users can trace sequences back
    with open(id_mapping_file, "w") as f:
        f.write("Internal_ID\tOriginal_ID\n")
        for s_id, orig_id in id_map.items():
            f.write(f"{s_id}\t{orig_id}\n")
    logging.info(f"Sequence ID mapping saved to {id_mapping_file}")

    return clean_nuc_fasta, prot_fasta

def align_proteins(prot_fasta, prefix, out_dir):
    aln_fasta = os.path.join(out_dir, "alignments", f"{prefix}.prot.aln.fasta")
    cmd = f"mafft --auto {shlex.quote(prot_fasta)} > {shlex.quote(aln_fasta)}"
    logging.info(f"Running MAFFT: {cmd}")
    run_command(cmd)
    return aln_fasta

def back_translate(nuc_fasta, prot_aln_fasta, prefix, out_dir):
    codon_aln_fasta = os.path.join(out_dir, "alignments", f"{prefix}.codon.aln.fasta")

    nuc_dict = SeqIO.to_dict(SeqIO.parse(nuc_fasta, "fasta"))
    prot_aln_records = list(SeqIO.parse(prot_aln_fasta, "fasta"))

    codon_aln_records = []

    for prot_rec in prot_aln_records:
        nuc_seq = str(nuc_dict[prot_rec.id].seq)
        prot_aln_seq = str(prot_rec.seq)

        codon_aln = []
        nuc_idx = 0
        for aa in prot_aln_seq:
            if aa == '-':
                codon_aln.append('---')
            else:
                codon_aln.append(nuc_seq[nuc_idx:nuc_idx+3])
                nuc_idx += 3

        codon_aln_seq = "".join(codon_aln)
        codon_aln_records.append(SeqRecord(Seq(codon_aln_seq), id=prot_rec.id, description=""))

    SeqIO.write(codon_aln_records, codon_aln_fasta, "fasta")

    # Also write a phylip format for PAML
    codon_aln_phylip = os.path.join(out_dir, "alignments", f"{prefix}.codon.aln.phy")

    # Write interleaved or sequential phylip manually for strict compatibility
    with open(codon_aln_phylip, 'w') as f:
        f.write(f" {len(codon_aln_records)} {len(codon_aln_records[0].seq)}\n")
        for rec in codon_aln_records:
            # Names are already safe and <30 chars, just left justify for phylip formatting
            name = str(rec.id).ljust(30)
            f.write(f"{name}  {str(rec.seq)}\n")

    return codon_aln_fasta, codon_aln_phylip

def build_tree(codon_aln_fasta, prefix, out_dir):
    tree_file = os.path.join(out_dir, "trees", f"{prefix}.tree")
    cmd = f"FastTree -nt {shlex.quote(codon_aln_fasta)} > {shlex.quote(tree_file)}"
    logging.info(f"Running FastTree: {cmd}")
    run_command(cmd)

    # Unroot the tree for PAML
    unrooted_tree = os.path.join(out_dir, "trees", f"{prefix}.unrooted.tree")
    # A simple way to ensure the tree is unrooted for PAML is to just pass it, PAML handles it
    # We will just copy it. PAML requires the tree file to end with a semicolon. FastTree does this.
    shutil.copy(tree_file, unrooted_tree)
    return unrooted_tree

def write_codeml_ctl(ctl_file, seq_file, tree_file, out_file, model_type):
    # model_type: 7 for M7, 8 for M8
    nssites = 7 if model_type == 7 else 8

    ctl_content = f"""seqfile = {seq_file}
treefile = {tree_file}
outfile = {out_file}

noisy = 9
verbose = 1
runmode = 0
seqtype = 1
CodonFreq = 2
clock = 0
aaDist = 0
model = 0
NSsites = {nssites}
icode = 0
Mgene = 0
fix_kappa = 0
kappa = 2
fix_omega = 0
omega = 0.4
fix_alpha = 1
alpha = 0
Malpha = 0
ncatG = 10
getSE = 0
RateAncestor = 0
Small_Diff = .5e-6
cleandata = 1
fix_blength = 0
"""
    with open(ctl_file, 'w') as f:
        f.write(ctl_content)

def run_paml(aln_phy, tree_file, prefix, out_dir):
    paml_dir = os.path.join(out_dir, "paml_logs", prefix)
    os.makedirs(paml_dir, exist_ok=True)

    # When running PAML, we chdir to paml_dir, so the paths in ctl must be absolute
    abs_aln = os.path.abspath(aln_phy)
    abs_tree = os.path.abspath(tree_file)
    abs_m7_out = os.path.abspath(os.path.join(paml_dir, "m7.out"))
    abs_m8_out = os.path.abspath(os.path.join(paml_dir, "m8.out"))

    # Run M7
    m7_ctl = os.path.join(paml_dir, "m7.ctl")
    write_codeml_ctl(m7_ctl, abs_aln, abs_tree, abs_m7_out, 7)

    # Run M8
    m8_ctl = os.path.join(paml_dir, "m8.ctl")
    write_codeml_ctl(m8_ctl, abs_aln, abs_tree, abs_m8_out, 8)

    curr_dir = os.getcwd()
    try:
        os.chdir(paml_dir)
        logging.info("Running PAML M7")
        run_command("codeml m7.ctl")
        logging.info("Running PAML M8")
        run_command("codeml m8.ctl")
    finally:
        os.chdir(curr_dir)

    return abs_m7_out, abs_m8_out

def parse_codeml_out(out_file):
    lnL = None
    with open(out_file, 'r') as f:
        for line in f:
            if line.startswith("lnL(ntime:"):
                # Example: lnL(ntime: 19  np: 21):  -2566.862663      +0.000000
                match = re.search(r"np:\s*\d+\):\s*([-+0-9.]+)", line)
                if match:
                    lnL = float(match.group(1))
    return lnL

def extract_beb_sites(m8_out):
    sites = []
    in_beb = False
    with open(m8_out, 'r') as f:
        lines = f.readlines()
        for i, line in enumerate(lines):
            if "Bayes Empirical Bayes (BEB) analysis" in line:
                in_beb = True
                continue
            if in_beb and "The grid" in line:
                break

            if in_beb:
                # Look for lines with probabilities like "*", "**"
                # Example:    45 M 0.963*  3.435 +- 1.234
                if "*" in line:
                    parts = line.strip().split()
                    if len(parts) >= 3:
                        site_pos = parts[0]
                        aa = parts[1]
                        prob_str = parts[2]
                        if prob_str.endswith("**"):
                            prob = float(prob_str[:-2])
                            sig = "** (>99%)"
                        elif prob_str.endswith("*"):
                            prob = float(prob_str[:-1])
                            sig = "* (>95%)"
                        else:
                            continue

                        sites.append({'site': site_pos, 'aa': aa, 'prob': prob, 'significance': sig})
    return sites

def perform_lrt(m7_out, m8_out, prefix, out_dir):
    lnL_m7 = parse_codeml_out(m7_out)
    lnL_m8 = parse_codeml_out(m8_out)

    report_file = os.path.join(out_dir, "reports", f"{prefix}_report.txt")

    with open(report_file, 'w') as f:
        f.write(f"=== Positive Selection Analysis Report for {prefix} ===\n\n")

        if lnL_m7 is None or lnL_m8 is None:
            msg = "Error: Could not extract likelihood values from codeml outputs."
            f.write(msg + "\n")
            logging.error(msg)
            return

        # df = 2 for M7 vs M8
        lrt_stat = 2 * (lnL_m8 - lnL_m7)
        # Using chi-square distribution with df=2
        p_val = chi2.sf(lrt_stat, df=2)

        f.write(f"Log Likelihood M7 (Null Model)       : {lnL_m7:.4f}\n")
        f.write(f"Log Likelihood M8 (Alternative Model): {lnL_m8:.4f}\n")
        f.write(f"LRT Statistic (2 * delta lnL)        : {lrt_stat:.4f}\n")
        f.write(f"Degrees of Freedom                   : 2\n")
        f.write(f"P-value                              : {p_val:.4e}\n\n")

        if p_val < 0.05:
            f.write("Result: SIGNIFICANT EVIDENCE of positive selection (p < 0.05).\n\n")
            beb_sites = extract_beb_sites(m8_out)
            if beb_sites:
                f.write("Positively Selected Sites (BEB > 95%):\n")
                f.write(f"{'Site':<10}{'Amino Acid':<15}{'Probability':<15}{'Significance'}\n")
                f.write("-" * 55 + "\n")
                for site in beb_sites:
                    f.write(f"{site['site']:<10}{site['aa']:<15}{site['prob']:<15}{site['significance']}\n")
            else:
                f.write("No specific sites reached >95% BEB significance despite significant LRT.\n")
        else:
            f.write("Result: No significant evidence of positive selection (p >= 0.05).\n")

    logging.info(f"Report generated: {report_file}")

def process_gene(fasta_file, out_dir):
    prefix = os.path.splitext(os.path.basename(fasta_file))[0]
    logging.info(f"Processing gene family: {prefix}")

    try:
        # 1. Validation and Translation
        clean_nuc, prot_fasta = validate_and_translate(fasta_file, prefix, out_dir)

        # 2. Protein Alignment
        prot_aln = align_proteins(prot_fasta, prefix, out_dir)

        # 3. Codon Alignment
        codon_aln, codon_phy = back_translate(clean_nuc, prot_aln, prefix, out_dir)

        # 4. Phylogeny
        tree_file = build_tree(codon_aln, prefix, out_dir)

        # 5. Run PAML
        m7_out, m8_out = run_paml(codon_phy, tree_file, prefix, out_dir)

        # 6. LRT and Reporting
        perform_lrt(m7_out, m8_out, prefix, out_dir)

    except Exception as e:
        logging.error(f"Failed processing {fasta_file}: {e}")

def main():
    parser = argparse.ArgumentParser(description="Automated Pipeline for Positive Selection Analysis using PAML (codeml).")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-i", "--input", help="Single input FASTA file containing homologous CDS sequences.")
    group.add_argument("-d", "--dir", help="Directory containing multiple FASTA files for batch processing.")
    parser.add_argument("-o", "--output", required=True, help="Output directory where all results will be saved.")

    args = parser.parse_args()

    # Create directory structure
    for sub in ["alignments", "trees", "paml_logs", "reports"]:
        os.makedirs(os.path.join(args.output, sub), exist_ok=True)

    if args.input:
        process_gene(args.input, args.output)
    elif args.dir:
        for f in os.listdir(args.dir):
            if f.endswith(".fasta") or f.endswith(".fa") or f.endswith(".fas"):
                fasta_path = os.path.join(args.dir, f)
                process_gene(fasta_path, args.output)

if __name__ == "__main__":
    main()
