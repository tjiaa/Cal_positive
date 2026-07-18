import os
import sys
import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def get_strain_name(record):
    """
    尝试从 record 的 features (source) 中提取菌株名或物种名。
    优先级: strain > organism > record.id
    """
    for feature in record.features:
        if feature.type == "source":
            if "strain" in feature.qualifiers:
                return feature.qualifiers["strain"][0].replace(" ", "_")
            elif "organism" in feature.qualifiers:
                return feature.qualifiers["organism"][0].replace(" ", "_")
    return record.id.replace(" ", "_")

def extract_genes(target_genes, gbk_dir, output_dir):
    """
    遍历目录下的所有 gbk 文件，提取指定的基因 CDS 序列，
    并为每个基因生成一个合并的 FASTA 文件。
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Dictionary to hold lists of SeqRecords for each target gene
    gene_records = {gene: [] for gene in target_genes}

    gbk_files = [f for f in os.listdir(gbk_dir) if f.endswith(('.gbk', '.gb', '.genbank'))]
    if not gbk_files:
        logging.warning(f"No GenBank files found in directory: {gbk_dir}")
        return

    for gbk_file in gbk_files:
        gbk_path = os.path.join(gbk_dir, gbk_file)
        logging.info(f"Processing {gbk_file}...")

        try:
            # Parse the GBK file (can contain multiple records, e.g., chromosomes/plasmids)
            for record in SeqIO.parse(gbk_path, "genbank"):
                strain_name = get_strain_name(record)

                for feature in record.features:
                    # We strictly look for CDS features to get coding sequences
                    if feature.type == "CDS":
                        gene_qualifier = feature.qualifiers.get("gene", [""])[0]

                        if gene_qualifier in target_genes:
                            # 自动处理正反链，提取精确的核苷酸序列
                            nuc_seq = feature.extract(record.seq)

                            # 序列ID格式：物种名/菌株名_基因名
                            seq_id = f"{strain_name}_{gene_qualifier}"
                            # 替换掉可能导致问题的特殊字符
                            seq_id = seq_id.replace(" ", "_").replace("/", "_").replace("(", "").replace(")", "")

                            seq_record = SeqRecord(
                                nuc_seq,
                                id=seq_id,
                                description=f"Extracted from {gbk_file}"
                            )

                            gene_records[gene_qualifier].append(seq_record)
        except Exception as e:
            logging.error(f"Error processing file {gbk_file}: {e}")

    # Write output fasta files
    for gene, records in gene_records.items():
        if records:
            out_fasta = os.path.join(output_dir, f"{gene}.fasta")
            SeqIO.write(records, out_fasta, "fasta")
            logging.info(f"Successfully extracted {len(records)} sequences for gene '{gene}' into {out_fasta}")
        else:
            logging.warning(f"No sequences found for target gene: {gene}")

def main():
    parser = argparse.ArgumentParser(description="Extract specific gene CDS sequences from multiple GenBank files.")
    parser.add_argument("genes", nargs="+", help="List of gene names to extract (e.g., thrA talB)")
    parser.add_argument("-d", "--dir", default=".", help="Directory containing .gbk files (default: current directory)")
    parser.add_argument("-o", "--output", default=".", help="Output directory for FASTA files (default: current directory)")

    args = parser.parse_args()

    extract_genes(args.genes, args.dir, args.output)

if __name__ == "__main__":
    main()
