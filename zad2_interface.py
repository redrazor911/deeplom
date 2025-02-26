import re, time, sys, argparse
from Bio import SeqIO

def process_files(fasta_file, prot_file, output_file, id_output_file="fasta_ids.txt"):
    start_time = time.time()
    try:
        with open(fasta_file) as f:
            protein_ids = {p.id.split('|')[1] for p in SeqIO.parse(f, "fasta") if len(p.id.split('|')) > 2}
            t = time.time()
            with open(id_output_file, 'w') as id_file:
                for pid in protein_ids:
                    id_file.write(f"{pid}\n")

            print(f"Read {len(protein_ids)} Protein IDs from {fasta_file} in {t-start_time:.2f}s. Saved to {id_output_file}")

    except Exception as e:
        print(f"Error reading {fasta_file}: {e}")
        sys.exit(1)

    found_ids = set()
    try:
        with open(prot_file, 'r', encoding='utf-8') as infile, \
             open(output_file, 'w', encoding='utf-8') as outfile:
            k = 0
            for line in (l.strip() for l in infile):
                k += 1
                if k%100000000==0:
                    print('analyzed 100m strings')
                parts = line.split('\t')
                if len(parts) == 2:
                    sub = parts[0].split('.')
                    if sub[0].strip() in protein_ids:
                        found_ids.add(sub[0])
                        outfile.write(f"{sub[0]}\t{parts[1]}\n")

    except Exception as e:
        print(f"Error reading {prot_file}: {e}")
        sys.exit(1)

    missing = protein_ids - found_ids
    if missing:
        print("Warning: Some Protein IDs not found:")
        for mid in missing: print(f"- {mid}")
    else:
        print("All Protein IDs from FASTA found.")

    total_time = time.time() - start_time
    print(f"Results written to {output_file} in {total_time:.2f}s")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process FASTA and protein files to extract matching IDs.")
    parser.add_argument("fasta_file", help="Path to the FASTA file")
    parser.add_argument("prot_file", help="Path to the prot.accession2taxid.FULL file")
    parser.add_argument("-o", "--output", dest="output_file", default="protein_ids.txt", help="Path to the output file")
    parser.add_argument("-i", "--ids", dest="fasta_ids_file", default="fasta_ids.txt", help="Path to save FASTA IDs")

    args = parser.parse_args()

    process_files(args.fasta_file, args.prot_file, args.output_file, args.fasta_ids_file)
