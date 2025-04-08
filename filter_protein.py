import re
import time
import sys
from Bio import SeqIO
import multiprocessing
import argparse
import os
import random  


def extract_protein_ids(fasta_file, id_output_file):
    start_time = time.time()
    try:
        with open(fasta_file) as f:
            protein_ids = {p.id.split('|')[1] for p in SeqIO.parse(f, "fasta") if len(p.id.split('|')) > 2}
            with open(id_output_file, 'w') as id_file:
                for pid in protein_ids:
                    id_file.write(f"{pid}\n")
            print(f"Read {len(protein_ids)} Protein IDs from {fasta_file} in {time.time() - start_time:.2f}s. Saved to {id_output_file}")
            return protein_ids
    except Exception as e:
        print(f"Error reading {fasta_file}: {e}")
        sys.exit(1)


def process_chunk(protein_ids, chunk, output_file, found_ids):
    local_found_ids = []
    try:
        for line in chunk:
            if line is None:
                continue
            (prot_id, tax_id) = line.split('\t')
            if '.' in prot_id:
                prot_id = prot_id.split('.')[0]
            if prot_id in protein_ids:
                local_found_ids.append((prot_id, tax_id))
    except Exception as e:
        print(f"Error processing chunk: {e}")
        return

    with open(output_file, 'a', encoding='utf-8') as outfile:
        for prot_id, tax_id in local_found_ids:
            outfile.write(f"{prot_id}\t{tax_id}\n")

    with found_ids.get_lock():
        found_ids.value.update(local_found_ids)


def determine_optimal_chunk_size(prot_file, target_memory_mb=250):
    file_size_bytes = os.path.getsize(prot_file)

    if file_size_bytes == 0:
        print("Warning: Protein file is empty.")
        return 0
    avg_line_size_bytes = estimate_average_line_size(prot_file, num_lines_to_sample=1000)

    if avg_line_size_bytes == 0:
        print("Warning: Average line size could not be determined. Returning default chunk_size.")
        return 1000

    chunk_size = target_memory_mb * 1024 * 1024 
    chunk_size = min(chunk_size, file_size_bytes)
    if chunk_size < avg_line_size_bytes:
        chunk_size = avg_line_size_bytes 


    max_chunk_size = 1000000000

    chunk_size = min(chunk_size, max_chunk_size)
    return int(chunk_size)


def estimate_average_line_size(prot_file, num_lines_to_sample=1000):
    total_sampled_size = 0
    num_sampled = 0
    try:
        with open(prot_file, 'r', encoding='utf-8') as f:
            for i, line in enumerate(f):
                total_sampled_size += len(line.encode('utf-8'))
                num_sampled += 1
                if num_sampled >= num_lines_to_sample:
                    break
    except Exception as e:
        print(f"Error estimating average line size: {e}")
        return 100

    if num_sampled == 0:
        return 100
    return float(total_sampled_size) / num_sampled



def process_files_parallel(fasta_file, prot_file, output_file, num_processes, id_output_file="fasta_ids.txt", max_lines=None):
    start_time = time.time()
    protein_ids = extract_protein_ids(fasta_file, id_output_file)
    manager = multiprocessing.Manager()
    found_ids = manager.Namespace()
    found_ids.value = set()
    found_ids.lock = manager.Lock()
    proc_lines = 0

    optimal_chunk_size = determine_optimal_chunk_size(prot_file, num_processes)

    if optimal_chunk_size is None:
        print("File is too small to chunk effectively.  Using a minimum chunk size.")
        chunk_size = 1000 
    else:
        print(f"Optimal chunk size: {optimal_chunk_size}")
        chunk_size = optimal_chunk_size

    with open(prot_file, 'r', encoding='utf-8') as infile:
        pool = multiprocessing.Pool(processes=num_processes)

        results = []
        while True:
            chunk = [next(infile, None) for _ in range(chunk_size)]
            if not any(chunk):
                break
            results.append(pool.apply_async(process_chunk, (protein_ids, chunk, output_file, found_ids)))
            proc_lines += len(chunk)
            if max_lines and proc_lines >= max_lines:
                break
        pool.close()
        pool.join()

    missing = set(protein_ids) - found_ids.value

    total_time = time.time() - start_time
    print(f"Results written to {output_file} in {total_time:.2f}s using {num_processes} processes.")


def benchmark(fasta_file, prot_file, output_file, num_processes_list, id_output_file="fasta_ids.txt", max_lines=None):
    print("Running benchmark...")
    results = []
    for num_processes in num_processes_list:
        start_time = time.time()
        process_files_parallel(fasta_file, prot_file, output_file, num_processes, id_output_file, max_lines)
        execution_time = time.time() - start_time
        results.append((num_processes, execution_time))
        print(f"Benchmark with {num_processes} processes completed in {execution_time:.2f} seconds.")

    print("Benchmark results:")
    for num_processes, execution_time in results:
        print(f"Processes: {num_processes}, Time: {execution_time:.2f}s")

def add_arguments(parser):
    parser.add_argument("-f", "--fasta_file", dest="fasta_file", default="uniprot_sprot.fasta",
                        help="Path to the FASTA file.")
    parser.add_argument("-p", "--prot_file", dest="prot_file", default="prot.accession2taxid.FULL",
                        help="Path to the protein file.")
    parser.add_argument("-o", "--output_file", dest="output_file", default="protein10_ids.txt",
                        help="Path to the output file.")
    parser.add_argument("-i", "--fasta_ids_file", dest="fasta_ids_file", default="fasta_ids.txt",
                        help="Path to the FASTA IDs file. If not provided will be created.")
    parser.add_argument("-proc", "--num_processes", dest="num_processes", type=int, default=4,
                        help="Number of processes to use (default: 4).")
    parser.add_argument("-b", "--benchmark", action="store_true", dest="b", help="Enable benchmark mode.")
    parser.add_argument("-m", "--max_lines", dest="max_lines", type=int, default=None,
                        help="Maximum number of lines to process from fasta file.")
    parser.add_argument("-mb", "--memory", dest="target_memory_mb", type=int, default=500,
                        help="The target memory usage (MB).  Defaults to 500MB")
    parser.add_argument("-b_proc", "--benchmark_processes", dest="benchmark_processes", nargs='+', type=int, default=[1, 2, 4, 8], help="List of process counts to use for benchmarking.  Example: --benchmark_processes 1 4 8")


def execute(args):
    if args.b:
        benchmark(args.fasta_file, args.prot_file, args.output_file, args.benchmark_processes, args.fasta_ids_file, args.max_lines)
    else:
        process_files_parallel(args.fasta_file, args.prot_file, args.output_file, args.num_processes, args.fasta_ids_file, args.max_lines)
