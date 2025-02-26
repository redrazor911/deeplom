import re, time, sys, argparse, threading
from Bio import SeqIO
from queue import Queue

def process_files(fasta_file, prot_file, output_file, id_output_file="fasta_ids.txt", num_threads=8, benchmark=False):
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

   
    lines_queue = Queue()
    try:
        with open(prot_file, 'r', encoding='utf-8') as infile:
            for line in infile:
                lines_queue.put(line.strip())
    except Exception as e:
        print(f"Error reading {prot_file}: {e}")
        sys.exit(1)

    
    found_ids = set()
    lock = threading.Lock()  


    def worker():
        while not lines_queue.empty():
            line = lines_queue.get()
            parts = line.split('\t')
            if len(parts) == 2:
                sub = parts[0].split('.')
                if sub[0].strip() in protein_ids:
                    with lock:  
                        found_ids.add(sub[0])
                        try:
                            with open(output_file, 'a', encoding='utf-8') as outfile:  # Open in append mode
                                outfile.write(f"{sub[0]}\t{parts[1]}\n")
                        except Exception as e:
                            print(f"Error writing to output file: {e}")
            lines_queue.task_done()  

   
    threads = []
    thread_start_time = time.time()
    for _ in range(num_threads):
        t = threading.Thread(target=worker)
        threads.append(t)
        t.start()


    lines_queue.join()


    for t in threads:
        t.join()

    thread_end_time = time.time()
    processing_time = thread_end_time - thread_start_time
    print(f"Processed {lines_queue.qsize()} lines with {num_threads} threads in {processing_time:.2f} seconds.")

    missing = protein_ids - found_ids
    if missing:
        print("Warning: Some Protein IDs not found:")
        for mid in missing: print(f"- {mid}")
    else:
        print("All Protein IDs from FASTA found.")

    total_time = time.time() - start_time
    print(f"Results written to {output_file} in {total_time:.2f}s")

    # Benchmark mode
    if benchmark:
        benchmark_results = {}
        for n_threads in [1, 2, 4, 8, 16, 32]: 
            start_bench_time = time.time()
            found_ids_bench = set()  
            lines_queue_bench = Queue() 
            try:
                with open(prot_file, 'r', encoding='utf-8') as infile:
                    for line in infile:
                        lines_queue_bench.put(line.strip())
            except Exception as e:
                print(f"Error reading {prot_file} in benchmark: {e}")
                sys.exit(1)


            def bench_worker():
                while not lines_queue_bench.empty():
                    line = lines_queue_bench.get()
                    parts = line.split('\t')
                    if len(parts) == 2:
                        sub = parts[0].split('.')
                        if sub[0].strip() in protein_ids: 
                            with lock:
                                found_ids_bench.add(sub[0])

                    lines_queue_bench.task_done()

            threads_bench = []

            for _ in range(n_threads):
                t = threading.Thread(target=bench_worker)
                threads_bench.append(t)
                t.start()

            lines_queue_bench.join()

            for t in threads_bench:
                t.join()

            end_bench_time = time.time()
            benchmark_results[n_threads] = end_bench_time - start_bench_time
            print(f"Benchmark with {n_threads} threads completed in {benchmark_results[n_threads]:.2f} seconds")

        print("\nBenchmark Results:")
        for threads, time_taken in benchmark_results.items():
            print(f"{threads} threads: {time_taken:.2f} seconds")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process FASTA and protein files with multi-threading.")
    parser.add_argument("fasta_file", help="Path to the FASTA file")
    parser.add_argument("prot_file", help="Path to the prot.accession2taxid.FULL file")
    parser.add_argument("-o", "--output", dest="output_file", default="protein_ids2.txt", help="Path to the output file")
    parser.add_argument("-i", "--ids", dest="fasta_ids_file", default="fasta_ids.txt", help="Path to save FASTA IDs")
    parser.add_argument("-t", "--threads", dest="num_threads", type=int, default=8, help="Number of threads to use (default: 8)")
    parser.add_argument("-b", "--benchmark", action="store_true", help="Run benchmark mode to evaluate performance with different thread counts")

    args = parser.parse_args()

    process_files(args.fasta_file, args.prot_file, args.output_file, args.fasta_ids_file, args.num_threads, args.benchmark)
