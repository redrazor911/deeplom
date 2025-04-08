import os
import glob
import time
import re

def read_virus_taxids(virus_file_path):
    virus_taxids = set()
    try:
        with open(virus_file_path, 'r', encoding='utf-8') as virus_file:
            for line in virus_file:
                line = line.strip()
                if line:
                    virus_taxids.add(line)
        print(f"Прочитано {len(virus_taxids)} вирусных TaxID из {virus_file_path}")
        return virus_taxids
    except FileNotFoundError:
        print(f"Ошибка: Файл не найден по пути: {virus_file_path}")
        exit()
    except UnicodeDecodeError:
        print(f"Ошибка: Проблемы с кодировкой файла {virus_file_path}. Попробуйте указать другую кодировку.")
        exit()
    except Exception as e:
        print(f"Произошла ошибка при чтении {virus_file_path}: {e}")
        exit()


def process_line(line, virus_taxid_set):
    if not line:
        return None
    try:
        protein_id, taxid = line.strip().split('\t')
        if taxid in virus_taxid_set:
            return line
    except ValueError:
        print(f"Skipping malformed line: {line.strip()}")
    return None


def extract_taxid_range_from_filename(filename):
    match = re.search(r'(\d+)-(\d+)', filename)
    if match:
        return int(match.group(1)), int(match.group(2))
    return None, None


def process_files_serial(prot_file, virus_file_path, output_file_path):
    start_time = time.time()
    total_lines = 0
    total_filtered_lines = 0

    # Find all protein files
    protein_files = glob.glob(os.path.join(prot_file, "*.txt"))
    if not protein_files:
        print(f"No protein files found in {prot_file}.")
        return

    try:
        with open(output_file_path, 'w', encoding='utf-8') as filtered_file, \
             open(virus_file_path, 'r', encoding='utf-8') as virus_file:

            for virus_taxid in virus_file:
                virus_taxid = virus_taxid.strip()
                if not virus_taxid:
                    continue # Skip empty lines

                # Find relevant protein file
                relevant_file = None
                for filename in protein_files:
                    min_taxid, max_taxid = extract_taxid_range_from_filename(filename)
                    if min_taxid is not None and max_taxid is not None:
                        if min_taxid <= int(virus_taxid) <= max_taxid:
                            relevant_file = filename
                            break  # Found a relevant file
                
                if not relevant_file:
                    print(f"No file found for virus taxid: {virus_taxid}")
                    continue

                with open(relevant_file, 'r', encoding='utf-8') as input_file:
                    for line in input_file:
                        total_lines += 1
                        if virus_taxid in line: #Упрощенная фильтрация (необязательно)
                            filtered_line = process_line(line, {virus_taxid})
                            if filtered_line:
                                filtered_file.write(filtered_line)
                                total_filtered_lines += 1


    except FileNotFoundError:
        print(f"Error: File not found")
    except UnicodeDecodeError:
        print(f"Error: UnicodeDecodeError. Please check the encoding of your files.")
    except Exception as e:
        print(f"An error occurred: {e}")
    finally:
        end_time = time.time()
        print(f"Total lines processed: {total_lines}")
        print(f"Total filtered lines written to file: {total_filtered_lines}")
        print(f"Execution time: {end_time - start_time:.2f} seconds")
        print(f"Results written to file: {output_file_path}")

def add_arguments(parser):
    parser.add_argument("prot_file", help="Path to the directory containing protein data files.")
    parser.add_argument("virus_file", help="Path to the file containing virus taxids.")
    parser.add_argument("output_file", help="Path to the output file.") # Output file is now a positional argument, no need for -o.

def execute(args):
    virus_taxids = read_virus_taxids(args.virus_file)
    if virus_taxids is None:
        print("Failed to read virus taxids. Exiting.", file=sys.stderr)
        return
    process_files_serial(args.prot_file, virus_taxids, args.output_file)

