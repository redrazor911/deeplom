import os

def split_protein_file_by_taxid_range(input_file, output_dir, range_size=100000, chunk_size=1000000):

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    file_handles = {}  

    try:
        with open(input_file, 'r', encoding='utf-8') as infile:
            while True:
                chunk = [next(infile, None) for _ in range(chunk_size)]
                if not any(chunk):
                    break  

                for line in chunk:
                    if line is None:
                        continue
                    line = line.strip()
                    if not line:  
                        continue

                    parts = line.split('\t')
                    if len(parts) != 2:
                        print(f"Skipping malformed line (incorrect number of parts): {line}")
                        continue

                    try:
                        prot_id, tax_id = parts
                        tax_id = int(tax_id) 
                    except ValueError as e:
                        print(f"Skipping malformed line (ValueError): {line}. Error: {e}")
                        continue


              
                    range_start = (tax_id // range_size) * range_size
                    range_end = range_start + range_size - 1
                    range_name = f"{range_start}-{range_end}"

              
                    if range_name not in file_handles:
                        output_file = os.path.join(output_dir, f"{range_name}.txt")
                        file_handles[range_name] = open(output_file, 'w', encoding='utf-8')

                    
                    file_handles[range_name].write(line + '\n')

    except Exception as e:
        print(f"An error occurred: {e}")

    finally:
       
        for handle in file_handles.values():
            if handle:
                handle.close()


def add_arguments(parser):
 
    parser.add_argument("input_file", help="Path to the input protein file.")
    parser.add_argument("output_dir", help="Path to the output directory.")
    parser.add_argument("--range_size", type=int, default=100000, help="Size of the taxid range for each output file (default: 100000).")
    parser.add_argument("--chunk_size", type=int, default=1000000, help="Size of the chunk to read from the input file (default: 1000000).")

def execute(args):

    split_protein_file_by_taxid_range(args.input_file, args.output_dir, args.range_size, args.chunk_size)

