import os

import os

def split_protein_file_by_taxid_range(input_file, output_dir, range_size=100000, chunk_size=1000000):
    """
    Reads a large protein file and splits it into multiple files based on tax_id ranges,
    handling files that do not fit in memory by processing in chunks.  Includes more robust error checking.

    Args:
        input_file (str): Path to the input protein file (protein10_ids.txt).
        output_dir (str): Path to the output directory where the split files will be created.
        range_size (int): The size of each tax_id range (default: 100000).
        chunk_size (int): The number of lines to read from the input file at once.
    """

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    file_handles = {}  # Store open file handles for each range

    try:
        with open(input_file, 'r', encoding='utf-8') as infile:
            while True:
                chunk = [next(infile, None) for _ in range(chunk_size)]
                if not any(chunk):
                    break  # End of file

                for line in chunk:
                    if line is None:
                        continue
                    line = line.strip()
                    if not line:  # Skip empty lines after stripping whitespace
                        continue

                    parts = line.split('\t')
                    if len(parts) != 2:
                        print(f"Skipping malformed line (incorrect number of parts): {line}")
                        continue

                    try:
                        prot_id, tax_id = parts
                        tax_id = int(tax_id)  # Convert tax_id to integer
                    except ValueError as e:
                        print(f"Skipping malformed line (ValueError): {line}. Error: {e}")
                        continue


                    # Determine the range for the current tax_id
                    range_start = (tax_id // range_size) * range_size
                    range_end = range_start + range_size - 1
                    range_name = f"{range_start}-{range_end}"

                    # Get or create the output file handle for this range
                    if range_name not in file_handles:
                        output_file = os.path.join(output_dir, f"{range_name}.txt")
                        file_handles[range_name] = open(output_file, 'w', encoding='utf-8')

                    # Write the line to the appropriate file
                    file_handles[range_name].write(line + '\n')

    except Exception as e:
        print(f"An error occurred: {e}")

    finally:
        # Close all open file handles
        for handle in file_handles.values():
            if handle:
                handle.close()


def add_arguments(parser):
    """Adds command line arguments to the parser."""
    parser.add_argument("input_file", help="Path to the input protein file.")
    parser.add_argument("output_dir", help="Path to the output directory.")
    parser.add_argument("--range_size", type=int, default=100000, help="Size of the taxid range for each output file (default: 100000).")
    parser.add_argument("--chunk_size", type=int, default=1000000, help="Size of the chunk to read from the input file (default: 1000000).")

def execute(args):
    """Executes the protein file splitting process."""
    split_protein_file_by_taxid_range(args.input_file, args.output_dir, args.range_size, args.chunk_size)

