import argparse
from skbio.io import read
import pandas as pd
import sys

def get_virus_taxids(names_file, nodes_file):
    try:
        names_df = read(names_file, format='taxdump', into=pd.DataFrame, scheme="names")
        nodes_df = read(nodes_file, format='taxdump', into=pd.DataFrame, scheme="nodes_slim")
        names_df = names_df[names_df['name_class'] == 'scientific name']
        nodes_df['lineage_names'] = [names_df.loc[get_lineage(taxid, nodes_df), 'name_txt'].tolist() for taxid in nodes_df.index]
        return sorted(nodes_df[nodes_df['lineage_names'].apply(lambda x: 'Viruses' in x)].index.tolist())
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        return []

def get_lineage(taxid, nodes_df):
    lineage = []
    current_id = taxid
    while current_id != 1 and current_id in nodes_df.index:
        lineage.append(current_id)
        current_id = nodes_df.loc[current_id, 'parent_tax_id']
    lineage.append(1)
    return lineage

def execute(args):
    taxids = get_virus_taxids(args.names_file, args.nodes_file)
    if taxids:
        with open(args.output, "w") as f:
            f.write("\n".join(map(str, taxids)) + "\n")
        print(f"Virus TaxIDs written to {args.output}")
    else:
        print("No virus TaxIDs found.", file=sys.stderr)

def add_arguments(parser):
    parser.add_argument("names_file", help="Path to names.dmp file")
    parser.add_argument("nodes_file", help="Path to nodes.dmp file")
    parser.add_argument("-o", "--output", help="Output file for virus TaxIDs", default="virus_taxids.txt")
