import argparse
import taxid_parser
import filter_protein
import index_protein
import filter_stage_2
def main():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command")
    
    parser_taxid_parser = subparsers.add_parser("taxid_parser")
    taxid_parser.add_arguments(parser_taxid_parser)
    
    parser_filter_protein = subparsers.add_parser("filter_protein")
    filter_protein.add_arguments(parser_filter_protein)
    
    parser_filter_protein = subparsers.add_parser("index_protein")
    index_protein.add_arguments(parser_index_protein)
    
    parser_filter_stage_2 = subparsers.add_parser("filter_stage_2")
    filter_stage_2.add_arguments(parser_filter_stage_2)
    
    args = parser.parse_args()
    if args.command == "taxid_parser":
        taxid_parser.execute(args)  
    elif args.command == "filter_protein":
        filter_protein.execute(args)
    elif args.command == "index_protein":
        index_protein.execute(args)  
    elif args.command == "filter_stage_2":
        filter_stage_2.execute(args)  
    else:
        parser.print_help()  
if __name__ == "__main__":
    main()
