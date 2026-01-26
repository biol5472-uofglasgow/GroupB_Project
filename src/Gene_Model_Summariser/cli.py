import argparse
from .GroupB_Project5 import main


def app():
    parser = argparse.ArgumentParser(description="A simple CLI tool.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-g','--gff', required=True, help='Path to GFF file')
    parser.add_argument('-f','--fasta', help='Path to optional reference FASTA file for sequence-derived metrics')
    parser.add_argument('-o','--outdir', help='Output directory for results', default='results/run_001/')
    args = parser.parse_args()
    
    # Call the main function from GroupB_Project5.py with the parsed arguments
    main(args.gff, args.fasta, args.outdir)
    
    return 0