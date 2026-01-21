import argparse


def app():
    parser = argparse.ArgumentParser(description="A simple CLI tool.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-g','--gff', required=True, help='Path to GFF file')
    parser.add_argument('-f','--fasta', help='Path to optional reference FASTA file for sequence-derived metrics')
    args = parser.parse_args()
    
    print('hello world')
    
    return 0