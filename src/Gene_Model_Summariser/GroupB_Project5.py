# Project 5
import pandas as pd
import numpy as np
import matplotlib
import os
import gffutils
import seaborn as sns
import Bio.SeqIO as SeqIO
import Bio.Seq as Seq
def main(gff_file, fasta_file=None):
    # temporary placeholder for the main functionality
    print(f"Processing GFF file: {gff_file}")
    if fasta_file:
        print(f"Using reference FASTA file: {fasta_file}")
    else:
        print("No reference FASTA file provided.")
    # Add further processing logic here