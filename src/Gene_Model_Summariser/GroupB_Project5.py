# Project 5
import pandas as pd
import numpy as np
import matplotlib
import os
import gffutils
import seaborn as sns
import logging
import Bio.SeqIO as SeqIO
import Bio.Seq as Seq
import sqlite3
from typing import Optional
from .fasta_validator import FastaChecker

# This is the main function for the Gene Model Summariser. 
def main(gff_file: str, fasta_file: Optional[str] = None) -> None:
    db = load_gff_database(gff_file)
    logger = setup_logger("gene_model_summariser.log")
    if fasta_file:
        fasta_checker = FastaChecker(fasta_file)
        if not fasta_checker.validate_fasta():
            logger.error("Invalid FASTA file provided. Exiting.")
            raise SystemExit(1)
        else:
            logger.info("FASTA file validated successfully.")
        fasta = fasta_checker.fasta_parse()

def setup_logger(log_file: str) -> logging.Logger:
    logger = logging.getLogger("GroupB_logger")
    logger.setLevel(logging.INFO)
    #prevent duplicates if run script multiple times
    if not logger.handlers:
        file = logging.FileHandler(log_file)
        file.setFormatter(logging.Formatter("%(levelname)s - %(asctime)s - %(message)s"))
        logger.addHandler(file)
    return logger


def load_gff_database(gff_file: str) -> gffutils.FeatureDB: # Create or connect to GFF database.
    db_path = gff_file.replace('.gff', '.db') # replace .gff with .db for database file name
    if not os.path.isfile(db_path): # if the gff.db file does not exist, create it
        try:
            db = gffutils.create_db(gff_file, dbfn=db_path, force=True, keep_order=True)
        except (sqlite3.OperationalError, ValueError):
            raise SystemExit(1)
    else: # if it does exist, connect to it
        try:
            db = gffutils.FeatureDB(db_path, keep_order=True) # connect to existing database
        except ValueError:
            raise SystemExit(1)
    return db # return the database object as db

class QC_flags:
    # Class to generate QC flags for gene models from parser data
    def __init__(self, db: gffutils.FeatureDB, fasta_file: Optional[str] = None) -> None:
        self.db = db
        self.fasta_file = fasta_file
        self.fasta = None
        if fasta_file:
            self.fasta = SeqIO.to_dict(SeqIO.parse(fasta_file, 'fasta'))
    
    def gc_content(self, sequence: str) -> float:
        # Function to calculate GC content of a given sequence
        if not sequence:
            return 0
        sequence = sequence.upper().rstrip()
        gc_count = sequence.count('G') + sequence.count('C')
        return gc_count / len(sequence) * 100
    
    def sequence_length(self, sequence: str) -> int:
        # Function to calculate the length of a given sequence
        if not sequence:
            return 0
        return len(sequence.rstrip())
    
    def N_content(self, sequence: str) -> tuple[int, float]:
        # Function to calculate the number of 'N' bases in a given sequence
        if not sequence:
            return 0, 0.0
        N_count = sequence.upper().count('N')
        try:
            N_cont_percent = (N_count / len(sequence)) * 100
        except ZeroDivisionError:
            N_cont_percent = 0.0
        return N_count, N_cont_percent

# Project 5: Gene Model Summariser
# Group B

#############################################################
#GFF PARSER TESTER
#############################################################



'''
Part 1 - parse the GFF

1.	Skip the hashtag lines 
2.	Validate it has 9 columns (if too little/too many, tell the user and continue to avoid crashing)
3.	Validate it is a gff file 
4.	Check every column is in the correct place and the file is robust (ie it follows the correct 9 column format) 
5.	Check every column has the correct data type and handle NA values/. Values correctly for each column - if any issues record a QC issue for flags and run.json. make sure to include a continue so the row doesn’t break 
6.	Parse the confirmed clean columns
7.	Parse attributes into a dictionary of contents so we can grab each section individually

2. Builds dicts so we can map the relationships 

1.	Map transcript_id to the gene_id (from transcript line)
2.	Need to parse the attributes line safely to make sure we link the geneID to the transcript_id
3.	Exon coutns - add a counter for the type column where it is exons - need to make sure the line is parsed correctly. Use parent=transcript_id here to map the exon back to the transcript to avoid overlaps 
4.	Has cds - easy just set true if type=cds and else false 

3. Grab the transcript and gene_ids

Gene lines set as type=gene with the gene_Id. 
Parse transcript lines to grab mRNA/transcript with the ID 
Primary transcript becomes a combination (ask hans/john any other suggestions) 
Exon/cds features attach via the parent class to transcripts

4. Calculations for the output file and QC - have we found any errors when calculating 

gene_id: from Parent on the transcript line
n_exons: from exon_count function 
has_cds: from has_cds function 
flags: grab from the flags_by_trancript function and resolve them 

5. Design the outputs in main()
.json? - speak with john/hans about the best way to do this
tsv should be easy enough inc the qc_flags.gff3
###########################################################
Testing logic?:
part 1 - parsing testing 
make sure the parser skips hashtag lines 
check strand is only +, - or . 
check start and end are integers and start < end
check seqid, source, type are strings
check score is float or NA
check phase is 0,1,2 or NA  
parsing the file should be easy enough using our old dataset
need to make sure each column in the gff is where it shoud be and correct data type
check each row is the correct length (9 columns)
check each data point is correct and normalise NA values
warn the user if not and store data for QC flags 
check for any missing data in the relationships
make sure each transcript links to a gene 
Parent=tx1;gene_id=geneA;ID=tx1 - make sure each row looks like this and warn user if it doesnt 
make sure to remove whitespace to avoid empty strings/crashing 
check for duplicates?
check there are no strange characters in the IDs
might want to check all the data is the same casing for consisrtency (lowercase everything?)

part 2: transcript --> gene testing 
check to make sure multiple transcripts link back to the same gene_id
could we check logically that each transcript has exons/cds?
check to make sure we dont have duplicate transcript ID appears twice with different gene IDs

part 3: Exon counting tests
check we are correctly counting exons/cds correctly 
Exon with multiple parents (Parent=tx1,tx2) - make sure both transcripts get the exon counted  
Two exons for the same transcript overlap in coordinates should count as 2 

part 4: has_cds tests
One CDS line with Parent=tx1 - expect has_cds=True test to return True 
no cds remains False 
CDS line with multiple parents - flag both lines to be true 

part 5: QC flag tests
missing gene_id in attributes - expect flag recorded for missing gene_id 

###### Extras on top (L4) ######
Part 6 - integrate FASTA file to get sequence lengths
Pushing to Bioconda
Push to nextflow and underastand how to use in pipeline
Make sure all dependencies are in place
Create environment.yml for conda
Make sure to test conda package locally before pushing
Set up bioconda recipe and submit PR
'''
