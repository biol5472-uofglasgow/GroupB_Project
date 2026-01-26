# Project 5 - Gene Model Summariser
import math
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
from .QC_check import QC_flags
from .gff_parser import GFF_Parser
from .gff_validator import check_db

# This is the main function for the Gene Model Summariser. 
def main(gff_file: str, fasta_file: Optional[str] = None, output_dir: str = ".") -> None:
    db = load_gff_database(gff_file)
    db_check = check_db(db)
    logger = setup_logger("gene_model_summariser.log")
    if db_check:
        tsv_results = GFF_Parser(db).tsv_output()
        if fasta_file:
            fasta_checker = FastaChecker(fasta_file)
            if not fasta_checker.validate_fasta():
                logger.error("Invalid FASTA file provided. Exiting.")
                raise SystemExit(1)
            fasta = fasta_checker.fasta_parse()
            results = QC_flags(db, fasta).gff_QC()
        else:
            results = QC_flags(db).gff_QC()
        output_results(tsv_results, results, output_dir)
    else:
        logger.error("GFF database validation failed. Exiting.")

# Join tsv_results and results on transcript IDs for output in singular results.tsv file. 
def output_results(tsv_data: dict, qc_data: dict, output_dir: str) -> None:
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    combined_data = []
    for transcript_id, tsv_metrics in tsv_data.items():
        qc_flags = qc_data.get(transcript_id, [])
        # Convert QC flags list to comma-separated string
        qc_flags_str = ','.join(qc_flags) if qc_flags else ''
        combined_entry = {**tsv_metrics, 'flags': qc_flags_str}
        combined_data.append(combined_entry)
    
    df = pd.DataFrame(combined_data)
    output_path = os.path.join(output_dir, "results.tsv")
    df.to_csv(output_path, sep='\t', index=False)

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
    gff_file = gff_file.replace('.gff3', '.gff').replace('.gff.gz', '.gff') # normalize file extension
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
