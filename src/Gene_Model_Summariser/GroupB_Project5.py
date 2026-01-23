# Project 5
from ast import List
import math
from venv import logger
from xml.parsers.expat import errors
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
        results = QC_flags(db, fasta).process_all_sequences()

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
