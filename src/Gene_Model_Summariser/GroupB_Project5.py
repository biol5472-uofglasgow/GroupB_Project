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
from pathlib import Path
from .fasta_validator import FastaChecker
from .QC_check import QC_flags
from .gff_parser import GFF_Parser
from .gff_validator import check_db
from .qc_flags_bed import TranscriptWithFlags, write_qc_bed
from .build_gff import build_gff
from .run_json_builder import make_run_json_file
from .run_json_builder import finalise_run_json_file


# This is the main function for the Gene Model Summariser. 
def main(gff_file: str, fasta_file: Optional[str] = None, output_dir: str = ".") -> None:
    """
    Main function for the Gene Model Summariser.
    gff_file: Path to the GFF file.
    fasta_file: Optional path to the FASTA file.
    output_dir: Directory where output files will be saved.
    """
    logger = setup_logger("gene_model_summariser.log") # Setup logger
    out_dir = Path(output_dir) #makes and set the output dir for run.json file 
    out_dir.mkdir(parents=True, exist_ok=True)
    run_path = None #add run_path - used later to check if run.json file made by finalise_run_json_file
    make_run_json_file(gff_file=Path(gff_file),fasta_file=Path(fasta_file) if fasta_file else None,
        output_dir=out_dir, results_filename="results.tsv",html_filename="results.html",run_filename="run.json")
    try:
        db = load_gff_database(gff_file) # Load or create GFF database
    except SystemExit:
        logger.error("Failed to load or create GFF database. Exiting.") # Log error if database loading fails
        raise SystemExit(1)
    try: 
        db_check = check_db(db) # Validate the GFF database
    except SystemExit:
        logger.error("GFF database validation encountered an error. Exiting.") # Log error if validation fails
        raise SystemExit(1)
    
    if db_check: # If database is valid, proceed
        tsv_results = GFF_Parser(db).tsv_output() # Parse GFF and generate TSV results
        if fasta_file: # If a FASTA file is provided, validate and parse it
            fasta_checker = FastaChecker(fasta_file, logger) # Create FastaChecker instance
            if not fasta_checker.validate_fasta(): # If the FASTA file is invalid, log error and exit
                logger.error("Invalid FASTA file provided. Exiting.") # Log error for invalid FASTA
                raise SystemExit(1)
            fasta = fasta_checker.fasta_parse() # Parse the FASTA file
            results = QC_flags(db, fasta).gff_QC() # Generate QC flags using both GFF and FASTA data
        else:
            results = QC_flags(db).gff_QC() # Generate QC flags using only GFF data
        output_results(tsv_results, results, output_dir, gff_file, db) # Output combined results to TSV file
    else:
        logger.error("GFF database validation failed. Exiting.") # Log error if GFF validation fails
        raise SystemExit(1)
    finalise_run_json_file(output_dir=out_dir, run_filename="run.json")

# Join tsv_results and results on transcript IDs for output in singular results.tsv file. 
def output_results(tsv_data: dict, qc_data: dict, output_dir: str, gff_file: str, db) -> None:
    """
    combine tsv_data and qc_data into a single TSV file in output_dir.
    tsv_data: Dictionary containing TSV metrics keyed by transcript IDs.
    qc_data: Dictionary containing QC flags keyed by transcript IDs.
    output_dir: Directory where the results.tsv file will be saved.
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    combined_data = [] # list to hold combined entries
    transcripts_with_flags = []

    gff_path = os.path.join(output_dir, 'qc_flags.gff3')
    with open(gff_path,'w') as gff_out:
        with open(gff_file, 'r') as gff_in:
            for line in gff_in:
                if line.startswith('##'):
                    gff_out.write(line)
    
        for transcript_id, tsv_metrics in tsv_data.items(): # iterate through tsv_data
            qc_flags = qc_data.get(transcript_id, []) # get corresponding QC flags
            qc_flags_str = ','.join(qc_flags) if qc_flags else '' # Convert QC flags list to comma-separated string
            combined_entry = {**tsv_metrics, 'flags': qc_flags_str} # merge dictionaries
            combined_data.append(combined_entry) # add to combined list
            
            # Create BED entry if transcript has flags
            if qc_flags:
                transcript = TranscriptWithFlags(
                    chrom=tsv_metrics.get('chrom'),
                    start=tsv_metrics.get('start'),
                    end=tsv_metrics.get('end'),
                    transcript_id=transcript_id,
                    qc_flags=set(qc_flags),
                    strand=tsv_metrics.get('strand')
                )
                transcripts_with_flags.append(transcript)

                    # Also write to GFF with QC flags
                build_gff(transcript_id, db, qc_flags_str, gff_out)
        
    # Write BED file if there are flagged transcripts
    if transcripts_with_flags:
        bed_path = Path(output_dir) / "qc_flagged.bed"
        write_qc_bed(transcripts_with_flags, bed_path)
    for transcript_id, tsv_metrics in tsv_data.items(): # iterate through tsv_data
        qc_flags = qc_data.get(transcript_id, []) # get corresponding QC flags
        qc_flags_str = ','.join(qc_flags) if qc_flags else '' # Convert QC flags list to comma-separated string
        combined_entry = {**tsv_metrics, 'flags': qc_flags_str} # merge dictionaries
        combined_data.append(combined_entry) # add to combined list
    
    df = pd.DataFrame(combined_data) # create DataFrame from combined data
    output_path = os.path.join(output_dir, "results.tsv") # define output file path
    df.to_csv(output_path, sep='\t', index=False) # save DataFrame to TSV file

def setup_logger(log_file: str) -> logging.Logger:
    """
    setup_logger: Configures and returns a logger that writes to the specified log file.
    log_file: Path to the log file where log messages will be written.
    """
    logger = logging.getLogger("GroupB_logger")
    logger.setLevel(logging.INFO)
    #prevent duplicates if run script multiple times
    if not logger.handlers:
        file = logging.FileHandler(log_file)
        file.setFormatter(logging.Formatter("%(levelname)s - %(asctime)s - %(message)s"))
        logger.addHandler(file)
    return logger


def load_gff_database(gff_file: str) -> gffutils.FeatureDB: # Create or connect to GFF database.
    """
    gff_file: Path to the GFF file with normalized extension. returns a gffutils FeatureDB object.
    db_path: Path to the database file (.db) derived from the GFF file.
    If the database file does not exist, it creates one from the GFF file.
    If it exists, it connects to the existing database.
    """
    # Create db_path by replacing extensions, but keep original gff_file for reading
    db_path = gff_file.replace('.gff3', '.db').replace('.gff.gz', '.db').replace('.gff', '.db')
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
