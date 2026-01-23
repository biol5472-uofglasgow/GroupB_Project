import sys
import os
import gffutils
import pytest
from pathlib import Path
from GroupB_Project5 import load_gff_database

def validate_gff_file(args):
    #convert string to path object -> cross platform compatibility ("\" on Linux/Mac and "/" on Windows)
    file_path = Path(args)
    
    #checking that the file exists using exists() method
    if not file_path.exists():
        sys.exit(f"Error: File not found: {args}")
    
    #checking that its actually a file using the is_file() method
    if not file_path.is_file():
        sys.exit(f"Error: Path is not a file: {args}")
    
    #checking file extension against allowed formats -> lower() to make it case insensitive
    #suffix method automatically extracts the file extension
    if file_path.suffix.lower() not in {".gff", ".gff3", ".gtf"}:
        sys.exit(f"Error: File format not supported: {args}. Must be GFF, GFF3, or GTF format.")

def validate_fasta_file(args):
    file_path = Path(args)
    
    if not file_path.exists():
        sys.exit(f"Error: File not found: {args}")

    if not file_path.is_file():
        sys.exit(f"Error: Path is not a file: {args}")
    #adding tester for fasta files
    if file_path.suffix.lower() not in {".fasta", ".fa", ".fna"}:
        sys.exit(f"Error: File format not supported: {args}. Must be FASTA format.")
    