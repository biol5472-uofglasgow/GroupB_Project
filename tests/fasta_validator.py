import sys
from pathlib import Path
import Bio.SeqIO as SeqIO

class FastaChecker:
    def __init__(self, fasta_file) -> None:
        self.fasta_file = fasta_file

    def validate_fasta(self):
        """
        Validates a FASTA file format using BioPython's SeqIO.
        Performs comprehensive checks including structure, content, and format.
        
        Args:
            file_path: Path to the FASTA file
            
        Returns:
            bool: True if valid FASTA format, False otherwise
        """
        try:
            file_path = Path(self.fasta_file)
            
            #check if file exists
            if not file_path.exists():
                print(f"Error: File {file_path} not found")
                return False
            
            #check if it's actually a file, not a directory
            if not file_path.is_file():
                print(f"Error: {file_path} is not a file")
                return False
            
            #check if the file is empty
            if file_path.stat().st_size == 0:
                print(f"Error: File {file_path} is empty")
                return False
            
            has_sequence = False
            valid_chars = set('ACGTN')
            seen_ids = set()
            sequence_count = 0
            
            #SeqIO.parse() will raise an exception if format is invalid
            for record in SeqIO.parse(str(file_path), "fasta"):
                has_sequence = True
                sequence_count += 1
                
                #check for empty header
                if not record.id or record.id.strip() == '':
                    print(f"Error: Empty header found at sequence {sequence_count}")
                    return False
                
                #check for duplicate sequence IDs
                if record.id in seen_ids:
                    print(f"Error: Duplicate sequence ID '{record.id}'")
                    return False
                seen_ids.add(record.id)
                
                #check for empty sequences
                if len(record.seq) == 0:
                    print(f"Error: Empty sequence for header '{record.id}'")
                    return False
                
                #check for whitespace in sequence
                seq_str = str(record.seq)
                if ' ' in seq_str or '\t' in seq_str or '\n' in seq_str:
                    print(f"Error: Whitespace found in sequence '{record.id}'")
                    return False
                
                #check for valid characters (case insensitive)
                seq_upper = seq_str.upper()
                if not all(c in valid_chars for c in seq_upper):
                    invalid_chars = set(c for c in seq_upper if c not in valid_chars)
                    print(f"Error: Invalid characters {invalid_chars} in sequence '{record.id}'")
                    return False
            
            if not has_sequence:
                print(f"Error: No sequence data found in {file_path}")
                return False
                
            return True
            
        except Exception as e:
            print(f"Error reading FASTA file: {e}")
            return False
