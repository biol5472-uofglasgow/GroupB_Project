from pathlib import Path
#importing Set for type hint -> e.g.: seen_ids: Set[str]
from typing import Set
import Bio.SeqIO as SeqIO




class FastaChecker:
    """
    Validates and parses FASTA files for gene model quality control.
    
    Provides comprehensive format validation and sequence parsing with
    detailed error logging for debugging pipeline failures.
    """

    def __init__(self, fasta_file, logger) -> None:
        """
        Initialize FASTA checker with file path and logger.
        
        Args:
            fasta_file: Path to FASTA file (str or Path)
            logger: Logger instance for error reporting
        """
        self.fasta_file = fasta_file
        self.logger = logger



    #accepting file path as Path or str
    def validate_fasta(self) -> bool:
        """
        Validates a FASTA file format using BioPython's SeqIO.
        
        Performs comprehensive checks:
        - At least one sequence present
        - No duplicate sequence IDs
        - No empty sequences or headers
        - Valid nucleotide characters (ACGTN, case-insensitive)
            
        Returns:
            bool: True if valid, False otherwise.
            
        Note:
            Logs specific errors for each validation failure but continues
            checking to report all issues in a single run.
        """
        logger = self.logger
        file_path = Path(self.fasta_file)
        valid = True

        #validation state tracking with type hints
        seen_ids: Set[str] = set()
        sequence_count: int = 0
        valid_chars: Set[str] = set('ACGTN')
        
        #parse sequences - BioPython will raise ValueError if malformed
        #try/except block here to catch exceptions where they happen
        try:
            for record in SeqIO.parse(str(file_path), "fasta"):
                sequence_count += 1
                
                #check for empty or whitespace-only header
                if not record.id or not record.id.strip():
                    logger.error(
                        f"Empty header at sequence {sequence_count}"
                    )
                    valid = False

                
                #check for duplicate IDs
                if record.id in seen_ids:
                    logger.error(f"Duplicate sequence ID: '{record.id}'")
                    valid = False
                seen_ids.add(record.id)
                
                #check for empty sequences
                if len(record.seq) == 0:
                    logger.error(
                        f"Empty sequence for ID: '{record.id}'"
                    )
                    valid = False

                #check for invalid characters (case-insensitive)
                seq_upper: str = str(record.seq).upper()
                invalid_chars: Set[str] = set(seq_upper) - valid_chars
                
                if invalid_chars:
                    logger.error(
                        f"Invalid characters {invalid_chars} in sequence '{record.id}'"
                    )
                    valid = False

        except ValueError as e:
            #biopython raises ValueError for malformed FASTA
            logger.error(
                f"Malformed FASTA format: {e}"
            )
            valid = False
        
        #ensure at least one sequence was found
        if sequence_count == 0:
            logger.error("No sequences found in file")
            valid = False
        
        return valid

    def fasta_parse(self):
        """
        Parse FASTA file into dictionary mapping sequence IDs to SeqRecord objects.
        
        Returns:
            dict: Dictionary of {sequence_id: SeqRecord} or None if parsing fails
            
        Note:
            Should only be called after validate_fasta() confirms file is valid.
        """
        logger = self.logger
        try:
            fasta = SeqIO.to_dict(SeqIO.parse(self.fasta_file, 'fasta'))
            return fasta
        except Exception as e:
            logger.error(f"Error parsing FASTA file: {e}")
            return None
