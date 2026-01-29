from pathlib import Path
#importing Set for type hint -> e.g.: seen_ids: Set[str]
from typing import Set
import Bio.SeqIO as SeqIO




class FastaChecker:
    def __init__(self, fasta_file, logger) -> None:
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
        
        Args:
            file_path: Path to the FASTA file
            
        Raises:
            FastaValidationError: If file fails validation checks
            
        Example:
            >>> validate_fasta("genome.fasta")
            >>> validate_fasta("bad.fasta")  #raises FastaValidationError
        """
        logger = self.logger
        file_path = Path(self.fasta_file)
        Valid = True

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
                    Valid = False

                
                #check for duplicate IDs
                if record.id in seen_ids:
                    logger.error(f"Duplicate sequence ID: '{record.id}'")
                    Valid = False
                seen_ids.add(record.id)
                
                #check for empty sequences
                if len(record.seq) == 0:
                    logger.error(
                        f"Empty sequence for ID: '{record.id}'"
                    )
                    Valid = False

                #check for invalid characters (case-insensitive)
                seq_upper: str = str(record.seq).upper()
                invalid_chars: Set[str] = set(seq_upper) - valid_chars
                
                if invalid_chars:
                    logger.error(
                        f"Invalid characters {invalid_chars} in sequence '{record.id}'"
                    )
                    Valid = False

        except ValueError as e:
            #biopython raises ValueError for malformed FASTA
            logger.error(
                f"Malformed FASTA format: {e}"
            )
            Valid = False
        
        #raise validation error if no sequences found
        if sequence_count == 0:
            logger.error("No sequences found in file")
            Valid = False
        
        return Valid

    def fasta_parse(self):
        logger = self.logger
        try:
            fasta = SeqIO.to_dict(SeqIO.parse(self.fasta_file, 'fasta'))
            return fasta
        except Exception as e:
            logger.error(f"Error parsing FASTA file: {e}")
            return None
