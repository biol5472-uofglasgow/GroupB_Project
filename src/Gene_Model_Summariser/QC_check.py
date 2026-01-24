import gffutils
from typing import Optional
from .gff_parser import GFF_Parser

class QC_flags:
    # Class to generate QC flags for gene models from parser data
    def __init__(self, db: gffutils.FeatureDB, fasta: Optional[dict] = None) -> None:
        self.db = db
        self.fasta = fasta
    
    def gc_content(self, sequence: str) -> float:
        '''
        Takes in a DNA sequence string from the fasta dictionary. 
        '''
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
    
    def process_all_sequences(self) -> dict[str, dict[str, float | int]]:
        """Process all sequences in the fasta dictionary and return metrics for each chromosome."""
        if not self.fasta:
            return {}
        
        results = {}
        for chrom_id, seq_record in self.fasta.items():
            sequence = str(seq_record.seq)
            n_count, n_percent = self.N_content(sequence)
            results[chrom_id] = {
                'gc_content': self.gc_content(sequence),
                'sequence_length': self.sequence_length(sequence),
                'n_count': n_count,
                'n_percent': n_percent
            }
        return results