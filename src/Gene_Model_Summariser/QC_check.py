import gffutils
from typing import Optional
from Bio.Seq import Seq
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
    
    def CDS_start_codon_check(self, sequence: str, cds_start, cds_end) -> bool:
        # Function to check if the CDS starts with a valid start codon (ATG)
        if not sequence or cds_start is None or cds_end is None:
            return False
        start_codon = sequence[cds_start-1:cds_start+2].upper()  # Adjust for 0-based indexing
        return start_codon == 'ATG'

    def CDS_stop_codon_check(self, sequence: str, cds_start, cds_end) -> bool:
        # Function to check if the CDS ends with a valid stop codon (TAA, TAG, TGA)
        if not sequence or cds_start is None or cds_end is None:
            return False
        stop_codon = sequence[cds_end-3:cds_end].upper()  # Adjust for 0-based indexing
        return stop_codon in {'TAA', 'TAG', 'TGA'}
    
    def gff_QC(self) -> dict[str, list[str]]:
        model = GFF_Parser(self.db).transcript_model()
        gff_flags = {}
        for transcript_id, features in model.items():
            gff_flags[transcript_id] = []
            '''counting exons'''
            exon_count = len(features['exon(s)'])
            if exon_count > 5:
                exon_flag = 'exon_count>5'
                gff_flags[transcript_id].append(exon_flag)
            exon_positions = [(exon.start, exon.end) for exon in features['exon(s)']]
            exon_positions.sort()
            for i in range(1, len(exon_positions)):
                if exon_positions[i][0] < exon_positions[i-1][1]:
                    overlaps = 'overlapping_exons'
                    gff_flags[transcript_id].append(overlaps)
                    break
            '''checking CDS start codon'''
            if self.fasta:
                chrom_id = features['gene'].chrom
                seq_record = self.fasta.get(chrom_id)
                if seq_record:
                    sequence = str(seq_record.seq)
                    strand = seq_record.strand
                    if strand == '-':
                        sequence = str(Seq(sequence).reverse_complement())
                    if features['CDS(s)']:
                        cds_count = 1
                        for feature in features['CDS(s)']:
                            cds_start = feature.start
                            cds_end = feature.end
                            if not self.CDS_start_codon_check(sequence, cds_start, cds_end):
                                start_codon_flag = f'invalid_start_codon'
                                gff_flags[transcript_id].append(start_codon_flag)
                            if not self.CDS_stop_codon_check(sequence, cds_start, cds_end):
                                stop_codon_flag = f'invalid_stop_codon'
                                gff_flags[transcript_id].append(stop_codon_flag)
                            cds_count += 1             
        return gff_flags
                        