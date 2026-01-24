import gffutils
from typing import Optional
from Bio.Seq import Seq
try:
    from .gff_parser import GFF_Parser
    from .fasta_validator import FastaChecker
except ImportError:
    from gff_parser import GFF_Parser
    from fasta_validator import FastaChecker

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
    
    def contains_N(self, sequence: str) -> bool:
        # Function to check if a sequence contains 'N' bases
        if not sequence:
            return False
        return 'N' in sequence.upper()
    
    def ambiguous_bases(self, sequence: str) -> bool:
        # Function to check for ambiguous bases other than A, C, G, T, N
        if not sequence:
            return False
        valid_bases = set('ACGTN')
        sequence = sequence.upper()
        for base in sequence:
            if base not in valid_bases:
                return True
        return False
    
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
            if self.fasta:
                # Get gene feature from database using gene ID
                gene_id = features['gene']
                gene_feature = self.db[gene_id]
                chrom_id = gene_feature.seqid
                seq_record = self.fasta.get(chrom_id)
                if seq_record:
                    chrom_sequence = str(seq_record.seq)
                    strand = gene_feature.strand
                    
                    if features['CDS(s)']:
                        # Sort CDS features by their start position
                        # Create list of (start_position, cds_feature) tuples
                        cds_with_positions = []
                        for cds in features['CDS(s)']:
                            cds_with_positions.append((cds.start, cds))
                        
                        # Sort by start position (first element of tuple)
                        cds_with_positions.sort()
                        
                        # Extract just the CDS features
                        cds_list = [cds for start, cds in cds_with_positions]
                        
                        # For negative strand, process in reverse order
                        if strand == '-':
                            cds_list = cds_list[::-1]
                        
                        # Build complete CDS sequence with phase adjustment
                        cds_seq = ''
                        for cds in cds_list:
                            # Extract CDS segment from chromosome
                            cds_segment = chrom_sequence[cds.start-1:cds.end]
                            
                            # Reverse complement if on negative strand
                            if strand == '-':
                                cds_segment = str(Seq(cds_segment).reverse_complement())
                            
                            # Apply phase offset (skip bases at start)
                            phase = int(cds.frame) if cds.frame != '.' else 0
                            if phase not in {0, 1, 2}:
                                gff_flags[transcript_id].append('invalid_CDS_phase')
                            cds_seq += cds_segment[phase:]
                        
                        # Check for issues in the complete CDS sequence
                        if self.contains_N(cds_seq):
                            gff_flags[transcript_id].append('N_in_CDS')
                        if self.ambiguous_bases(cds_seq):
                            gff_flags[transcript_id].append('ambiguous_bases_in_CDS')
                        
                        #check length multiple of 3
                        if len(cds_seq) % 3 != 0:
                            gff_flags[transcript_id].append('CDS_not_multiple_of_3')
                        
                        # Check start codon (first 3 bases of CDS)
                        if len(cds_seq) >= 3:
                            start_codon = cds_seq[:3].upper()
                            if start_codon != 'ATG':
                                gff_flags[transcript_id].append('invalid_start_codon')
                        
                        # Check stop codon (last 3 bases of CDS)
                        if len(cds_seq) >= 3:
                            stop_codon = cds_seq[-3:].upper()
                            if stop_codon not in {'TAA', 'TAG', 'TGA'}:
                                gff_flags[transcript_id].append('invalid_stop_codon')
                    else:
                        gff_flags[transcript_id].append('no_CDS')
                        
        return gff_flags

gff_file = r"C:\Users\jtspy\Desktop\Python\BCPyAssessment\PlasmoDB-54_Pfalciparum3D7.gff"
db_path = gff_file.replace('.gff', '.db').replace('.gff3', '.db').replace('.gff.gz', '.db')
db = gffutils.create_db(gff_file, dbfn=db_path, force=True, keep_order=True)
fasta_file = r"C:\Users\jtspy\Desktop\Python\BCPyAssessment\PlasmoDB-54_Pfalciparum3D7_Genome.fasta"
fasta_checker = FastaChecker(fasta_file)
fasta = fasta_checker.fasta_parse()
qc = QC_flags(db, fasta)
results = qc.gff_QC()
print(results)             