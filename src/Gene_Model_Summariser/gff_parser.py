import gffutils
from typing import Optional

class GFF_Parser:

    def __init__(self, db: gffutils.FeatureDB) -> None:
        self.db = db

    def get_genes(self) -> list[gffutils.Feature]:
        """Retrieve all gene features from the GFF database."""
        # Try common gene feature types
        genes = list(self.db.features_of_type('gene'))
        if not genes:
            genes = list(self.db.features_of_type('protein_coding_gene'))
        return genes
    
    def get_transcripts(self, gene_id: str) -> list[gffutils.Feature]:
        """Retrieve all transcript features for a given gene ID."""
        transcripts = list(self.db.children(gene_id, featuretype='mRNA', order_by='start'))
        return transcripts
    
    def get_exons(self, transcript_id: str) -> list[gffutils.Feature]:
        """Retrieve all exon features for a given transcript ID."""
        exons = list(self.db.children(transcript_id, featuretype='exon', order_by='start'))
        return exons
    
    def get_cds(self, transcript_id: str) -> list[gffutils.Feature]:
        """Retrieve all CDS features for a given transcript ID."""
        cds_features = list(self.db.children(transcript_id, featuretype='CDS', order_by='start'))
        return cds_features
    
    def count_exons(self, transcript_id: str) -> int:
        """Count the number of exons for a given transcript ID."""
        exons = list(self.db.children(transcript_id, featuretype='exon'))
        exon_count = len(exons)
        return exon_count
    
    def check_cds(self, transcript_id: str) -> bool:
        """Check if a given transcript has associated CDS features."""
        has_cds = False
        cds_features = list(self.db.children(transcript_id, featuretype='CDS'))
        if cds_features:
            has_cds = True
        return has_cds
    
    def tsv_output(self) -> dict:
        """putting together into one clean model"""
        output_dict = {}
        genes = self.get_genes()
        for gene in genes:
            if gene.id:
                transcripts = self.get_transcripts(gene.id)
            else:
                continue
            for transcript in transcripts:
                if transcript.id:
                    transcript_id = transcript.id
                    exon_count = self.count_exons(transcript.id)
                    has_cds = self.check_cds(transcript.id)
                    output_dict[transcript.id] = {
                        'gene_id': gene.id,
                        'transcript_id': transcript_id,
                        'exon_count': exon_count,
                        'has_cds': has_cds,
                        'chrom': transcript.chrom,
                        'start': transcript.start,
                        'end': transcript.end,
                        'strand': transcript.strand
                    }
        return output_dict
    
    def transcript_model(self) -> dict:
        """Generates a summary model of transcripts with exon counts and CDS presence."""
        model = {}
        genes = self.get_genes()
        for gene in genes:
            if gene.id:
                transcripts = self.get_transcripts(gene.id)
                for transcript in transcripts:
                    if transcript.id:
                        exons = self.get_exons(transcript.id)
                        cds_features = self.get_cds(transcript.id)
                        model[transcript.id] = {
                            'gene': gene.id,
                            'exon(s)': exons,
                            'CDS(s)': cds_features
                        }
        return model