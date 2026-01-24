import gffutils
from typing import Optional

class GFF_Parser:

    def __init__(self, db: gffutils.FeatureDB) -> None:
        self.db = db

    def get_genes(self) -> list[gffutils.Feature]:
        """Retrieve all gene features from the GFF database."""
        genes = list(self.db.features_of_type('gene'))
        return genes
    
    def get_transcripts(self, gene_id: str) -> list[gffutils.Feature]:
        """Retrieve all transcript features for a given gene ID."""
        transcripts = list(self.db.children(gene_id, featuretype='mRNA', order_by='start'))
        return transcripts
    
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
    
    def put_together(self) -> dict:
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
                    exon_count = self.count_exons(transcript.id)
                    has_cds = self.check_cds(transcript.id)
                    output_dict[transcript.id] = {
                        'gene_id': gene.id,
                        'exon_count': exon_count,
                        'has_cds': has_cds
                    }
        return output_dict