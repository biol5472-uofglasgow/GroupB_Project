import pytest
import gffutils
from pathlib import Path
from Gene_Model_Summariser.gff_parser import GFF_Parser

@pytest.fixture
def gff_db_fixture(tmp_path):
    """
    Create a GFF database from a fixture GFF file for testing.
    """
    gff_fixture = Path(__file__).parent / "Fixtures" / "models.gff3"
    db_path = tmp_path / "test.db"
    db = gffutils.create_db(str(gff_fixture), dbfn=str(db_path), force=True, keep_order=True)
    return db

class TestGFFParser:

    def test_initialization(self, gff_db_fixture):
        """
        Test that the parser initializes correctly with a GFF database.
        """
        parser = GFF_Parser(gff_db_fixture)
        assert parser.db is not None
    
    def test_get_genes(self, gff_db_fixture):
        """
        Test retrieval of gene features from the GFF database.

        models.gff3 fixture contains 1 gene feature called "gene1" with featuretype "gene". get_genes() also tries to find 'protein_coding_gene' if no 'gene' features found.
        """
        parser = GFF_Parser(gff_db_fixture)
        genes = parser.get_genes()
        assert len(genes) == 1
        assert genes[0].id == "gene1"
        assert genes[0].featuretype == "gene"
    
    def test_get_transcripts(self, gff_db_fixture):
        """
        Test retrieval of transcript features for a given gene ID.

        models.gff3 fixture contains 2 mRNA features (tx1, tx2) under gene1.
        """
        parser = GFF_Parser(gff_db_fixture)
        transcripts = parser.get_transcripts("gene1")
        assert len(transcripts) == 2
        assert "tx1" in [t.id for t in transcripts]
        assert "tx2" in [t.id for t in transcripts]
    
    def test_get_exons(self, gff_db_fixture):
        """
        Test retrieval of exon features for a given transcript ID.

        models.gff3 fixture: tx1 has 3 exons, tx2 has 2 exons.
        """
        parser = GFF_Parser(gff_db_fixture)
        exons_tx1 = parser.get_exons("tx1")
        assert len(exons_tx1) == 3

        exons_tx2 = parser.get_exons("tx2")
        assert len(exons_tx2) == 2
    
    def test_get_cds(self, gff_db_fixture):
        """
        Test retrieval of CDS features for a given transcript ID.

        models.gff3 fixture: tx1 has 3 CDS features, tx2 has no CDS features (edge case).
        """
        parser = GFF_Parser(gff_db_fixture)
        cds_tx1 = parser.get_cds("tx1")
        assert len(cds_tx1) == 3

        cds_tx2 = parser.get_cds("tx2")
        assert len(cds_tx2) == 0

    def test_count_exons(self, gff_db_fixture):
        """
        Test counting of exon features for a given transcript ID. Test is similar to test_get_exons, but counting is the intended purpose for count_exons().

        models.gff3 fixture: tx1 has 3 exons, tx2 has 2 exons.
        """
        parser = GFF_Parser(gff_db_fixture)
        assert parser.count_exons("tx1") == 3
        assert parser.count_exons("tx2") == 2