import pytest
import gffutils
from pathlib import Path
from Gene_Model_Summariser.QC_check import QC_flags

@pytest.fixture
def gff_db_fixture(tmp_path):
    """
    Create a GFF database from a fixture GFF file for testing.
    """
    gff_fixture = Path(__file__).parent / "Fixtures" / "models.gff3"
    db_path = tmp_path / "test.db"
    db = gffutils.create_db(str(gff_fixture), dbfn=str(db_path), force=True, keep_order=True)
    return db


@pytest.fixture
def fasta_file_fixture():
    """
    Return path to fixture FASTA file.
    """
    from Bio.SeqIO import to_dict, parse
    fasta_path = Path(__file__).parent / "Fixtures" / "ref.fasta"
    return to_dict(parse(str(fasta_path), 'fasta'))

class TestQCCheck:

    def test_gff_QC(self, gff_db_fixture):
        """
        Test the gff_QC method to ensure it generates QC flags correctly.
        """

        flags = QC_flags(gff_db_fixture).gff_QC()
        assert "exon_count>5" in flags.get("tx3", [])
        assert "overlapping_exons" in flags.get("tx3", [])
        assert "no_CDS" in flags.get("tx2", [])
    
    def test_gff_QC_with_fasta(self, gff_db_fixture, fasta_file_fixture):
        """
        Test the gff_QC method with FASTA data to ensure it generates QC flags correctly.
        """

        qc_checker = QC_flags(gff_db_fixture, fasta_file_fixture)
        flags = qc_checker.gff_QC()
        
        # Test tx3 has multiple QC flags
        tx3_flags = flags.get("tx3", [])
        assert "exon_count>5" in tx3_flags
        assert "overlapping_exons" in tx3_flags
        assert "invalid_CDS_phase" in tx3_flags
        assert "N_in_CDS" in tx3_flags
        assert "invalid_start_codon" in tx3_flags
        assert "invalid_stop_codon" in tx3_flags
        
        # Test tx2 has no_CDS flag
        assert "no_CDS" in flags.get("tx2", [])
