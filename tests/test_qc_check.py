import pytest
import gffutils
from pathlib import Path
from Gene_Model_Summariser.fasta_validator import FastaChecker
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