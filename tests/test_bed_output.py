from pathlib import Path
from Gene_Model_Summariser.qc_flags_bed import TranscriptWithFlags, write_qc_bed

def test_bed_output_only_flagged(tmp_path: Path) -> None:
    """
    Test that BED output includes only transcripts with QC flags.

    Verifies filtering logic: transcripts without flags should be excluded.

    Args:
        tmp_path: Pytest fixture providing temporary directory for the test files.
    """
    #create test data with one flagged and one unflagged transcript
    transcripts: list[TranscriptWithFlags] = [
        TranscriptWithFlags("chr1", 1000, 5000, "t1", {"NO_CDS"}, "+"), 
        TranscriptWithFlags("chr1", 6000, 8000, "t2", set(), "+"),
    ]

    output: Path = tmp_path / "test.bed"
    write_qc_bed(transcripts, output)

    lines: list[str] = output.read_text().strip().split("\n")
    #expected: header line + 1 data line (only t1, which has flags)
    assert len(lines) == 2
    assert "t1" in lines[1]
    assert "t2" not in lines[1]

def test_bed_format_correct(tmp_path: Path) -> None:
    """
    Test that BED output conforms to BED9 format specification.
    
    Validates all nine required BED fields are present and correctly formatted:
    chrom, chromStart, chromEnd, name, score, strand, thickStart, thickEnd, itemRgb.
    
    Args:
        tmp_path: Pytest fixture providing temporary directory for test files
    """
    transcripts: list[TranscriptWithFlags] = [
        TranscriptWithFlags('chr1', 1000, 5000, 't1', {'NO_CDS', 'SHORT'}, '+')
    ]

    output: Path = tmp_path / "test.bed"
    write_qc_bed(transcripts, output)

    lines: list[str] = output.read_text().strip().split("\n")
    data_line: str = lines[1]
    fields: list[str] = data_line.split("\t")

    #verify BED9 format -> 9 fields
    assert len(fields) == 9
    #verify individual fields against expected values.
    assert fields[0] == 'chr1', "Chromosome field incorrect"
    assert fields[1] == '1000', "Start position field incorrect"
    assert fields[2] == '5000', "End position field incorrect"
    assert 'NO_CDS' in fields[3], "Name field should contain NO_CDS flag"
    assert 'SHORT' in fields[3], "Name field should contain SHORT flag"
    assert fields[4] == '0', "Score field should be 0"
    assert fields[5] == '+', "Strand field incorrect"
    assert fields[6] == '1000', "ThickStart should match start position"
    assert fields[7] == '5000', "ThickEnd should match end position"
    assert ',' in fields[8], "RGB color field should contain comma-separated values"