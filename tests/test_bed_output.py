from pathlib import Path
from Gene_Model_Summariser.qc_flags_bed import TranscriptWithFlags, write_qc_bed

def test_bed_output_only_flagged(tmp_path: Path) -> None:
    """
    Docstring for test_bed_output_only_flagged
    
    :param tmp_path: Description
    """
    transcripts: list[TranscriptWithFlags] = [
        TranscriptWithFlags("chr1", 1000, 5000, "t1", {"NO_CDS"}, "+"), 
        TranscriptWithFlags("chr1", 6000, 8000, "t2", set(), "+"),
    ]

    output: Path = tmp_path / "test.bed"
    write_qc_bed(transcripts, output)

    lines: list[str] = output.read_text().strip().split("\n")
    assert len(lines) == 2
    assert "t1" in lines[1]
    assert "t2" not in lines[1]

def test_bed_format_correct(tmp_path: Path) -> None:
    """
    Docstring for test_bed_format_correct
    
    :param tmp_path: Description
    """
    transcripts: list[TranscriptWithFlags] = [
        TranscriptWithFlags('chr1', 1000, 5000, 't1', {'NO_CDS', 'SHORT'}, '+')
    ]

    output: Path = tmp_path / "test.bed"
    write_qc_bed(transcripts, output)

    lines: list[str] = output.read_text().strip().split("\n")
    data_line: str = lines[1]
    fields: list[str] = data_line.split("\t")

    assert len(fields) == 9
    assert fields[0] == 'chr1'
    assert fields[1] == '1000'
    assert fields[2] == '5000'
    assert 'NO_CDS' in fields[3]
    assert 'SHORT' in fields[3]
    assert fields[4] == '0'
    assert fields[5] == '+'
    assert fields[6] == '1000'
    assert fields[7] == '5000'
    assert ',' in fields[8]