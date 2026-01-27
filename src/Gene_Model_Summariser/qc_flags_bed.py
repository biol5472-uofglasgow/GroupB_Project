@dataclass
class TranscriptWithFlags:
    chrom: str
    start: int
    end: int
    transcript_id: str
    qc_flags: set[str]
    strand: str

def write_qc_bed(transcripts: List[TranscriptWithFlags], output_path: Path) -> None:
    """
    Write flagged transcripts to a BED file.
    Only includes transcripts with QC flags.
    """
    #defining colors for the flags -> RGB format
    flag_colors = {
        "flag1": "255,0,0",
        "flag2": "0,255,0",
        "flag3": "0,0,255",
    }

    with output_path.open("w") as file:
        #BED track header with itemRgb enabled for color visualization
        file.write("track name='QC Flagged Transcripts' description='Transcripts with QC flags' itemRgb='On'\n")

        for transcript in transcripts:
            if not transcript.qc_flags:
                continue

            flags_str = ",".join(sorted(transcript.qc_flags))
            name = f"{transcript.transcript_id}/{flags_str}"

            #selecting flag color based on the first flag
            color = flag_colors.get(list(transcript.qc_flags)[0], "128,128,128")
            thick_start = transcript.start
            thick_end = transcript.end

            file.write(f"{transcript.chrom}\t{transcript.start}\t{transcript.end}\t{name}\t0\t{transcript.strand}\t{thick_start}\t{thick_end}\t{color}\n")