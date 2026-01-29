from dataclasses import dataclass
from pathlib import Path

@dataclass
class TranscriptWithFlags:
    chrom: str
    start: int
    end: int
    transcript_id: str
    qc_flags: set[str]
    strand: str

def write_qc_bed(transcripts: list[TranscriptWithFlags], output_path: Path) -> None:
    """
    Write flagged transcripts to a BED file for visualization in genome browsers.
    Only includes transcripts with QC flags (transcripts that failed quality checks).
    
    Uses a two-pass approach:
    1. First pass: Scan all transcripts to discover which QC flags exist in the dataset
    2. Second pass: Assign colors dynamically and write BED file
    
    This makes the function general-purpose - it works with any QC flags without
    needing to hardcode flag names in advance.
    
    Args:
        transcripts: List of TranscriptWithFlags objects to process
        output_path: Path where the BED file will be written
    
    Output format:
        BED9 format with RGB colors for each QC flag type
    """
    #first pass -> collect all unique flags from the dataset
    all_flags = set()
    for transcript in transcripts:
        all_flags.update(transcript.qc_flags)
    
    #define a color palette to cycle through
    color_palette = [
        "255,0,0",      # red
        "255,165,0",    # orange
        "255,0,255",    # magenta
        "0,0,255",      # blue
        "0,255,0",      # green
        "255,255,0",    # yellow
        "0,255,255",    # cyan
        "128,0,128",    # purple
    ]
    
    #assign colors to flags alphabetically
    #using modulo to cycle through colors if more flags than colors
    sorted_all_flags = sorted(all_flags)
    flag_colors: dict[str, str] = {}
    for i, flag in enumerate(sorted_all_flags):
        flag_colors[flag] = color_palette[i % len(color_palette)]
    
    #second pass -> write BED file
    with output_path.open("w") as file:
        file.write("track name='QC Flagged Transcripts' description='Transcripts with QC flags' itemRgb='On'\n")
        
        for transcript in transcripts:
            if not transcript.qc_flags:
                continue
            
            sorted_flags = sorted(transcript.qc_flags)
            flags_str: str = "|".join(sorted_flags)
            name: str = f"{transcript.transcript_id}/{flags_str}"
            color: str = flag_colors.get(sorted_flags[0], "128,128,128")
            #set thick start and end to match BED format requirements
            thick_start: int = transcript.start
            thick_end: int = transcript.end
            
            #score is set to 0 as it's not used here
            file.write(f"{transcript.chrom}\t{transcript.start}\t{transcript.end}\t{name}\t0\t{transcript.strand}\t{thick_start}\t{thick_end}\t{color}\n")