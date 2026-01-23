#functions for .json file
import json
#writing the run.json file 
#what do i want to inclde in the run.json file?:
#inputs (gff file, fasta file) 
#outputs (transcript_summary.tsv, qc_flags.gff3, OPTINONAL : added.gff3, removed.gff3, changed.gff3/summary.tsv)
#gff file counts (genes, transcripts, exons, cds features)
#distributions (transcripts per gene, exons per transcript)
#qc (flag counts) - when these are defined we can return these afterwards
#parameters (min_exon_length, max_exon_length, min_cds_lengt, max_exon_length) 

def write_json_file(output_path: str, data: dict) -> None:
    with open(output_path, 'w') as run_json_file:
        json.dump(data, run_json_file, indent=4)
    print(f"JSON file written to {output_path}")
    