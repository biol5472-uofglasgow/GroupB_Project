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
#look into using genome2json package? https://pypi.org/project/genome2json/ - seems to convert gff+fasta to json 


def write_json_file(output_path: str, data: dict) -> None:

