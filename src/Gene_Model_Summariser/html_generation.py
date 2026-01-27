#pillar 3 - HTML Report Generation

'''
PILLAR 3 – STATIC HTML REPORT CONTENT
Provenance / Run Details (top of report)
Tool name and version
Timestamp report was generated
Command used and parameters (including QC thresholds)
Input GFF3 file path
Input GFF3 file hash (e.g. SHA256, if stored)
Links to raw Pillar 1 output files:
../transcript_summary.tsv
../run.json
../qc_flags.gff3 or ../qc_flags.bed (if present)
Main Transcript Summary Table (from transcript_summary.tsv)
Derived Summary Metrics (from transcript_summary.tsv)
Total number of genes (unique gene_id values)
Total number of transcripts (number of rows in the TSV)
Transcripts per gene: mean, median, and maximum
Percentage of transcripts with has_cds = true
Percentage of transcripts with QC flags (flags column not empty)
Top 3 most common QC flags (with counts)
Visualisations (generated using matplotlib)
Bar chart showing the distribution of transcripts per gene
X-axis: number of transcripts per gene
Y-axis: number of genes
Histogram showing the distribution of exon counts across transcripts
Values taken from the n_exons column
Bar chart showing counts per QC flag type
One bar per flag
Y-axis: number of transcripts with that QC issue
Side-by-side bar chart comparing flagged vs unflagged transcripts
QC Flag Definitions Table
Table listing each QC flag
Description of what each flag means
'''

from pathlib import Path
from jinja2 import Environment, FileSystemLoader, select_autoescape
import pandas as pd
import json


def generate_html_report(tsv_output: dict) -> str:  
    # tsv_output will be renamed once Pillar 1 tsv_output dict is finished and finalised

    # Get the directory of the current file and set as template folder for Jinja2
    pillar3_folder = Path(__file__).resolve().parent  

    # Set up Jinja2 environment from the folder
    env = Environment(
        loader=FileSystemLoader(str(pillar3_folder)),
        autoescape=select_autoescape(["html", "xml"]),
    )

    # Load the HTML template file from the pillar3 folder
    template_name = "groupB.html.j2"
    try:
        template = env.get_template(template_name)
    except Exception as e:
        raise FileNotFoundError(
            f"Could not find template '{template_name}' in {pillar3_folder}"
        ) from e

    # tsv_output will be available in Jinja as {{ data }}
    html_output = template.render(data=tsv_output)

    return html_output



    #building a function to open and extract data from tsv and json files (will be tested once these are finished)
    #once the directory containing this information is built, i will change the path from pillar1_dir to the correct path
            
    def load_pillar1_outputs(pillar1_dir: Path) -> tuple[pd.DataFrame, dict]:
    pillar1_dir = Path(pillar1_dir) #ensure pillar1_dir is a Path object

    tsv_path = pillar1_dir / "transcript_summary.tsv" #construct the full path to the transcript summary TSV file
    
    json_path = pillar1_dir / "run.json" #construct the full path to the run.JSON file

    df = pd.read_csv(tsv_path, sep="\t") #read the transcript summary TSV into a pandas DataFrame
    
    run_info = json.loads(json_path.read_text(encoding="utf-8")) #load the contents of run.json into a Python dictionary
    return df, run_info 














