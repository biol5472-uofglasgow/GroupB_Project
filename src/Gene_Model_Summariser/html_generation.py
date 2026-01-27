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

Summary Metrics Section:
Total number of genes
Total number of transcripts
Mean, median, and maximum transcripts per gene
Number and percentage of transcripts with CDS
Number and percentage of transcripts with QC flags  


Visualisations (generated using matplotlib): 

Bar chart showing the distribution of transcripts per gene(X-axis: number of transcripts per gene/Y-axis: number of genes)
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


#the function to generate the HTML report using Jinja2 templating (will be saved into a separate HTML generation file(groupB.html.j2) once finalised)
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


    #function to compute summary metrics from the transcript summary DataFrame
    def compute_summary_metrics(df: pd.DataFrame) -> dict:
    # Function to compute summary metrics from the transcript summary DataFrame
    total_genes = int(df["gene_id"].nunique()) #calculate the total number of unique gene IDs
    total_transcripts = int(len(df)) #calculate the total number of transcripts (rows in the DataFrame)
    
    # calculate transcripts per gene statistics: mean, median, maximum 
    transcript_per_gene = df.groupby("gene_id")["transcript_id"].nunique()  #how many genes have how many transcripts
    transcript_mean = float(transcript_per_gene.mean()) if len(transcript_per_gene) else 0.0 #calculate mean transcripts per gene
    transcript_median = float(transcript_per_gene.median()) if len(transcript_per_gene) else 0.0 #calculate median transcripts per gene
    transcript_max = int(transcript_per_gene.max()) if len(transcript_per_gene) else 0.0 #calculate maximum transcripts per gene

    # calculate percentage of transcripts with has_cds = true
    has_cds_count = int(df["has_cds"].sum()) #count how many transcripts have has_cds = true
    has_cds_percent = (has_cds_count / total_transcripts) * 100 if total_transcripts > 0 else 0.0 #calculate percentage of transcripts with has_cds = true

    # calculate percentage of transcripts with QC flags (flags column not empty)
    flagged_transcripts_count = int(df[df["flags"].notna() & (df["flags"] != "")].shape[0]) #count transcripts with non-empty flags
    flagged_transcripts_percent = (flagged_transcripts_count / total_transcripts) * 100 if total_transcripts > 0 else 0.0 #calculate percentage of flagged transcripts
    
    return {
        "total_genes": total_genes,
        "total_transcripts": total_transcripts,
        "transcript_mean": transcript_mean,                 
        "transcript_median": transcript_median,
        "transcript_max": transcript_max,
        "has_cds_count": has_cds_count,
        "has_cds_percent": has_cds_percent,
        "flagged_transcripts_count": flagged_transcripts_count,
        "flagged_transcripts_percent": flagged_transcripts_percent,
    }

#function to generate bar chart data for QC flag types
def compute_qc_flag_count(df: pd.DataFrame) -> dict[str, int]: 
    #function to compute the distribution of QC flag types from the transcript summary DataFrame
    flag_counts = {}  #dictionary to hold counts of each flag type

    for flags in df["flags"].dropna(): #iterate over non-null flags in the DataFrame
        for flag in flags.split(","): #split multiple flags by comma
            flag = flag.strip() #remove leading/trailing whitespace
            if flag: #if the flag is not empty
                flag_counts[flag] = flag_counts.get(flag, 0) + 1 #increment the count for this flag type

    return flag_counts #return the dictionary of flag counts


#function to generate histogram data for exon counts 
def compute_exon_count_for_histogram(df: pd.DataFrame) -> dict[int, int]: 
    # Function to compute the distribution of exon counts from the transcript summary DataFrame
    exon_count_distribution = {}  #dictionary to hold counts of each exon count

    for n_exons in df["n_exons"].dropna():  #iterate over non-null exon counts in the DataFrame
        exon_count = int(n_exons)  #convert exon count to integer
        exon_count_distribution[exon_count] = exon_count_distribution.get(exon_count, 0) + 1  #increment the count for this exon count

    return exon_count_distribution  #return the dictionary of exon count distribution

#function to generate bar chart data for transcripts per gene distribution
def compute_transcripts_per_gene_distribution(df: pd.DataFrame) -> dict[int, int]:
    # Function to compute the distribution of transcripts per gene from the transcript summary DataFrame
    transcripts_per_gene = df.groupby("gene_id")["transcript_id"].nunique()  #group by gene_id and count unique transcripts

    distribution = {}  #dictionary to hold counts of each transcripts per gene value

    for count in transcripts_per_gene:  #iterate over the counts of transcripts per gene
        count = int(count)  #convert count to integer
        distribution[count] = distribution.get(count, 0) + 1  #increment the count for this transcripts per gene value

    return distribution  #return the dictionary of transcripts per gene distribution

#function to compute counts of flagged vs unflagged transcripts
def compute_flagged_vs_unflagged(df: pd.DataFrame) -> dict[str, int]:
    flagged = int((flags_clean != "").sum()) #count transcripts with non-empty flags
    unflagged = int(len(df) - flagged) #calculate unflagged transcripts by subtracting flagged from total
    return {"flagged": flagged, "unflagged": unflagged} #return the counts as a dictionary



