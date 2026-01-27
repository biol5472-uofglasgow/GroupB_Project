'''
Docstring for Gene_Model_Summariser.html_generation
This module contains functions to generate HTML reports from the outputs of Pillar 1 (transcript summary TSV and run.json).
It uses Jinja2 templating to create the HTML structure and pandas to process the TSV data
'''

from pathlib import Path
from jinja2 import Environment, FileSystemLoader, select_autoescape
import pandas as pd
import json
import matplotlib.pyplot as plt

####################################################################################################################################################################################
#function used to generate the HTML report using Jinja2 templating
############################################################################################################################################################################################
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
        raise FileNotFoundError(f"Could not find template '{template_name}' in {pillar3_folder}") from e

    # tsv_output will be available in Jinja as {{ data }}
    html_output = template.render(data=tsv_output)

    return html_output

####################################################################################################################################################################################
#building a function to open and extract data from tsv and json files
####################################################################################################################################################################################

def load_pillar1_outputs(pillar1_dir: Path) -> tuple[pd.DataFrame, dict]:
    
    pillar1_dir = Path(pillar1_dir) #ensure pillar1_dir is a Path object

    tsv_path = pillar1_dir / "transcript_summary.tsv" #construct the full path to the transcript summary TSV file
    
    json_path = pillar1_dir / "run.json" #construct the full path to the run.JSON file

    df = pd.read_csv(tsv_path, sep="\t") #read the transcript summary TSV into a pandas DataFrame

    run_info = json.loads(json_path.read_text(encoding="utf-8")) #load the contents of run.json into a Python dictionary
    return df, run_info 

####################################################################################################################################################################################
#functions to compute various metrics from the transcript summary DataFrame
####################################################################################################################################################################################

    #function to compute summary metrics from the transcript summary DataFrame
def compute_summary_metrics(df: pd.DataFrame) -> dict:
    # Function to compute summary metrics from the transcript summary DataFrame
    total_genes = int(df["gene_id"].nunique()) #calculate the total number of unique gene IDs
    total_transcripts = int(len(df)) #calculate the total number of transcripts (rows in the DataFrame)
    
    # calculate transcripts per gene statistics: mean, median, maximum 
    transcript_per_gene = df.groupby("gene_id")["transcript_id"].nunique() #how many genes have how many transcripts
    transcript_mean = float(transcript_per_gene.mean()) if len(transcript_per_gene) else 0.0 #calculate mean transcripts per gene
    transcript_median = float(transcript_per_gene.median()) if len(transcript_per_gene) else 0.0 #calculate median transcripts per gene
    transcript_max = int(transcript_per_gene.max()) if len(transcript_per_gene) else 0 #calculate maximum transcripts per gene

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

#function to create a summary metrics table DataFrame from the computed metrics - this will be used in the HTML report
def summary_metrics_table(metrics: dict) -> pd.DataFrame:
    rows = [
        ("Total genes", metrics["total_genes"]),
        ("Total transcripts", metrics["total_transcripts"]),
        ("Mean transcripts per gene", metrics["transcript_mean"]),
        ("Median transcripts per gene", metrics["transcript_median"]),
        ("Max transcripts per gene", metrics["transcript_max"]),
        ("Transcripts with CDS (count)", metrics["has_cds_count"]),
        ("Transcripts with CDS (%)", metrics["has_cds_percent"]),
        ("Flagged transcripts (count)", metrics["flagged_transcripts_count"]),
        ("Flagged transcripts (%)", metrics["flagged_transcripts_percent"]),
    ]
    return pd.DataFrame(rows, columns=["metric", "value"]) 

####################################################################################################################################################################################
#functions to compute data for visualisations from the transcript summary DataFrame
####################################################################################################################################################################################

#function to generate histogram data for exon counts 
def compute_exon_count_for_histogram(df: pd.DataFrame) -> dict[int, int]: 
    # Function to compute the distribution of exon counts from the transcript summary DataFrame
    exon_count_distribution = {}  #dictionary to hold counts of each exon count
    for count_exon in df["exon_count"].dropna():  #iterate over non-null exon counts in the DataFrame
        exon_count = int(count_exon)  #make sure exon_count is an integer
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
    flags_clean = df["flags"].astype(str).str.strip() #clean flags column by converting to string and stripping whitespace
    flagged = int((flags_clean != "").sum()) #count transcripts with non-empty flags
    unflagged = int(len(df) - flagged) #calculate unflagged transcripts by subtracting flagged from total
    return {"flagged": flagged, "unflagged": unflagged} #return the counts as a dictionary

#dictionary defining QC flag descriptions
QC_FLAG_DEFINITIONS = {
    "exon_count>5": "Transcript has more than 5 exons.",
    "overlapping_exons": "At least two exons overlap in genomic coordinates.",
    "invalid_CDS_phase": "A CDS feature has an invalid phase/frame (not 0/1/2).",
    "N_in_CDS": "CDS sequence contains one or more 'N' bases.",
    "ambiguous_bases_in_CDS": "CDS contains bases outside A/C/G/T/N.",
    "CDS_not_multiple_of_3": "Total CDS length is not divisible by 3.",
    "invalid_start_codon": "CDS does not start with ATG.",
    "invalid_stop_codon": "CDS does not end with TAA/TAG/TGA.",
    "CDS_too_short": "CDS length is < 3 bases.",
    "no_CDS": "Transcript has no CDS features.",
}

QC_FLAG_NAMES = list(QC_FLAG_DEFINITIONS.keys()) #list of all QC flag names

# function to count how many transcripts have each flag
def compute_qc_flag_count_per_transcript(df: pd.DataFrame) -> dict[str, int]:
    flag_counts: dict[str, int] = {} # dictionary to hold counts of each flag type
    for flags in df["flags"].fillna("").astype(str): # iterate over flags, treating NaN as empty string
        flag_set = {f.strip() for f in flags.split(",") if f.strip() != ""} # split and clean flags
        for flag in flag_set: # count each unique flag once per transcript
            flag_counts[flag] = flag_counts.get(flag, 0) + 1 # increment the count for this flag type
    return flag_counts

####################################################################################################################################################################################
#functions to build visualisations from the compute data functions above
####################################################################################################################################################################################

#function to plot transcripts per gene distribution bar chart
def plot_exon_count_histogram(exon_count_distribution: dict[int, int], outpath: Path) -> str:
    exon_counts = sorted(exon_count_distribution.keys()) #sorted list of exon counts
    transcript_counts = [exon_count_distribution[x] for x in exon_counts] #corresponding transcript counts

    plt.figure(figsize=(8, 4)) #set figure size
    plt.bar(exon_counts, transcript_counts) #create bar chart
    plt.xlabel("Number of exons") #label x-axis
    plt.ylabel("Number of transcripts") #label y-axis
    plt.title("Exon count distribution") #set chart title
    plt.xticks(exon_counts)  # show each exon count on the x-axis (long labels in horizontal will disrupt the chart)
    plt.tight_layout() #adjust layout to prevent clipping
    plt.savefig(outpath, dpi=200) #save the figure 
    plt.close() #close the file and return 
    return outpath.name  # return image filename for embedding in HTML


#function to plot transcripts per gene distribution bar chart
def plot_transcripts_per_gene_distribution(distribution: dict[int, int], outpath: Path) -> str:
    transcripts_per_gene = sorted(distribution.keys()) #sorted list of transcripts per gene counts
    gene_counts = [distribution[x] for x in transcripts_per_gene] #corresponding gene counts

    plt.figure(figsize=(8, 4)) #set figure size
    plt.bar(transcripts_per_gene, gene_counts) #create bar chart
    plt.xlabel("Transcripts per gene") #label x-axis
    plt.ylabel("Number of genes") #label y-axis
    plt.title("Transcripts per gene distribution") #set chart title
    plt.xticks(transcripts_per_gene)  # show each integer count on the x-axis (long labels in horizontal will disrupt the chart)
    plt.tight_layout() #adjust layout to prevent clipping
    plt.savefig(outpath, dpi=200) #save the figure 
    plt.close() #close the file and return 
    return outpath.name  # return image filename for embedding in HTML

#function to plot flagged vs unflagged transcripts bar chart
def plot_flagged_vs_unflagged(counts: dict[str, int], outpath: Path) -> str:
    labels = ["flagged", "unflagged"] #labels for the two categories
    values = [counts["flagged"], counts["unflagged"]] #corresponding counts

    plt.figure(figsize=(6, 4)) #set figure size
    plt.bar(labels, values) #create bar chart
    plt.xlabel("Transcript status") #label x-axis
    plt.ylabel("Number of transcripts") #label y-axis
    plt.title("Flagged vs unflagged transcripts") #set chart title
    plt.tight_layout() #adjust layout to prevent clipping
    plt.savefig(outpath, dpi=200) #save the figure 
    plt.close() #close the file and return 
    return outpath.name  # return image filename for embedding in HTML

#function to plot qc flag counts per transcript bar chart
def plot_qc_flag_counts_per_transcript(flag_counts: dict[str, int], outpath: Path) -> str:
    flags = QC_FLAG_NAMES  # fixed order from your definitions
    counts = [flag_counts.get(f, 0) for f in flags] # corresponding counts

    plt.figure(figsize=(10, 4)) #set figure size
    plt.bar(flags, counts) #create bar chart
    plt.xlabel("QC flag") #label x-axis
    plt.ylabel("Number of transcripts") #label y-axis
    plt.title("QC flag counts (unique per transcript)") #set chart title
    plt.xticks(rotation=45, ha="right") #rotate x-axis labels for readability (long labels in horizontal will disrupt the chart)
    plt.tight_layout() #adjust layout to prevent clipping
    plt.savefig(outpath, dpi=200) #save the figure 
    plt.close() #close the file and return 
    return outpath.name  # return image filename for embedding in HTML

####################################################################################################################################################################################
#putting it all together
#load -> compute -> save plots -> build report data
####################################################################################################################################################################################
if __name__ == "__main__":
    # 1) Load Pillar 1 outputs
    pillar1_dir = Path("path/to/pillar1_outputs") #this needs changed to include Johns output for the pillar1_dir (will check when finished)
    df, run_info = load_pillar1_outputs(pillar1_dir)

    # 2) Compute plot data
    exon_count_distribution = compute_exon_count_for_histogram(df)
    transcripts_per_gene_distribution = compute_transcripts_per_gene_distribution(df)
    flagged_vs_unflagged_counts = compute_flagged_vs_unflagged(df)
    qc_flag_counts_per_transcript = compute_qc_flag_count_per_transcript(df)
    metrics = compute_summary_metrics(df)
    metrics_table = summary_metrics_table(metrics).to_dict(orient="records") 


    # 3) Make a figures directory + save plots into the directory 
    figures_dir = pillar1_dir / "figures" 
    figures_dir.mkdir(parents=True, exist_ok=True)

    exon_count_distribution_plot_file = plot_exon_count_histogram(exon_count_distribution, figures_dir / "exon_count_distribution.png")
    transcript_per_gene_plot_file = plot_transcripts_per_gene_distribution(transcripts_per_gene_distribution, figures_dir / "transcripts_per_gene_distribution.png")
    flagged_plot_file = plot_flagged_vs_unflagged(flagged_vs_unflagged_counts, figures_dir / "flagged_vs_unflagged.png")
    qc_per_transcript_plot_file = plot_qc_flag_counts_per_transcript(qc_flag_counts_per_transcript, figures_dir / "qc_flags_per_transcript.png")

    # 4) Build the dictionary that Jinja2 will use
    report_data = {
        "run_info": run_info,
        "metrics": metrics,
        "metrics_table": metrics_table,
        "plots": {
            "exon_count": exon_count_distribution_plot_file,
            "transcripts_per_gene": transcript_per_gene_plot_file,
            "flagged_vs_unflagged": flagged_plot_file,
            "qc_flags_per_transcript": qc_per_transcript_plot_file,
        },
    }

    # 5) Render + write HTML
    html = generate_html_report(report_data)
    (pillar1_dir / "report.html").write_text(html, encoding="utf-8")