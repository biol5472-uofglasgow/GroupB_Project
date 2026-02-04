'''
Docstring for Gene_Model_Summariser.html_generation
This module contains functions to generate HTML reports from the outputs of Pillar 1 (transcript summary TSV and run.json).
It uses Jinja2 templating to create the HTML structure and pandas to process the TSV data
'''

from pathlib import Path
from datetime import datetime
from jinja2 import Environment, FileSystemLoader, select_autoescape
import pandas as pd
import json
import matplotlib.pyplot as plt

####################################################################################################################################################################################
#function used to generate the HTML report using Jinja2 templating
############################################################################################################################################################################################
#the function to generate the HTML report using Jinja2 templating (will be saved into a separate HTML generation file(groupB.html.j2) once finalised)
def generate_html_report(report_data: dict, template_dir: Path) -> str:
    template_dir = Path(template_dir)

    env = Environment(
        loader=FileSystemLoader(str(template_dir)),
        autoescape=select_autoescape(["html", "xml"]),
    )

  #load the HTML Jinja2 template from template_dir
    template_name = "groupB.html.j2"
    try:
        template = env.get_template(template_name)
    except Exception as e:
        raise FileNotFoundError(
            f"Could not find template '{template_name}' in {template_dir}"
        ) from e

    #report_data dict is available in the template as {{ data }}
    html_output = template.render(data=report_data)

    return html_output


####################################################################################################################################################################################
#building functions to open and extract data from tsv and json files
#loader for run.json for HTML 
####################################################################################################################################################################################

def load_outputs(output_dir: str | Path) -> tuple[pd.DataFrame, dict]:
    output_dir = Path(output_dir)

    tsv_path = output_dir / "transcript_summary.tsv"
    json_path = output_dir / "run.json"

    if not tsv_path.exists():
        raise FileNotFoundError(f"Missing transcript summary TSV: {tsv_path}")
    if not json_path.exists():
        raise FileNotFoundError(f"Missing run metadata JSON: {json_path}")

    df = pd.read_csv(tsv_path, sep="\t")

    required = {"gene_id","transcript_id","exon_count","has_cds","flags"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"transcript_summary.tsv missing columns: {sorted(missing)}")


    with json_path.open("r", encoding="utf-8") as f:
        run_info = json.load(f)

    return df, run_info

####################################################################################################################################################################################
#parser for run.json for HTML header 
####################################################################################################################################################################################


def parse_runjson_time(ts: str) -> datetime:
    return datetime.fromisoformat(ts)

def build_provenance(run_info: dict) -> dict:
    tool_name = run_info.get("tool", {}).get("name")
    tool_version = run_info.get("tool", {}).get("version")

    start_script = run_info.get("timestamp", {}).get("start")
    end_script = run_info.get("timestamp", {}).get("end")

    duration_seconds = None
    if start_script and end_script:
        try:
            duration_seconds = (parse_runjson_time(end_script) - parse_runjson_time(start_script)).total_seconds()
        except Exception:
            duration_seconds = None

    return {
        "tool_name": tool_name,
        "tool_version": tool_version,
        "start": start_script,
        "end": end_script,
        "duration_seconds": duration_seconds}

####################################################################################################################################################################################
#functions to compute various metrics from the transcript summary DataFrame
####################################################################################################################################################################################

#function to compute summary metrics from the transcript summary DataFrame
def compute_summary_metrics(df: pd.DataFrame) -> dict:
    # Function to compute summary metrics from the transcript summary DataFrame
    total_genes = int(df["gene_id"].nunique())  # calculate the total number of unique gene IDs
    total_transcripts = int(len(df))  # calculate the total number of transcripts (rows in the DataFrame)

    # calculate transcripts per gene statistics: mean, median, maximum
    transcript_per_gene = df.groupby("gene_id")["transcript_id"].nunique()  # number of transcripots per gene
    transcript_mean = round(float(transcript_per_gene.mean()), 2) if len(transcript_per_gene) else 0.0  # calculate mean transcripts per gene (2 d.p.)
    transcript_median = round(float(transcript_per_gene.median()), 2) if len(transcript_per_gene) else 0.0  # calculate median transcripts per gene (2 d.p.)
    transcript_max = int(transcript_per_gene.max()) if len(transcript_per_gene) else 0  # calculate maximum transcripts per gene

    # calculate percentage of transcripts with has_cds = true
    # count how many transcripts have has_cds = true
    has_cds_bool = df["has_cds"].astype(str).str.lower().isin(["true"])  # make sure tsv can be compatible strings and bool
    has_cds_count = int(has_cds_bool.sum())
    has_cds_percent = round((has_cds_count / total_transcripts) * 100, 2) if total_transcripts > 0 else 0.0  # calculate percentage of transcripts with has_cds = true (2 d.p.)

    # calculate percentage of transcripts with QC flags (flags column not empty)
    flags_clean = df["flags"].fillna("").astype(str).str.strip()  # clean flags to make sure stripped of white space and remove any NANs with '' for counting
    flagged_transcripts_count = int((flags_clean != "").sum())  # count flags
    flagged_transcripts_percent = round((flagged_transcripts_count / total_transcripts) * 100, 2) if total_transcripts > 0 else 0.0  # calculate percentage flagged (2 d.p.)

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
def compute_exon_count(df: pd.DataFrame) -> dict[int, int]: 
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
    flags_clean = df["flags"].fillna("").astype(str).str.strip() #clean flags column by converting to string and stripping whitespace
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

#function to plot exon count distribution (bar chart)
def plot_exon_count_histogram(exon_count_distribution: dict[int, int], outpath: Path) -> str:
    exon_counts = sorted(exon_count_distribution.keys()) #sorted list of exon counts
    transcript_counts = [exon_count_distribution[x] for x in exon_counts] #corresponding transcript counts

    plt.figure(figsize=(8, 4)) #set figure size
    plt.bar(exon_counts, transcript_counts) #create bar chart
    plt.xlabel("Number of exons") #label x-axis
    plt.ylabel("Number of transcripts") #label y-axis
    plt.title("Exon count distribution") #set chart title
    plt.xticks(exon_counts)  # show each exon count on the x-axis
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
#function to run functions and save the data to a dictionary to be extracted for plotting
#also built the table for the HTML file 
def compute_report_stats(df: pd.DataFrame) -> dict:
    #summary numbers for the report
    metrics = compute_summary_metrics(df)
    metrics_table_records = summary_metrics_table(metrics).to_dict(orient="records")

    #data needed to make each plot
    exon_count_histogram_data = compute_exon_count(df)
    transcripts_per_gene_bar_data = compute_transcripts_per_gene_distribution(df)

    flagged_vs_unflagged_bar_data = compute_flagged_vs_unflagged(df)
    qc_flag_counts_per_transcript_data = compute_qc_flag_count_per_transcript(df)

    #package everything into one dictionary
    report_stats = {}
    report_stats["summary_metrics"] = metrics
    report_stats["summary_metrics_table"] = metrics_table_records

    report_stats["plot_inputs"] = {}
    report_stats["plot_inputs"]["exon_count_histogram_data"] = exon_count_histogram_data
    report_stats["plot_inputs"]["transcripts_per_gene_bar_data"] = transcripts_per_gene_bar_data
    report_stats["plot_inputs"]["flagged_vs_unflagged_bar_data"] = flagged_vs_unflagged_bar_data
    report_stats["plot_inputs"]["qc_flag_counts_per_transcript_data"] = qc_flag_counts_per_transcript_data

    return report_stats


#used to save the report figures
#takes in the data from the report_stats and majke the output folder for the figures for the HTML to embed them
def save_report_figures(plot_inputs: dict, output_dir: Path) -> dict[str, str]:
    # Create a folder called "figures" inside the output directory
    figures_dir = output_dir / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True) 

    #Plot 1: Exon count distribution - make the plot and save it as a .png 
    exon_count_plot_filename = plot_exon_count_histogram(plot_inputs["exon_count_histogram_data"],
        figures_dir / "exon_count_distribution.png")

    #Plot 2: Transcripts per gene distribution - make the plot and save it as a .png  
    transcripts_per_gene_plot_filename = plot_transcripts_per_gene_distribution(plot_inputs["transcripts_per_gene_bar_data"],
        figures_dir / "transcripts_per_gene_distribution.png")

    #Plot 3: Flagged vs unflagged transcripts -make the plot and save it as a .png 
    flagged_vs_unflagged_plot_filename = plot_flagged_vs_unflagged(plot_inputs["flagged_vs_unflagged_bar_data"],
        figures_dir / "flagged_vs_unflagged.png")

    #Plot 4: QC flags for one count per transcript - make the plot and save it as a .png
    qc_flags_plot_filename = plot_qc_flag_counts_per_transcript(plot_inputs["qc_flag_counts_per_transcript_data"],
        figures_dir / "qc_flags_per_transcript.png")

    #return the image filenames - these are saved to the same root directory as the transcript_summary.tsv but in a seperate directory
    #these figs will then be embedded into the HTML file 
    return {
        "exon_count_plot": exon_count_plot_filename,
        "transcripts_per_gene_plot": transcripts_per_gene_plot_filename,
        "flagged_vs_unflagged_plot": flagged_vs_unflagged_plot_filename,
        "qc_flags_per_transcript_plot": qc_flags_plot_filename}

#####################################################################################################################
#Once all built in Python, put into data dictionary in Jinja2 format
####################################################################################################################
#used to build a string with the report data, summary metrics and links to raw .tsv and .json files 
#as well as qc flag definitions
def build_report_data(report_stats: dict, figures: dict) -> dict:
    return {
        "summary_metrics": report_stats["summary_metrics"],
        "summary_metrics_table": report_stats["summary_metrics_table"],
        "qc_flag_definitions": QC_FLAG_DEFINITIONS,
        "qc_flag_names": QC_FLAG_NAMES,
        "figures": figures,
        "artefacts": {"results_tsv": "transcript_summary.tsv", "run_json": "run.json"},
    }


# this is used for the CLI endpoint to generate the report
def run_report(output_dir: Path, template_dir: Path | None = None) -> Path:
    output_dir = Path(output_dir)  # output directory

    # If template_dir not provided, use the Gene_Model_Summariser package directory
    if template_dir is None:
        template_dir = Path(__file__).resolve().parent  # .../Gene_Model_Summariser/
    else:
        template_dir = Path(template_dir)

    df, run_info = load_outputs(output_dir)  # reads transcript_summary.tsv + run.json and load them in
    report_stats = compute_report_stats(df)  # compute stats for the reports
    figures = save_report_figures(report_stats["plot_inputs"], output_dir)  # save figures into this report

    # generates the report_data dictionary for Jinja2 loading
    report_data = build_report_data(report_stats, figures)

    # load data in from run.json
    report_data["run_info"] = run_info
    report_data["provenance"] = build_provenance(run_info)

    html = generate_html_report(report_data, template_dir=template_dir)  # loads groupB.html.j2 from template_dir

    # load HTML report into output dir
    out_html = output_dir / "report.html"
    out_html.write_text(html, encoding="utf-8")
    return out_html