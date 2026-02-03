# GroupB_Project

A command-line tool to summarize gene models and output basic QC metrics.

**Project Selection:** Project 5  
**Group Members:** John Hardin, Hans Henrik Norberg, Dom Thompson  
**Built for:** BIOL5472 Course at the University of Glasgow

## Current Version
v1.1.0

## About
### What it is
This tool is used to analyse a GFF file, and optionally a FASTA file, and output QC metrics for flagged transcripts. 

### Outputs
This tool produces 6 outputs:
1. transcript_summary.tsv:
- This is your basic output file that outputs a ror each transcript in the format: gene_id, transcript_id, n_exons, has_cds, chrom, start, end, strand, flags
2. qc_flags.gff3:
- This file is a copy of relevant transcript features with QC flags embedded as an attribute. 
- Output in the standard GFF file format and loadable in genome browsers.
3. qc_flagged.bed:
- Outputs a BED file containing the loci of flagged items for quick browsing
- Loadable in genome browsers with color coordination
4. report.html:
- An HTML file that provides a visual report of the flagged and unflagged transcripts
5. figures:
- a folder containing all of the figures made/used in the HTML report
6. run.json
- a run file that contains a record of the tool, timestamp, inputs, fasta file (if provided), outputs, and HTML result file. 

### Assumptions
This program makes a few assumptions when processing QC flags that should be considered when using this tool
1. Only accepts ATG as a valid start codon
   * will flag as invalid_start_codon if a different start codon is used
2. Terminal CDS must end in a stop codon
   * If the terminal CDS does not end in a stop codon, it will be flagged as invalid_stop_codon
3. Assumes phase values are correct within the GFF file
   * If the phase value is incorrect or incomplete, this will affect the previous flags
Violations of these assumptions are recorded as QC flags but may reflect annotation conventions rather than biological errors


## Installation

### Conda Installer : 
1. Create and activate the conda environment:
```bash
conda env create -f environment.yml
conda activate biol5472_groupB
```
2. Install the tool:
```bash
pip install -e . --no-deps
```

### Docker:
1. Pull from docker hub
```bash
docker pull beyondourminds/gene-summariser:latest
```

### Pip

1. Clone the repository:
```bash
git clone https://github.com/biol5472-uofglasgow/GroupB_Project.git
cd GroupB_Project
```

2. Create a virtual environment:
```bash
python -m venv .venv
```

3. Activate the virtual environment:

**Windows (PowerShell):**
```powershell
.venv\Scripts\Activate.ps1
```

**Windows (Command Prompt):**
```cmd
.venv\Scripts\activate.bat
```

**Linux/Mac:**
```bash
source .venv/bin/activate
```
### For Developers
4. Install the package in editable mode:
```bash
pip install -e .
```

### For End Users

4. Simply install the package:
```bash
pip install .
```

Or install directly from GitHub:
```bash
pip install git+https://github.com/biol5472-uofglasgow/GroupB_Project.git
```

## Usage

### Input Commands
1. -g or --gff (required)
- takes in the path for the gff file you wish to use
2. -f or --fasta (optional)
- takes in the path for the fasta file you wish to use
3. -o or --outdir (optional)
- Takes in the desired directory for output
- If no arguments provided, defaults to the directory of the inputted gff file as results/run_# where # is the current run number

### Conda
```bash
GroupB-tool --gff data/models.gff --fasta data/ref.fasta --outdir results/
```

### Docker
For docker, the input file directory must be mounted using the -v command as shown:
```bash
docker run -v FilePathToData/data:/data beyondourminds/gene-summariser:latest -g /data/gffFile.gff -f /data/fastaFile.fasta
```
If data files are in your current working directory
```bash
docker run -v $(pwd):/data beyondourminds/gene-summariser:latest -g /data/gffFile.gff -f /data/fastaFile.fasta
```

### For help and available options:
Run tool with no provided arguments, or provide the --help command

## Dependencies

- pandas >= 2.3
- gffutils >= 0.12
- Python >= 3.13
- matplotlib >= 3.10
- seaborn >= 0.13
- biopython >= 1.85
- jinja2 >= 3.1.0
