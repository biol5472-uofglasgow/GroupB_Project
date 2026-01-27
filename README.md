# GroupB_Project

A command-line tool to summarize gene models and output basic GC metrics.

**Project Selection:** Project 5  
**Group Members:** John Hardin, Hans Henrik Norberg, Dom Thompson  
**Built for:** BIOL5472 Course at the University of Glasgow

## Current Version
v1.0.0

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
docker pull beyondourminds/groupb-tool:latest
```

### For Development

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

4. Install the package in editable mode:
```bash
pip install -e .
```

### For End Users

Simply install the package:
```bash
pip install .
```

Or install directly from GitHub (once published):
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
- If no arguments provided, defaults to the directory of the inputted gff file as results/run_00# where # is the current run number

### Conda
```bash
groupb.py --gff data/models.gff --fasta data/ref.fasta --outdir results/
```

### Docker
For docker, the input file directory must be mounted using the -v command as shown:
```bash
docker run -v FilePathToData/data:/data beyondourminds/gene-summariser:latest -g /data/gffFile.gff -f /data/fastaFile.fasta
```
if data files are in your current working directory
```bash
docker run -v $(pwd):/data beyondourminds/gene-summariser:latest -g /data/gffFile.gff -f /data/fastaFile.fasta
```

### For help and available options:
```bash
GroupB-tool --help
```

## Dependencies

- pandas >= 2.3
- gffutils >= 0.12
- Python >= 3.13
- matplotlib >= 3.10
- seaborn >= 0.13
- biopython >= 1.85
