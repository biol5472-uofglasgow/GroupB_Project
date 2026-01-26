# GroupB_Project

A command-line tool to summarize gene models and output basic GC metrics.

**Project Selection:** Project 5  
**Group Members:** John Hardin, Hans Henrik Norberg, Dom Thompson  
**Built for:** BIOL5472 Course at the University of Glasgow

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

3. Run the tool
```bash
groupb.py --gff data/models.gff --fasta data/ref.fasta --outdir results/
```

## Docker
1. Pull from docker hub
```bash
docker pull beyondourminds/groupb-tool:latest
```

2. Run the tool
```bash
docker run -v FilePathToData/data:/data beyondourminds/groupb-tool:latest -g /data/gffFile.gff -f /data/fastaFile.fasta

# or if data files are in your current working directory

docker run -v $(pwd):/data beyondourminds/groupb-tool:latest -g /data/gffFile.gff -f /data/fastaFile.fasta
```

### For Development

1. Clone the repository:
```bash
git clone <repository-url>
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
pip install git+https://github.com/yourusername/GroupB_Project.git
```

## Usage

### Canonical Run Command -- to be updated as project is built
```bash
GroupB-tool
```

For help and available options:
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
