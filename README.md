# GroupB_Project

A command-line tool to summarize gene models and output basic GC metrics.

**Project Selection:** Project 5  
**Group Members:** John Hardin, Hans Henrik Norberg, Dom Thompson  
**Built for:** BIOL5472 Course at the University of Glasgow

## Installation

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
