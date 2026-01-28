# functions for .json file builder
# writing the run.json file

import json
from pathlib import Path
from datetime import datetime
import logging
from typing import Any, Optional
TOOL_NAME = "Gene model summariser (transcript/gene QC summary) - Group B"
TOOL_VERSION = "1.0.0"

logger = logging.getLogger(__name__)

#Uses the runner's local timezone automatically
def whats_the_time_mr_wolf() -> str:
    return datetime.now().astimezone().isoformat(timespec="milliseconds")

#get the metadata to add to run.json 
def file_meta(path: str | Path | None) -> dict: #take in file as a string/path
    if path is None:
        return {"path": None, "bytes": None} #fasta is optional - return this 

    p = Path(path) #save path as a string or a path 
    if not p.is_file():
        logger.error("Expected file not found: %s", p) #no file present, tell the user and return it as empty string for the run.json
        return {"path": str(p), "bytes": None}

    return {"path": str(p), "bytes": p.stat().st_size} #get the path name as with path as string and size of file saved as bytes 


#write the json_file structure 
def write_json_file(output_path: str | Path, data: dict[str, Any]) -> None:
    with open(output_path, "w", encoding="utf-8") as file:
        json.dump(data, file, indent=2)
        file.write("\n")

#run the json file including writing out the tool name/version, start time/input files and outputs 
def build_run_json(start_time: str, gff_file: Path,fasta_file: Optional[Path], output_dir: Path,
    results_filename: str = "results.tsv", html_filename: str = "results.html") -> dict[str, Any]:
    
    output_dir = Path(output_dir) #ensure output_dir is a Path object
    results_path = output_dir / results_filename #grab results.tsv path (will make sure to give full directory link when finished)
    html_path = output_dir / html_filename #grab html_filename (will make sure to give full directory link when finished)

    return {
        "tool": {"name": TOOL_NAME, "version": TOOL_VERSION},
        "timestamp": {"start": start_time, "end": None},
        "inputs": {"gff": file_meta(gff_file),"fasta": file_meta(fasta_file)},
        "outputs": {"results_tsv": {"path": str(results_path), "bytes": None}, "results_html": {"path": str(html_path), "bytes": None},
        },
    }

#add what time the file ends after running / size of the output files once built 
def update_end_time_and_output_sizes(output_dict: dict[str, Any]) -> dict[str, Any]:
    #set the end time
    output_dict["timestamp"]["end"] = whats_the_time_mr_wolf()

    # check each output file and store its size
    script_outputs = ["results_tsv", "results_html"]

    #for the script outputs file, loop over each to grab its path and calculate its size 
    for output_key in script_outputs:
        file_path_as_text = output_dict["outputs"][output_key]["path"]
        file_path = Path(file_path_as_text)

        # if the file exists, get its size. if not, keep bytes as None
        if file_path.is_file():
            output_dict["outputs"][output_key]["bytes"] = file_path.stat().st_size
        else:
            output_dict["outputs"][output_key]["bytes"] = None

    return output_dict


#create an initial run.json at the start of the pipeline. will be imported into main()
#this writes tool metadata, a start timestamp, input file metadata,and placeholder output metadata
def make_run_json_file(gff_file: Path, fasta_file: Optional[Path], output_dir: Path, results_filename: str = "results.tsv", 
                       html_filename: str = "results.html",run_filename: str = "run.json") -> Path:
    
    output_dir = Path(output_dir) #ensure output_dir is a Path
    run_path = output_dir / run_filename #full path to the run.json file

    start_time = whats_the_time_mr_wolf() #record the start time
    #build the run.json dictionary structure
    run_dict = build_run_json(start_time=start_time, gff_file=Path(gff_file),fasta_file=Path(fasta_file),
        output_dir=output_dir,results_filename=results_filename,html_filename=html_filename,)

    #write initial run.json (end time unknown yet) as not finished yet
    write_json_file(run_path, run_dict)

    #return the path
    return run_path 

#once eveerything has run in main(), capture the output files, time they were made and how big they are 
def finalise_run_json_file(output_dir: Path, run_filename: str | Path = "run.json") -> Path:
    output_dir = Path(output_dir) #ensure output_dir is a Path
    run_path = output_dir / run_filename #full path to the run.json file

    #load the previously written run.json created by make_run_json_file
    #although these functions are next to eachother, they will sandwich the main logic to store start/end time and input/output files
    with open(run_path, "r", encoding="utf-8") as file:
        run_dict = json.load(file)

    #update end time and compute output file sizes
    run_dict = update_end_time_and_output_sizes(run_dict)
    write_json_file(run_path, run_dict)
    return run_path

