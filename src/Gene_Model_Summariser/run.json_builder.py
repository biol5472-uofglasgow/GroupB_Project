# functions for .json file builder
# writing the run.json file

import json
from pathlib import Path
from datetime import datetime
import logging
from typing import Any
TOOL_NAME = "Gene model summariser (transcript/gene QC summary) - Group B"
TOOL_VERSION = "1.0.0"

logger = logging.getLogger(__name__)

#Uses the runner's local timezone automatically
def whats_the_time_mr_wolf() -> str:
    return datetime.now().astimezone().isoformat(timespec="seconds")

#get the metadata to add to run.json 
def file_meta(path: str | Path) -> dict: #take in file as a string/path
    p = Path(path) #save path as a string or a path 
    if not p.is_file():
        logger.error("Expected file not found: %s", p) #no file present, tell the user and return it as empty string for the run.json
        return {"path": str(p), "bytes": None}

    return {"path": str(p), "bytes": p.stat().st_size} #get the path name as with path as string and size of file saved as bytes 


#write the json_file structure 
def write_json_file(output_path: str | Path, data: dict[str, Any]) -> None:
    with open(output_path, "w", encoding="utf-8") as file:
        json.dump(data, file, indent=2, sort_keys=True)
        file.write("\n")
