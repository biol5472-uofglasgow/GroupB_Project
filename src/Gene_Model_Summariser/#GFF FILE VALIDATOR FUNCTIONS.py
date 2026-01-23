#GFF FILE VALIDATOR FUNCTIONS
#use line_length_checker when actually parsing the gff file to check for any errors 
#check each feature in the gff db for the following:
# skip hashtag lines/blank lines to check every valid line has 9 columns 

def line_length_checker(line: str, line_number: int):
    stripped = line.strip()
    if stripped == "" or stripped.startswith("#"):
        return None
    columns = stripped.split("\t")
    if len(columns) != 9:
        logger.error(
            f"Line {line_number}: Expected 9 tab-separated columns, found {len(columns)}. Line was: {stripped}")
        return None
    return columns

#required fields (seqid, source, type)
def validate_required_fields(feature) -> bool:
    if feature.seqid is None or str(feature.seqid).strip() == "":
        logger.error(f"Missing seqid for feature {feature.id}")
        return False
    if feature.source is None or str(feature.source).strip() == "":
        logger.error(f"Missing source for feature {feature.id}")
        return False
    if feature.featuretype is None or str(feature.featuretype).strip() == "":
        logger.error(f"Missing type for feature {feature.id}")
        return False
    return True

#coordinates (start, end)
def validate_coordinates(feature) -> bool:
    if feature.start is None or feature.start == "":
        logger.error(f"Feature missing start on feature {feature.id}")
        return False
    if feature.end is None or feature.end == "":
        logger.error(f"Feature missing end on feature {feature.id}")
        return False
    else:
        try:
            if int(feature.start) > int(feature.end):
                logger.error(f"Start is bigger than end: {feature.start} > {feature.end}")
                return False
        except (TypeError, ValueError):
            logger.error(f"Invalid start/end values: start={feature.start}, end={feature.end}")
            return False
    return True

#strand (+, -, .)
def validate_strand(feature) -> bool:
    if feature.strand in {"+", "-", "."}:
        return True
    else:
        logger.error(f"Invalid strand value for feature {feature.id}: {feature.strand}")
        return False

#score (float or .)
def validate_score(feature) -> bool:
    if feature.score is None or feature.score == ".":
        return True
    try:
        float(feature.score)
        return True
    except (TypeError, ValueError):
        logger.error(f"Invalid score value for feature {feature.id}: {feature.score}")
        return False

#phase (0, 1, 2 or .)
def validate_phase(feature) -> bool:
    # phase often appears as ".", None, "0"/"1"/"2", or 0/1/2 (int)
    if feature.phase in {None, ".", "0", "1", "2", 0, 1, 2}:
        return True
    else:
        logger.error(f"Invalid phase value for feature {feature.id}: {feature.phase}")
        return False

#attributes (key=value pairs) - need to write the code to check this properly

def check_db(db) -> bool:
    data_ok = True
    for feature in db.all_features():
        # 1) required fields
        if validate_required_fields(feature) is False:
            data_ok = False
        # 2) coordinates
        if validate_coordinates(feature) is False:
            data_ok = False
        # 3) strand
        if validate_strand(feature) is False:
            data_ok = False
        # 4) score
        if validate_score(feature) is False:
            data_ok = False
        # 5) phase
        if validate_phase(feature) is False:
            data_ok = False
    # True means everything passed, False means at least one feature failed
    return data_ok