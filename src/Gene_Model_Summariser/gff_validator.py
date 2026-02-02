#GFF FILE VALIDATOR FUNCTIONS
#use line_length_checker when actually parsing the gff file to check for any errors 
#check each feature in the gff db for the following:
# skip hashtag lines/blank lines to check every valid line has 9 columns 

import logging

logger = logging.getLogger("GroupB_logger")

def line_length_checker(line: str, line_number: int):
    stripped = line.strip()
    if stripped == "" or stripped.startswith("#"):
        return None
    columns = stripped.split("\t")
    if len(columns) != 9:
        logger.error(f"Line {line_number}: Expected 9 tab-separated columns, found {len(columns)}. Line was: {stripped}")
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
    if feature.frame in {None, ".", "0", "1", "2", 0, 1, 2}:
        return True
    else:
        logger.error(f"Invalid phase value for feature {feature.id}: {feature.frame}")
        return False

# check attributes features
def validate_attributes(feature) -> bool:
    attributes = getattr(feature, "attributes", None)
    feature_id = getattr(feature, "id", "unknown")

    if not isinstance(attributes, dict):
        logger.error(f"Invalid attributes for feature {feature_id}")
        return False

    for key, value in attributes.items():
        if not key or value in (None, "", []):
            logger.error(f"Invalid attribute {key} for feature {feature_id}")
            return False
    return True

#validator - rerturns TRrue is everything passes otherwise false
#logs a summary of how many features(rows) and checks(one of the functions above) did not pass the tests 
def check_db(db) -> bool:
    total = 0
    failed_features = 0
    failed_checks = 0

    for feature in db.all_features():
        total += 1
        feature_failed = False

        # 1) required fields
        if not validate_required_fields(feature):
            failed_checks += 1
            feature_failed = True

        # 2) coordinates
        if not validate_coordinates(feature):
            failed_checks += 1
            feature_failed = True

        # 3) strand
        if not validate_strand(feature):
            failed_checks += 1
            feature_failed = True

        # 4) score
        if not validate_score(feature):
            failed_checks += 1
            feature_failed = True

        # 5) phase
        if not validate_phase(feature):
            failed_checks += 1
            feature_failed = True

        # 6) attributes
        if not validate_attributes(feature):
            failed_checks += 1
            feature_failed = True

        if feature_failed:
            failed_features += 1

    if failed_features == 0:
        logger.info(f"GFF validation passed: {total} features checked, 0 failures.")
        return True

    logger.error(
        f"GFF validation failed: {total} features checked, "
        f"{failed_features} features had errors ({failed_checks} failed checks)."
    )
    return False