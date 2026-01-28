import gffutils

def build_gff(transcript_id,db: gffutils.FeatureDB, qc_flags: str, gff_out) -> None:
    feature = db[transcript_id]
    if 'QC_flags' in feature.attributes:
        feature.attributes['QC_flags'].append(qc_flags)
    else:
        feature.attributes['QC_flags'] = [qc_flags]
    gff_out.write(str(feature) + '\n')