'''
Part 1 - parse the GFF - done 
Data checker - done needs implemented 
Fasta parsers - saved as a dict 

DT PLAN;;

this is what i plan to do for this today and tomorrow::

LOOK INTO JINJA2 FOR HTML REPORTING

1. add details at the top in a bullet point 
tool name + version
timestamp
command used / parameters (especially QC thresholds)
input GFF3 path + hash (if stored)
links to:
../transcript_summary.tsv
../run.json
../qc_flags.gff3 or ../qc_flags.bed (if exists)

2. extract tsv information and store as bulletpoints:
total genes (unique gene_ids)
total rows (total transcripts in the gff)
transcripts per gene - mean/median/max
percentage with has_cds=true 
percentage QC flags (flags not null) 
top 3 most common flags 

3. use matplotlib to do the following:
1. store the data for transcript per gene distribution from n_transcripts_per_gene by grouping unique gene_Ids and then showing as a bar count
2. histogram of n_exons acrtoss transcripts
3. counts per QC flag type - how many of each QC issue did we see
4. side by side bar comparisons of flagged QC v unflagged QC 

4. table defining flag definitions (what does each flag mean?  





3. Grab the transcript and gene_ids

Gene lines set as type=gene with the gene_Id. 
Parse transcript lines to grab mRNA/transcript with the ID 
Primary transcript becomes a combination (ask hans/john any other suggestions) 
Exon/cds features attach via the parent class to transcripts 

4. Calculations for the output file and QC - have we found any errors when calculating 

gene_id: from Parent on the transcript line
n_exons: from exon_count function 
has_cds: from has_cds function 
flags: grab from the flags_by_trancript function and resolve them 

5. Design the outputs in main()
.json? - speak with john/hans about the best way to do this
tsv should be easy enough inc the qc_flags.gff3
3 pillars today - extension to look nice. run.json-define in python as a dict
how to write values of a dict to a json file. 


###########################################################
Testing logic?:
part 1 - parsing testing 
make sure the parser skips hashtag lines 
check strand is only +, - or . 
check start and end are integers and start < end
check seqid, source, type are strings
check score is float or NA
check phase is 0,1,2 or NA  
parsing the file should be easy enough using our old dataset
need to make sure each column in the gff is where it shoud be and correct data type
check each row is the correct length (9 columns)
check each data point is correct and normalise NA values
warn the user if not and store data for QC flags 
check for any missing data in the relationships
make sure each transcript links to a gene 
Parent=tx1;gene_id=geneA;ID=tx1 - make sure each row looks like this and warn user if it doesnt 
make sure to remove whitespace to avoid empty strings/crashing 
check for duplicates?
check there are no strange characters in the IDs
might want to check all the data is the same casing for consisrtency (lowercase everything?)

part 2: transcript --> gene testing 
check to make sure multiple transcripts link back to the same gene_id
could we check logically that each transcript has exons/cds?
check to make sure we dont have duplicate transcript ID appears twice with different gene IDs

part 3: Exon counting tests
check we are correctly counting exons/cds correctly 
Exon with multiple parents (Parent=tx1,tx2) - make sure both transcripts get the exon counted  
Two exons for the same transcript overlap in coordinates should count as 2 

part 4: has_cds tests
One CDS line with Parent=tx1 - expect has_cds=True test to return True 
no cds remains False 
CDS line with multiple parents - flag both lines to be true 

part 5: QC flag tests
missing gene_id in attributes - expect flag recorded for missing gene_id 

###### Extras on top (L4) ######
Part 6 - integrate FASTA file to get sequence lengths
Pushing to Bioconda
Push to nextflow and underastand how to use in pipeline
Make sure all dependencies are in place
Create environment.yml for conda
Make sure to test conda package locally before pushing
Set up bioconda recipe and submit PR
'''