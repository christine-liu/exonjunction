# exonjunction

R script for visualizing exon-exon junctions of APP (amyloid precursor protein) in genomic reads. 

Needs an exon reference file (provided) (exonsAPP.bed)

Lines that contain paths that need to be edited by the user denoted by #***

Input file is in bed12 format - fastq file must be aligned using STAR option --outSAMattributes All. Resulting bam file should be filtered to only include lines that do not contain "jI:B:i,-1" (indicates no junction detected). Bam file is then converted to bed12 for use with this R script.
