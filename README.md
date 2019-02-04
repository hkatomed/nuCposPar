# nuCposPar
Codes for nuCpos parameter construction

This is a set of codes for generation of parameters used in nuCpos.


## Hardware and Software requirements ###############################
We used an Apple Mac Pro (OS 10.9.4) in which R (3.3.3) had been installed. 
To run the R scripts, R packages listed below must be loadable.

1. Biostrings (>=2.42.1)
2. TTR (>=0.23.4)

Note: final size of the nuCposPar directory will be about 2.2 GB.


## Required files ###################################################
Before running scripts, you need to obtain the files (1-4) listed below and store them in appropriate directories.

## Files that should be stored in "genome" 
1. Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz
    ftp://ftp.ensembl.org/pub/release-84/fasta/saccharomyces_cerevisiae/dna/
2. Schizosaccharomyces_pombe.ASM294v2.dna.toplevel.fa.gz
    ftp://ftp.ensemblgenomes.org/pub/fungi/release-42/fasta/schizosaccharomyces_pombe/dna/
3. Mus_musculus.NCBIM37.67.dna_rm.toplevel.fa.gz
    ftp://ftp.ensembl.org/pub/release-67/fasta/mus_musculus/dna/

## File that should be stored in "Mm_map" 
4. GSM2183909_unique.map_95pc.txt.gz
    ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2183nnn/GSM2183909/suppl/


## Parameter construction ###########################################
Run the R and sh scripts as below. 

cd nuCposPar
Rscript scripts/mm_genome.R
Rscript scripts/sc_genome.R
Rscript scripts/sp_genome.R
bash scripts/mm_unique_split.sh
Rscript scripts/mm_regions.R
Rscript scripts/mm_models.R
Rscript scripts/yeast_models.R
Rscript scripts/lHBA_models.R
Rscript scripts/model_merge.R


## Output ############################################################
sysdata.rda will contain all the parameters used in nuCpos package.

