# Color-blindess-gene-research
This is a high-level research project of OPN1LW, OPN1MW, and OPN1SW and their relationships to color blindness. Project was conducted and completed in July 2023.

# Project Steps:
1. Perform general research of each gene in UCSC Genome Browser (https://genome.ucsc.edu/cgi-bin/hgGateway). Includes obtaining chromosome position, gene coordinates, and count of protein-coding regions.
2. Obtain reference protein sequences for each gene.
3. Identify pathogenic single-nucleotide polymorphisms (SNPs) from dbSNP (https://www.ncbi.nlm.nih.gov/snp/).
4. Review 1000 Genomes Project data for allelic frequencies of pathogenic SNPs.
5. Summarize gene expression levels from GEUVADIS RNA-Sequencing Project sample files. 

# Summary of Findings:
1. Reference DNA sequence lengths:
     OPN1LW: 598 nucleotides
     OPN1MW: 2,289 nucleotides
     OPN1SW: 1,033 nucleotides
2. Based on the clinical significance categories listed in dbSNP, no pathogenic SNPs were included in the 1000 Genomes Project variant files. Only non-pathogenic/unknown clinical significance variants were analyzed.
3. The GEUVADIS RNA-Sequencing Project patient samples did not capture gene expression for OPN1MW. Gene expression values for OPN1LW were null. Only expression values for OPN1SW were analyzed.
   
![Non-Pathogenic SNPS Average Allelic Frequency](https://github.com/smbr15/Color-blindess-gene-research/blob/main/nonpath%20SNP%20freq.png)


An ANOVA test was performed to compare the average allelic frequencies of OPN1LW, OPN1MW, and OPN1SW variants. This produced a p-value of 0.3233. Given a 5% level of significance, these three average allelic frequencies are not significantly different. 


![Average OPN1SW Gene Expression by Population](https://github.com/smbr15/Color-blindess-gene-research/blob/main/opn1sw%20population%20gene%20expression.png)


![Average OPN1SW Gene Expression by Sex](https://github.com/smbr15/Color-blindess-gene-research/blob/main/opn1sw%20sex%20gene%20expression.png)


The average gene expression of OPN1SW was 1.1245 for males and 1.0849 for females. Since OPN1SW occurs on chromosome 7, it is not surprising that these are similar. It would have been interesting to see if these differed for OPN1LW and OPN1MW since they occur on chromosome X, but the data was not available. 
