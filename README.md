# E-BAG
Enumeration Based Ambigous Genotypes (HLA genes)

A python3 based script for finding potential ambiguous genotypes for HLA-typing, based on the idea of Gene Feature Enumeration (Mack SJ. 2015). HLA genotyping can be ambiguous due to the lack of phasing especially for genes like DPB1 that IMGT contains alleles sequenced only for few gene features (UTRs, Exons, Introns).  E-BAG enumerates the gene features of a requested HLA gene from IMGT and use that enumeration in order to report possible ambiguous genotypes with shuffled gene features for an input allele pair(genotype).

1) For the input HLA allele pair and the input set of gene features searches for all the alleles having the same features   enumeration
2) Calculates all allele pairs
3) Outputs only the possible ambiguous allele pairs according to the input pair enumeration for the gene features requested to
   search.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Create a directory with the .xml file from IMGT. 
Use xml2enumeration.py to create the enumeration file for a requested HLA-gene
on the command line type:
python3 scriptName xmlFile genename(for example DPB1)

Then use test4ambiguousGenotypes.py to search for possible ambiguous genotypes for your genepair and the gene features you want to examine
on the command line type:
python3 scriptName alleleName1 alleleName2 features(example:5UTR,E1,I1,E2,3UTR) inputfile outputfile
type allele names without the HLA-. example DPB1*13:01:01')
