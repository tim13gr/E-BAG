# E-BAG
Enumeration Based Ambigous Genotypes (HLA genes)

A python3 based script for finding potential ambiguous genotypes for HLA-typing, based on the idea of Gene Feature Enumeration (Mack SJ. 2015). HLA genotyping can be ambiguous due to the lack of phasing especially for genes like DPB1 that IMGT contains alleles sequenced only for few gene features (UTRs, Exons, Introns).  E-BAG enumerates the gene features of a requested HLA gene from IMGT and use that enumeration in order to report possible ambiguous genotypes with shuffled gene features for an input allele pair(genotype).

1) For the input HLA allele pair and the input set of gene features searches for all the alleles having the same features   enumeration
2) Calculates all allele pairs
3) Outputs only the possible ambiguous allele pairs according to the input pair enumeration for the gene features requested to
   search.
