# Analysis of biological data

__Differential gene expression__. Coming soon.

__Gene Ontology (GO) term enrichment__. Gene Ontology is a major bioinformatics initiative to unify the representation of gene and gene product attributes across all species. It provides a system for hierarchically classifying genes or gene products to terms organized in a graph structure called an ontology. The terms are grouped into three categories: __Molecular Function (MF)__ (describing the molecular activity of a gene), __Biological Process (BP)__ (describing the larger cellular or physiological role carried out by the gene, coordinated with other genes) and __Cellular Component (CC)__ (describing the location in the cell where the gene product executes its function). Each gene can be described (annotated) with multiple terms. Gene Ontology enrichment involves testing for a biased presence of certain GO terms population of genes (Genes of Interest (GOI)) compared to all genes (universe).

There are several programs to test for GO term enrichment. However, my program of choice is [topGO](http://bioconductor.org/packages/release/bioc/html/topGO.html) mostly because it is straightfoward to use and has good user documentation (See script `topGO.r`). I would avoid drawing too large of conclusions from GO term enrichment...
