Workflow for differential expression analysis in proteomics data

Here it is reperted a R workflow for proteomics data analysis from protein expression projects, stored on the PRIDE database and reported on the COVID-19 Data portal.
This is an R pipeline to analyze protein expression data, built on lung cell lines infected by SARS-CoV-2 variants: ​B.1, Delta, and Omicron BA.1 (Mezler et al. 2023) https://www.ebi.ac.uk/pride/archive/projects/PXD037265.
This pipeline can obtain DEPs for each variant, starting from normalized protein expression data, and it enables to obtain LogFC and FDR values highly overlapping with DEPs in Mezler et al. 2023, as well as building an input file to overlay the DEPs on C19DM.
Furthrmore, it is shared on WorkflowHub, to mark with BY-COVID and cite it on IDtk.
