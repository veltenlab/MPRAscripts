# MPRAscripts

Collection of snakemake pipelines and perl/Rmd scripts to process CRS-Barcode assocation libraries and Barcode counting libraries from [lentiMPRA](https://www.nature.com/articles/s41596-020-0333-5)

The pipeline consists out of 3 separate steps
1) CRS-BC association of the MPRA libraries
2) BC counting for a MPRA experiment
3) Screen report and parameter grid search


Step 1: 1_CRS_BC_association

Step 2: 2_BC_count

Step 3: 3_Report_screen



Important information for the manuscript:
We performed name changes in the manuscript compared to the naming in the processing

Processing name: Manuscript name
LibraryA_minPInitial: Library A HSC
LibraryA_minPOptimized: Library A K562
LibraryH_minP: Library B minP
LibraryH_minCMV: Library B minCMV
LibraryB_minP: Library C

ScreenID fo the different experiments:
Library A HSC -> 2_HSC_LibraryA_minPInitial
Library B HSC -> 4_HSC_LibraryH_minP
Library C HSC -> 6_HSC_LibraryB_minP

Library A K562 -> 8_K562_Transient_LibraryA_minPOptimized
Library B K562 -> 9_K562_Transient_LibraryH_minP
Library B minCMV K562 -> 10_K562_Transient_LibraryH_minCMV
Library C K562 -> 11_K562_Transient_LibraryB_minP
Library B 7d PI K562 -> 12_K562_Integrated_LibraryH_minP