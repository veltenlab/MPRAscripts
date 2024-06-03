# MPRAscripts
 
Snakemake pipeline for processing [lentiMPRA](https://www.nature.com/articles/s41596-020-0333-5) data. The scripts run different perl/Rmd scripts to process CRS-BC association libraries and perform barcode counting.

## The pipeline consists out of 3 separate steps
- **[CRS-BC association](1_CRS_BC_association)**: Creates dictionary of CRS-BS association
- **[BC counting](2_BC_count)**: Counts the amount of barcodes for RNA/DNA samples of an MPRA experiment
- **[Parameter search](3_Report_screen/scripts)**: Performs a parameter grid search to estimate the best processing parameters for the MPRA experiment


We processed lentiMPRA for three different GRE libraries. Sequences and design information for each of them is linked [here](Library_design/).<br>
To illustrate the input for each step of the pipeline, we process [Library A](Library_design/design_files/LibASeqDesign.fa) as an example. 

__Important information for the processing of our libraries:__ <br> 
We have a common way on how our files/fastq's are called, for example `1MDNA1` or `3ERNA2`. The first 2 letters describe the cell type, afterwards its stated if it is DNA or RNA and the last letter depicts the replicate number.

## Step 1: 1_CRS_BC_association

In the [config](1_CRS_BC_association/config.yml) we select the different parameters for the processing of the association. 

```BASH
Design:
  [LibASeqDesign] 
# All options: LibASeqDesign, LibHSeqDesign, LibBSeqDesign

FastqDir:
  [LibraryA_minPInitial] 
# All options: LibraryA_minPInitial, LibraryA_minPOptimized, LibraryH_minP, LibraryH_minCMV, LibraryB_minP

LibraryName:
  [LibAminPInitial] 
# All options: LibAminPInitial, LibAminPOpt, LibraryH_minP, LibBminP

#NUmber of threads for bwa
threads_bwa: 8

# here enter path to directory with fastq files
InDir: "/path/to/crs-bc-fastq/folder" 

# output directory of the CRS-BC association
OutDir: "/path/to/project/folder/MPRAscripts/1_CRS_BC_association/output_snake" 

#  Which design folders to use
DesignFolder: "/path/to/project/folder/MPRAscripts/Library_design" 
```
<br>

Afterwards the script indexes the design files `bwa index {input.input}`, <br>
maps the reads `bwa mem -a -t {params.threads} {input.design} {input.FWD} {input.REV} | samtools view -b > {output.bam}`, <br>
links CRS with BC `perl scripts/MapCrsBCBamFilter.pl {input.BC} {input.bam} {output}` <br>
and creates a final report for this CRS-BC association `CrsBCReport.Rmd`.

Off note: Since LibB contains identical sequences (differ only in being reverse complement) we needed to modify the [perl script](1_CRS_BC_association/LibB).

## Step 2: 2_BC_count

In this [config](2_BC_count/config.yml) file me specify different directories, as well the parameters on how the barcodes for a MPRA experiment shall be counted and if down-sampling should be performed.

```BASH
#Directories
AssDir: "/path/to/project/folder/MPRAscripts/1_CRS_BC_association/output_snake" 
InDir: "/path/to/mpra-fastq/folder"
OutDir: "/path/to/project/folder/MPRAscripts/2_BC_count/output_snake"

# Parameters for Libraries
Design:
  [LibASeqDesign] 
# All options: LibASeqDesign, LibHSeqDesign, LibBSeqDesign

ScreenID:
  [2_HSC] 
# HSC: 2_HSC (LibA), 4_HSC (LibB), 7_HSC (LibC) 
# K562: 8_K562_Transient (LibA), 9_K562_Transient (LibB), 10_K562_Transient (LibB minCMV), 11_K562_Transient (LibC), 12_K562_Integrated (LibB 7d PI)

FastqDir:
  [LibraryA_minPInitial] 
# All options: LibraryA_minPInitial, LibraryA_minPOptimized, LibraryH_minP, LibraryH_minCMV, LibraryB_minP

LibraryName:
  [LibAminPInitial] 
# All options: LibAminPInitial, LibAminPOpt, LibraryH_minP, LibBminP

# Parameters for count_BC.pl
BCcounting:
  MINRUMI: 20 # Generate a  entry for how many BCs have this much UMIs fron 1 to the number
  MAXRUMI: 10000000 # Set high treshhold in case some sequence is highly covered
  downs: 0 # 0 = no downsampling or 1 = downsampling
  downsto: 0.5 # percentage of data to downsample to
```

Core of this step is [count_BC.pl](2_BC_count/scripts/count_BC.pl), which does the following steps in short:

<pre>
#0. Read in the filtered assignments from the previous step
#1. Check that the fwd and the rev. BC read ar the same
#2. Check if it matches any of the assigned BCs perfectly.
#3. (optional) if not, use kmers to rapidly find a match and then match to confirm - currently count these to see how big the problem is.
     #the problem here is that even with an editing distance of 1 there can be several matches etc. require an extact match, but use the larger barcode list.
#4. Remember the UMI that this BCs has and the number of reads for that UMI
#5. Merge similar UMIs based (remove more lowly covered UMIs)
#6. Report statistics - how covered is each UMI (sequencing saturation)
#7. Output number of UMIs per BC, assignment
</pre>


## Step 3: 3_Report_screen

In the last step we perform a parameter grid search to evaluate the best performance.
This is done for each MPRA experiment individually, so this the [Snakemake pipeline](3_Report_screen/scripts/2_HSC_LibraryA_minPInitial) for Library A.

The [config file](3_Report_screen/scripts/2_HSC_LibraryA_minPInitial/config.yml) is split into 2 components
The first components consists out of the following:

```BASH
InDir: "/path/to/project/folder/MPRAscripts/2_BC_count/output_snake"
OutDir: "/path/to/project/folder/MPRAscripts/3_Report_screen/output_snake"

#Screen name
Screen: [2_HSC_LibraryA_minPInitial]

InitialScreen:
  cores: 7
  bccrsass: "/path/to/project/folder/MPRAscripts/1_CRS_BC_association/output_snake/LibraryA_minPInitial/LibASeqDesign_LibAminPInitial_mapped_filtered.csv.gz"
  normalizetoinput: FALSE
```

With this we create `experiment_matrix.rds` by running the [following script](3_Report_screen/scripts/processScreen.R). The script needs to know which CRS-BC association to use and if normalization of the MPRA data shall be performed with the DNA samples of the experiment or with the DNA from the CRS-BC association. The final `experiment_matrix.rds` consists out of a list with an entry for each cell type. <br>

In the second step we specifiy the parameters for the grid search. <br>

__These are the key parameters for the grid search:__
<pre>
  pseudocount: [TRUE, FALSE] # If pseudocounting shall be performed
  usecol: ["UMI1","UMI2",...] #Number of UMI´s per GRE to use
  crsfilter: [25,30,...] # How many unqiue barcodes need to be there in order to pass filter
  AssReadsFilter: [3,4,...] # Number of reads in the CRS-BC association file
  dnacrsisumi: ["Non-Binarize","Binarize"] # If DNA events should be binarized and seen as only qualitative and not quantitative measure.
</pre>

And these are the total settings for the second part:

```BASH
Screenvar:
  pseudocount: [TRUE, FALSE] 
  usecol: ["UMI1","UMI2","UMI3","UMI4","UMI5"] 
  crsfilter: [10,15,20,25,30,35,40,45,50,60,70,80,90,100] 
  AssReadsFilter: [1,2,3,4,5] 
  dnacrsisumi: ["Non-Binarize","Binarize"] 
  mindnaumi: [1]
  minrnaumi: [0]
  pseudocounttype: ["DNA"] 
  pseudocountvalue: [1]
  medupper: [3]
  medlower: [2]


Screenparams:
  cores: 7
  returnexperiments: FALSE
  bccrsass: "/path/to/project/folder/MPRAscripts/1_CRS_BC_association/output_snake/LibraryA_minPInitial/LibASeqDesign_LibAminPInitial_mapped_filtered.csv.gz"
  metainfo: "/path/to/project/folder/MPRAscripts/Library_design/meta_files/LibASeqDesign_meta_info.csv"
  metacontrolsinfo: "/path/to/project/folder/MPRAscripts/Library_design/meta_files/LibASeqDesign_controls_meta_info.csv"
  normalizetoinput: FALSE
  rnacrsisumi: FALSE
  barcodeWhitelist: FALSE
  homopolymer: 0
  filtercontrolGRE: FALSE
  plot: TRUE
  returnreport: TRUE
  savereport: TRUE
  filtercolumns: TRUE 
  Trp53Lib: "/path/to/project/folder/MPRAscripts/Library_design/meta_files/LibASeqDesign_meta_info.csv"


Screenloop:
  pseudocount: [TRUE, FALSE] 
  crsfilter:
    TRUE: [10,15,20,25,30,35,40,45,50,60,70,80,90,100]
    FALSE: [10,15,20,25,30,35,40,45,50,60,70,80,90,100]
```

The [script](3_Report_screen/scripts/processAll.R) will create a folder for each unique combination of parameters.
The script will also load a sequecing statistic of the MPRA experiment, so for the example of Library A its stored [here](2_BC_count/data). 

## Important information of the processing compared to the manuscript
We performed name changes in the manuscript (left column) compared to the naming in the processing (right column). <br>
<pre>
Library A HSC -> LibraryA_minPInitial
Library A K562 -> LibraryA_minPOptimized
Library B minP -> LibraryH_minP
Library B minCMV -> LibraryH_minCMV
Library C -> LibraryB_minP
</pre>

ScreenID fo the different experiments:
<pre>
Library A HSC -> 2_HSC_LibraryA_minPInitial
Library B HSC -> 4_HSC_LibraryH_minP
Library C HSC -> 6_HSC_LibraryB_minP

Library A K562 -> 8_K562_Transient_LibraryA_minPOptimized
Library B K562 -> 9_K562_Transient_LibraryH_minP
Library B minCMV K562 -> 10_K562_Transient_LibraryH_minCMV
Library C K562 -> 11_K562_Transient_LibraryB_minP
Library B 7d PI K562 -> 12_K562_Integrated_LibraryH_minP
</pre>

If you have any questions at all, do not hesitate to write to us at lars.velten@crg.eu.

Please cite our manuscript if you are using our scripts/data
 <br>
 
