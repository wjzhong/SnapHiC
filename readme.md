# SnapHiC-G: Single nucleus analysis pipeline for Hi-C data using global background
#### (Latest updates: Aug 15, 2022, version 0.1.1)
## Identifying chromatin loops from single cell Hi-C data using global background
### Contents:
1. [Installation](#installation)
2. [Required input files](#required-input-files)
3. [Optional input files](#optional-input-files)
4. [Running SnapHiC-G](#running-snaphicg)
5. [Output files](#the-output-file)
6. [Testing SnapHiC-G](#testing-snaphicg)
7. [Recommendations for parallel computing](#recommendations-for-parallel-computing)
8. [Running for different resolutions](#running-different-resolutions)
9. [More details on running SnapHiC-G](#more-detail-on-running-snaphicg)
10. [Contact us](#contact-us)
11. [Citation](#citation)

<h3 id=installation>1. Installation</h3>

You can download the singularity image/recipe of SnapHiC-G [here](TODO). If you don't have access to singularity or prefer your own installation, please follow the following steps (Note that SnapHiC-G is developed on top of the SnapHiC. SnapHiC-G requires the same python packages or modules as the SnapHiC. If you have previously installed SnapHiC, there is no need to install the required packages or modules again.):    
1. Install python version >= 3.6. 
2. Make sure MPI (e.g. [open-mpi](https://www.open-mpi.org/)) is installed on your system and is on your path (`which mpicc` can return a path). 
3. Use pip to install the required modules: 
```
pip install -r requirements.txt
```
4. Install Java version >= 8 and include it in your path. Java is required to generate .hic files for visualization in Juicebox.

<h3 id=required-input-files>2. Required input files</h3>   

**The following four input files are mandatory for SnapHiC-G.**. 

1. "Tab-separated" or "tab-separated and gzipped" files containing the mapped read pairs (contacts) for each single cell. These contact files can be generated from raw fastq files following the methods described in previous publications ([PMID: 31384045](https://pubmed.ncbi.nlm.nih.gov/31384045/), [PMID: 31501549](https://pubmed.ncbi.nlm.nih.gov/31501549/), and [PMID: 28682332](https://pubmed.ncbi.nlm.nih.gov/28682332/)), or other single cell Hi-C data preprocessing pipelines such as Dip-C (https://github.com/tanlongzhi/dip-c). In each file, one line represents one contact with 2 columns for chromosome name and 2 columns for the mapped positions (bp).   
2. chrom.sizes file for the genome of interest, which can be downloaded from [here](https://hgdownload.soe.ucsc.edu/downloads.html). Files for mm10, hg19 and hg38 are included in the `ext` directory.  
3. A binned bed file of the genomic regions that are excluded from loop calling. In our study, we defined mappability for each 10KB bin based on the restriction enzyme MboI, as described in our previous study [PMID: 23023982](https://pubmed.ncbi.nlm.nih.gov/23023982/). We removed all 10KB bins with mappability <=0.8, and all 10KB bins overlapped with the ENCODE blacklist regions. Filtered regions for mm10, hg19 and hg38 at 10KB resolution with the restriction enzyme MboI, i.e., bins with mappability <=0.8 and bins overlapped with the ENCODE blacklist regions, are included in the `ext` directory. Local genomic features (including mappability scores) for different reference genomes, different bin resolutions, and different restriction enzymes can be downloaded [here](http://enhancer.sdsc.edu/yunjiang/resources/genomic_features/).

4. A transcription start site (TSS) file for the genome of interest. Files for mm10, hg19 and hg38 downloaded from UCSC refGene annotation are included in the `ext` directory. 

<h3 id=optional-input-files>3. Optional input files</h3>

**The following two input files are optional for SnapHiC-G.**. 

1. A promoter file for genes. Files for the promoter of the ODC cell are included in the `ext` directory. If there is no input of promoter file, the default is +/-500bp of TSS as promoter.  

2. An enhancer file for genes. Files for the enhancer of the ODC cell are included in the `ext` directory. If user does not input enhancer file, SnapHiC-G keeps bin pairs that at least one bin has the promoter. If user inputs enhancer file, SnapHiC-G keeps promoter-promoter bin pairs, and promoter-enhancer bin pairs.


<h3 id=running-snaphicg>4. Running SnapHiC-G</h3>

We strongly recommend using an HPC environment where you can request multiple nodes/processors and allocate memory. Alternatively, you can still run SnapHiC-G using multithreaded and single-processor environments, provided enough memory and runtime.  

1. Put all the contact files, one for each cell, into the same directory.  

2. Edit the *run_step1.sh* and *run_step2.sh* files. If you use an HPC cluster with a job scheduler such as PBS or SLURM (we strongly recommend), specify the required nodes, processors, memory, and load the required modules (python3.6+, MPI, java8+). If you use a regular compute node without job scheduler, or a desktop computer, you can skip this step (you still need big memory and the computation can be slow).    

3. Set the following variables in the *run_step1.sh* and *run_step2.sh* files:  
&Tab;`SnapHiC_G_dir`="/path/to/directory/where/SnapHiC-G/is/located/" (path to the directory containing the *snap.py* file in SnapHiC-G).  
&Tab;`parallelism`="parallel" (it can take one of the three options: **parallel**, **threaded**, or **single-proc**. Use **parallel** if you use an HPC with job scheduler, **threaded** if you use multiple processors without job scheduler, and **singl-proc** otherwise).    
&Tab;`number_of_processors`=15 (if you use **threaded** or **parallel**, please specify the number of processors).  
&Tab;`indir`="/path/where/the/contact/files/are/stored" (contact files should be tab separated. They can be gzipped).  
&Tab;`suffix`="contacts.txt.gz" (File name suffix for the contact files, which is used to distinguish input files if the same input directory contains other files).  
&Tab;`outdir`="/path/where/the/output/files/will/be/saved".  
&Tab;`chrs`="2 4" (Two integers indicating the column numbers of the chromosomes in the contact files. Starting from 1).  
&Tab;`pos`="3 5" (Two integers indicating the column numbers of the read mapped positions in the contact files. Starting from 1).  
&Tab;`chrlen`="ext/mm10.chrom.sizes" (chrom.sizes file. You can download it from [here](https://hgdownload.soe.ucsc.edu/downloads.html)).  
&Tab;`genome`="mm10" (Name of the reference genome. SnapHiC-G currently accepts mm10, hg19 and hg38. It is used to determine the number of autosomal chromosomes, and to generate .hic file for visualization in Juicebox).   
&Tab;`filter_file`="ext/mm10_filter_regions.txt" (A binned bed file of the genomic regions that are excluded from loop calling, such as the low mappability regions and the ENCODE blacklist regions. We provide these files for mm10, hg19 and hg38 at 10KB resolution with the restriction enzyme MboI in the *ext* directory).   
&Tab;`tss_file`="ext/mm10.refGene.transcript.TSS.061421.txt" (A TSS file for the genome of interest. We provide these files for mm10, hg19 and hg38 in the *ext* directory).
&Tab;`promoter_file`="ext/ODC.promoters.txt" (A promoter file for genes).  
&Tab;`enhancer_file`="ext/ODC.enhancers.txt" (An enhancer file for genes).  
&Tab;`steps`="" (steps to run the pipeline. Recommended (1) "bin rwr" at first, (2) then  "hic interaction postprocess").   
&Tab;`prefix`="ODC" (Name of the dataset. It is the prefix of the output files, including *.postprocessed.all_candidates.bedpe*, and *allChr.[hic/cool]*. See details in [Output files](#the-output-file)). 

4. Execute the run file with the modified variables, or submit it to the job scheduler. 

<h3 id=the-output-file>5. Output files</h3> 

The SnapHiC-G-identified chromatin loops are stored in the file: *<outdir>/postprocessed/\*.postprocessed.all_candidates.bedpe*. (* is the name of the dateset, which is provided via the **--prefix** argument. If user does not provide the **--prefix** argument, * will be replaced with *combined*). This tab-separated file includes the following 11 columns:  
- chr1, x1, x2, chr2, y1, y2: the start and end position of a loop. 
- outlier_count: the number of cells with >1.96 normalized contact probability (with respect to global background) at a loop.  
- pvalue, tstat, fdr_dist: statistical measures (P-values from the one sample t-test, t-statistics from the one sample t-test, and false discovery rate for all bin pairs at the same 1D genomic distance) computed for each loop. We recommend using fdr_dist as the measure of the loop strength.
- case_avg: across all cells, the average normalized contact probability of each loop (case). 

In addition, SnapHiC-G provides a file containing all the loop candidates in the file: *<outdir>/postprocessed/\*.postprocessed.all_candidates.bedpe*. 

SnapHiC-G also generates .hic and .cool files stored at *<outdir>/hic/allChr.[hic/cool]*. SnapHiC-G first computes the % of outlier cells (i.e., the proportion of cells with normalized contact probability >1.96), and then takes the integer ceiling of 100*(% of outlier cells) to create a count matrix. SnapHiC-G then uses the Juicer software and the cooler software to convert the count matrix into .hic file and .cool file, respectively. User can skip this step by including **--no-hic** or **--no-cool** argument in the run file. 

In addition to the output files described above, SnapHiC-G generates multiple intermediate output files. Each of the five steps in the process (bin, rwr, hic, interaction and postprocess) creates a separate sub-directory in the `outdir`. The first two steps (bin and rwr) generate output files for each chromosome in each cell. The remaining three steps (hic, interaction and postprocess) combine all cells and generate one file for each chromosome. 

<h3 id=testing-snaphicg>6. Testing SnapHiC-G</h3>

To test SnapHiC-G, please download this sample dataset [here](http://renlab.sdsc.edu/abnousa/snapHiC/test/input/Ecker/ODC_100.tar.gz) (1.2GB), which contains the contact files of 100 randomly selected oligodendrocytes from Lee et al study ([PMID: 31501549](https://pubmed.ncbi.nlm.nih.gov/31501549/), Ref: hg19). The final list of XXX[TODO] loops can be downloaded from [here](TODO: odc_100.postprocessed.all_candidates.bedpe) (TODO: filesize). The complete subdirectory of output files for this sample dataset can be downloaded from [here](TODO: add the results /ODC_100). 

After downloading this sample dataset [here](http://renlab.sdsc.edu/abnousa/snapHiC/test/input/Ecker/ODC_100.tar.gz) (1.2GB), untar it to generate 100 contact files, each containing the contacts for one cell. Run files for this sample dataset are in the *test_run_scripts* directory. You can set the input and output directories, load the required modules, and submit the run files.

<h3 id=recommendations-for-parallel-computing>7. Recommendations for parallel computing</h3>

You can use as many processors as possible for the first two steps ('bin' and 'rwr'), as long as each processor has sufficient memory. For 10Kb bin resolution, we recommend 30-40GB of memory per processor. The remaining three steps ('hic', 'interaction' and 'postprocess') are performed for each chromosome separately. Therefore, we recommend using as many processors as the number of chromosomes in the genome. Providing sufficient memory for each processor for the remaining three steps is important. For 2Mb maximal genomic distance, 10Kb bin resolution, and for human and mouse genome, we recommend providing 20GB of memory per processor. We recommend using less than five processors per node. Here are the settings we used in our study.  
| genome | #cells | binsize | distance | run step | nodes | processors per node | memory per node | runtime |  
| --- | --- | --- | --- | --- | --- | --- | --- | --- |  
| mm10 | 100 | 10KB | 1Mb | 1 (bin rwr) | 15 | 3 | 96GB | ??? hrs |  
| mm10 | 100 | 10KB | 1Mb | 2 (hic inter. postproc.) | 10 | 2 | 96GB | ??? hrs |  
| mm10 | 200 | 10KB | 1Mb | 1 (bin rwr) | 15 | 3 | 96GB | ??? hrs |  
| mm10 | 200 | 10KB | 1Mb | 2 (hic inter. postproc.) | 10 | 2 | 96GB | ??? hrs |  
| mm10 | 300 | 10KB | 1Mb | 1 (bin rwr) | 15 | 3 | 96GB | ??? hrs |  
| mm10 | 300 | 10KB | 1Mb | 2 (hic inter. postproc.) | 10 | 2 | 96GB | ??? hrs |  
| mm10 | 400 | 10KB | 1Mb | 1 (bin rwr) | 15 | 3 | 96GB | ??? hrs |  
| mm10 | 400 | 10KB | 1Mb | 2 (hic inter. postproc.) | 10 | 2 | 96GB | ??? hrs |  

<h3 id=running-different-resolutions>8. Running for different resolutions </h3> 

In order to run the SnapHiC-G algorithm for other binsizes, we provide the filter regions and chromosome sizes files for 20kb, and 25kb resolutions in the "ext" directory. The user should direct to correct *filter regions* files and the *chromosome sizes* files in the run files. The user should also specify the resolution with the *binsize* argument when running the SnapHiC-G (e.g. --binsize 20000).  
 
<h3 id=more-detail-on-running-snaphicg>9. More details on running SnapHiC-G</h3>

SnapHiC-G consists of five steps: (1) bin: binning, (2) rwr: imputing contact probability via random walk with restart (RWR), (3) hic: combining all cells to generate .hic file, (4) interaction: identifying loop candidates based on the global background model, and (5) postprocess: filter loop candidates based on statistical measures. User can specify which step or steps to run as a command line argument. To specify the steps, use the argument `--step` followed by any of the five options: 'bin', 'rwr', 'hic', 'interaction', and 'postprocess'. If user don't specify any steps, SnapHiC-G will run all five steps. 

We strongly recommend running SnapHiC-G in the following two parts, which are described in the run files. 

Part 1:
```
python3 snap.py -i <input-directory> -o <output-directory> -c <column-numbers-for-chromosomes-in-data> -p <column-numbers-for-read-positions-in-data> -l <chromosome-lengths-file> -g <genome: one-of-hg/mm> [--parallel] --steps "bin rwr"
```

followed by Part 2:
```
python3 snap.py -i <input-directory> -o <output-directory> -l <chromosome-lengths-file> -g <genome: one-of-hg/mm> [--parallel] --steps "hic interaction postprocess"
```

First of all, when scHi-C dataset consists of a large number of cells, user can run the computation in parallel as much as possible by using a large number of processors per node. However, since the available memory is limited, the computation might fail due to the out of memory error (which SnapHiC-G currently is not able to detect). By waiting until the RWR step is finished, user can make sure that all chromosomes in all cells are processed. **If user restarts the RWR step, SnapHiC-G can recognize the chromosomes/cells that have already been processed and will not repeat the RWR computation for them.** 

In case you run out of memory and you suspect that some of the RWR computations might be corrupted, you can use the script (utils/validate_rwr.py) to assess whether RWR computation finished successfully in all cells. The input of the **validate_rwr.py** script is the output directory of the SnapHiC-G run. Please use command: `python utils/validate_rwr.py /path/to/output/of/your/run`. This script will generate a file named **missing.txt** in the <output>/rwr directory, which contains the name of RWR chromosomes/cells in which the RWR computation is corrupted. You can delete these files and restart the RWR step in SnapHiC-G.  

In addition, for steps 3 and 4, we recommend providing at least 20GB of memory per processor for datasets containig 100-500 cells. We also recommend not using more than 5 processors per node for step 3, since it involves simultainious access to files, and depending on your system, there might be a limit on the allowed maximum number of open files. By breaking down the computation into these two parts, user can optimally adjust the memory and processor per node according to the HPC resource and the number of cells in scHi-C dataset. 

Additional arguments that can be used, which can be found by *python snap.py --help*. 


<h3 id=contact-us>10. Contact us</h3>

For any questions, comments and suggestions regarding SnapHiC-G, please submit an issue to Github with the details of your system and run. You can also send email to Wujuan Zhong (<zhongwujuan@gmail.com>), Weifang Liu (<weifangl@live.unc.edu>), Yun Li (<yunli@med.unc.edu>) or Ming Hu (<hum@ccf.org>).

<h3 id=citation>11. Citation</h3>

If you use SnapHiC-G, please cite our pre-print at bioRxiv: TBA

TBA 
