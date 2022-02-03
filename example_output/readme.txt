The file "odc_100.postprocessed.all_candidates.bedpe" is the example output file from the 5th step ("postprocess") of SnapHiC-G. 

Because the first three steps ("bin", "rwr" and "hic") of SnapHiC-G are the same as SnapHiC, this output file ("odc_100.postprocessed.all_candidates.bedpe") can be generated using the script "example_interaction_postprocess.sh" as below:

1. download and unzip the file "hic.tar.gz" from the SnapHiC output at the link: http://renlab.sdsc.edu/abnousa/snapHiC/test/output/Ecker/ODC_100/
2. put the "hic" folder under the directory specified in the variable 'indir' of the script
3. modify other variables such as 'snapHiC_G_dir', 'suffix', 'outdir', 'chrlen', 'filter_file', 'tss_file' and other settings to successfully run the script

Note that the script "example_interaction_postprocess.sh" only runs the 4th and 5th steps of the SnapHiC-G method. For your own dataset, you need to run the full five steps ("bin", "rwr", "hic", "interaction", and "postprocess") to get the SnapHiC-G results. 
