FindSV core
===========
Run

To analyse one bam file and put the output in the output_folder type:

        python FindSV.py --bam file.bam --output output_folder

To analyse a folder containing bam files type:

        python FindSV.py --folder input_folder --output output_folder

or skip the --output flag to put the output in the directory given py the config file.

Installation
============
Conda install
prerequisites
    Conda

        use the command:  
                python FindSV.py --install

        choose the install mode that suits your application the best. For exmple:
                python FindSV.py --install UPPMAX

to run FindSV-core on the uppmax system.

Restart module
============
for info on how to restart samples analysed by the FIndSV-core pipeline, type:

        python FindSV.py --restart


Settings
=========

        Reference_dir
            The path to the reference directory is set using the reference_dir flag.
            The reference needs to be split per chromosome.

Other settings
=============

        general:
            TMPDIR: path to the directory used as scratch drive(not implemented yet)

            account: the slurm account

            output: optional defult output folder


        calling:

            FindTranslocation:
  
                minimum supporting pairs: the minimum number o pairs to call a variant.
    
    
         CNVnator
  
            bin size: the base pair sie of each bin used to search for CNVs


        annotation:

            internal frequency DB
  
                the minimum overlap to count a variant as a hit in the DB
    
            Genmod
  
                rank model path
    
                    the path to the genmod rankmodel
      
            VEP
  
                cache_dir
    
                    the same as the vep --dir option; this option is required on uppmax


