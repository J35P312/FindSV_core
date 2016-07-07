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
prerequisites: Conda, python 2.7.11

The FindSV pipeline is setup using the install command. This command generates a config file. FindSV will not allow the user to start the analysis before the installation is finished.

To install FindSV, use the following command:

                python FindSV.py --install

choose the install mode that suits your application the best. For exmple:

                python FindSV.py --install --UPPMAX

this will setup FindSV to run on the uppmax system. FindSV may be installed using one out of four distinct modes:

        UPPMAX:
                The most convenient setting, this mode will configure FindSV to run on the UPPMAX cluster
                python FindSV.py --install --UPPMAX
                
        auto:
                Installs everything, this mode may be run in three sub settings, to run CNVnator the ROOT library must be setup. To do so, the auto option is run using ne out of three different settings:
                standard:
                        tries to match the server OS with a precompiled version of ROOT.
                        python FindSV.py --install --auto
                no_root
                        skips the installation of root, FindSV assumes root is installed according to the ROOT manual:
                        python FindSV.py --install --auto --no_root
                
                compile_root
                        If FindSV is unnable to find a suitable pre compiled root package, and the user has not installed root already. FindSV will download the source code of root and compile it.
                        python FindSV.py --install --auto --compile_root
                
        conda:
                Installs the conda environments required by FindSV, the pah to callers and other tools is added manually
                python FindSV.py --install --conda
        
        manual:
                nohing is installed, an empty config file is generated, the user must configure this file manually.
                python FindSV.py --install --manual


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


