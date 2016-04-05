import argparse
import sys
import yaml
import os
import FindSV_modules
import subprocess
import re
import fnmatch
import tracking_module
import setup
import test_modules

#the module used to perform the variant calling
def calling(args,config,output,scripts,programDirectory):
    prefix=args.prefix

    #run the callers
    input_vcf=""
    sbatch_ID=[]
    for caller in config["FindSV"]["calling"]:
        general_config=config["FindSV"]["general"]
        caller_config=config["FindSV"]["calling"][caller]
        #run the FT script
        if caller == "FT":
            job_name="FT_{}".format(prefix)
            process_files=os.path.join(output,"slurm/calling/",job_name)
            #generate the header
            FT=scripts["FindSV"]["header"].format(account=general_config["account"],time="30:00:00",name=job_name,filename=process_files)
            output_prefix=os.path.join(output,"{}_FT".format(prefix))
            #generate the body
            FT += scripts["FindSV"]["calling"][caller].format(output=output_prefix,FT_path=caller_config["FT_path"],bam_path=args.bam,minimum_suporting=caller_config["minimum_supporting_pairs"])
            sbatch_ID.append(submitSlurmJob( os.path.join(output,"slurm/calling/FT_{}.slurm".format(prefix)),FT) )
            input_vcf += "{}_inter_chr_events.vcf ".format(output_prefix)
            input_vcf += "{}_intra_chr_events.vcf ".format(output_prefix)
        #run cnvnator
        elif caller =="CNVnator":
            job_name="CNVnator_{}".format(prefix)
            process_files=os.path.join(output,"slurm/calling/",job_name)
            CNVNator=scripts["FindSV"]["header"].format(account=general_config["account"],time="30:00:00",name=job_name,filename=process_files)
            output_prefix=os.path.join(output,prefix)
            #if the user want to use uppmax settings, load CNVNator module, otherwise load rootsys to path, if none is given, assume that rootsys is permanently added to path
            if not general_config["UPPMAX"] == "":
                CNVNator +=scripts["FindSV"]["UPPMAX"].format(modules="bioinfo-tools CNVnator")
                caller_config["CNVnator_path"]="cnvnator"
                caller_config["CNVnator2vcf_path"]="cnvnator2VCF.pl"
            elif not caller_config["ROOTSYS"] =="":
                CNVNator +=scripts["FindSV"]["ROOTSYS"].format( rootdir=caller_config["ROOTSYS"] )
            CNVNator += scripts["FindSV"]["calling"][caller].format(output=output_prefix,CNVnator_path=caller_config["CNVnator_path"],bam_path=args.bam,bin_size=caller_config["bin_size"],reference_dir=caller_config["reference_dir"],CNVnator2vcf_path=caller_config["CNVnator2vcf_path"])
            input_vcf += "{}_CNVnator.vcf ".format(output_prefix)
            sbatch_ID.append(submitSlurmJob( os.path.join(output,"slurm/calling/CNVnator_{}.slurm".format(prefix)) ,CNVNator) )

    return(input_vcf,sbatch_ID)

#the module used to perform the combining of multiple callers
def combine_module(args,config,output,scripts,programDirectory,input_vcf,sbatch_ID):
    general_config=config["FindSV"]["general"]
    #combine module; combine all the caller modules into one VCF
    prefix=args.prefix
    output_prefix=os.path.join(output,prefix)
    job_name="combine_{}".format(prefix)
    process_files=os.path.join(output,"slurm/combine/",job_name)
    merge_VCF_path=os.path.join(programDirectory,"internal_scripts","mergeVCF.py")
    contig_sort=os.path.join(programDirectory,"internal_scripts","contigSort.py")
    combine=scripts["FindSV"]["header"].format(account=general_config["account"],time="3:00:00",name=job_name,filename=process_files)
    #if we are on Uppmax, the samtools module is loaded, otherise it is assumed to be correctly installed
    combine += scripts["FindSV"]["afterok"].format(slurm_IDs=":".join(sbatch_ID))
    if not general_config["UPPMAX"] == "":
        combine +=scripts["FindSV"]["UPPMAX"].format(modules="bioinfo-tools samtools")
    #if we are not on uppmax and the samtools conda module is installed
    elif not config["FindSV"]["conda"]["samtools"] == "":
        combine +=scripts["FindSV"]["conda"].format(environment="samtools_FINDSV")
    outputVCF=output_prefix+"_FindSV.vcf"
    combine += scripts["FindSV"]["combine"]["combine"].format(output=output_prefix,merge_vcf_path=merge_VCF_path,input_vcf=input_vcf,contig_sort_path=contig_sort,bam_path=args.bam,output_vcf=outputVCF)
    combine_ID=submitSlurmJob( os.path.join(output,"slurm/combine/combine_{}.slurm".format(prefix)) , combine)
    
    return(outputVCF,combine_ID)

#the module used to perform the annotation
def annotation(args,config,output,scripts,programDirectory,outputVCF,combine_ID):
    general_config=config["FindSV"]["general"]
    prefix=args.prefix
    output_prefix=os.path.join(output,prefix)
    #annotation module; filter and annotate the samples
    annotation_config=config["FindSV"]["annotation"]
    job_name="annotation_{}".format(prefix)
    process_files=os.path.join(output,"slurm/annotation/",job_name)
    annotation = scripts["FindSV"]["header"].format(account=general_config["account"],time="10:00:00",name=job_name,filename=process_files)
    annotation += scripts["FindSV"]["annotation"]["header"].format(combine_script_id=combine_ID)
    
    #if uppmax modules are chosen, the vep module is loaded
    if not general_config["UPPMAX"] == "":
        annotation +=scripts["FindSV"]["UPPMAX"].format(modules="bioinfo-tools vep")
    #otherwise if the vep conda environment is installed, we will use it
    elif not config["FindSV"]["conda"]["vep"] == "":
        annotation +=scripts["FindSV"]["conda"].format(environment="VEP_FINDSV")
    
    #add frequency database annotation
    if not annotation_config["DB"]["DB_script_path"] == "" and not annotation_config["DB"]["DB_path"] == "":
        inputVCF=outputVCF
        outputVCF=output_prefix+"_frequency.vcf"
        annotation += scripts["FindSV"]["annotation"]["DB"].format(query_script=annotation_config["DB"]["DB_script_path"],output=output_prefix,db_folder_path=annotation_config["DB"]["DB_path"],input_vcf=inputVCF,output_vcf=outputVCF)
    cache_dir=""
    
    #add vep annotation
    if not annotation_config["VEP"]["cache_dir"] == "":
        cache_dir=" --dir {}".format(annotation_config["VEP"]["cache_dir"])
    
    #DO not use local system vep if uppmax or conda is chosen
    if not general_config["UPPMAX"] == "" or not config["FindSV"]["conda"]["vep"] == "":
        inputVCF=outputVCF
        outputVCF=output_prefix+"_vep.vcf"
        annotation += scripts["FindSV"]["annotation"]["UPPMAX_VEP"].format(vep_path=annotation_config["VEP"]["VEP.pl_path"],output=output_prefix,port=annotation_config["VEP"]["port"],cache_dir=cache_dir,input_vcf=inputVCF,output_vcf=outputVCF)
    #if we do not use uppmax or conda and a path to the vep script is added in the config, then use that vep script(otherwise skip vep annotation)
    elif not annotation_config["VEP"]["VEP.pl_path"] == "":
        inputVCF=outputVCF
        outputVCF=output_prefix+"_vep.vcf"
        annotation += scripts["FindSV"]["annotation"]["VEP"].format(vep_path=annotation_config["VEP"]["VEP.pl_path"],output=output_prefix,port=annotation_config["VEP"]["port"],cache_dir=cache_dir,input_vcf=inputVCF,output_vcf=outputVCF)

    #add genmod annotation
    if not annotation_config["GENMOD"]["GENMOD_rank_model_path"] == "":
        #use the genmod conda module if the user wishes to do so
        if not config["FindSV"]["conda"]["samtools"] == "":
            annotation +=scripts["FindSV"]["conda"].format(environment="GENMOD_FINDSV")
        inputVCF=outputVCF
        outputVCF=output_prefix+"_genmod.vcf"
        annotation += scripts["FindSV"]["annotation"]["GENMOD"].format(genmod_score_path=annotation_config["GENMOD"]["GENMOD_rank_model_path"],output=output_prefix,input_vcf=inputVCF,output_vcf=outputVCF)
    #create a final cleaned vcf
    inputVCF=outputVCF
    outputVCF=output_prefix+"_cleaned.vcf"
    clean_VCF_path=os.path.join(programDirectory,"internal_scripts","cleanVCF.py")
    annotation += scripts["FindSV"]["annotation"]["cleaning"].format(output=output_prefix,VCFTOOLS_path=clean_VCF_path,input_vcf=inputVCF,output_vcf=outputVCF)
    return(outputVCF,submitSlurmJob( os.path.join(output,"slurm/annotation/annotation_{}.slurm".format(prefix)) , annotation))
    
    



#this function prints the scripts, submits the slurm job, and then returns the jobid
def submitSlurmJob(path,message):
    slurm=open( path ,"w")
    slurm.write(message)
    slurm.close()

    process = "sbatch {0}".format(path)
    p_handle = subprocess.Popen(process, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)    
    p_out, p_err = p_handle.communicate()
    try:
        return( re.match(r'Submitted batch job (\d+)', p_out).groups()[0] );
    except:
        return("123456")

#read the config file, prefer command line option to config
def readconfig(path,command_line):
    config={}
    if path:
        config_file=path
    else:
        programDirectory = os.path.dirname(os.path.abspath(__file__))
        config_file=os.path.join(programDirectory,"config.txt")
    with open(config_file, 'r') as stream:
        config=yaml.load(stream)
    return(config)

def main(args):
    programDirectory = os.path.dirname(os.path.abspath(__file__))
    config=readconfig(args.config,args);
    caller_slurm_ID=[];
    caller_output=[];
    #fetch the scripts
    scripts=FindSV_modules.main()
    #initiate the output location
    if(args.output):
        output=args.output
    else:
        output=config["FindSV"]["general"]["output"]
    
    tracking=True;
    if config["FindSV"]["general"]["tracking"] == "" or args.no_tracking:
        tracking=False
    #create a folder to keep the output sbatch scripts and logs
    if not os.path.exists(os.path.join(output,"slurm/calling/")):
        os.makedirs(os.path.join(output,"slurm/calling/"))
    if not os.path.exists(os.path.join(output,"slurm/combine/")):
        os.makedirs(os.path.join(output,"slurm/combine/"))
    if not os.path.exists(os.path.join(output,"slurm/annotation/")):
        os.makedirs(os.path.join(output,"slurm/annotation/"))
    if not os.path.exists(os.path.join(output,"tracker.yml")) or not tracking:
        tracking_module.generate_tracker(output)
    #the prefix of the output is set to the prefix of the bam file
    prefix=args.bam.split("/")[-1]
    args.prefix=prefix.replace(".bam","")
    
    
    with open(os.path.join(output,"tracker.yml"), 'r') as stream:
        tracker=yaml.load(stream)
    
    #run the callers
    if not args.prefix in tracker["FindSV"]["CNVnator"] and not args.prefix in tracker["FindSV"]["FindTranslocations"]:    
        outputVCF,sbatch_ID=calling(args,config,output,scripts,programDirectory)
        caller_output=outputVCF.split()

        tracking_module.add_sample(args.prefix,args.bam,[caller_output[0],caller_output[1]],sbatch_ID[0],"FindTranslocations",output)
        tracking_module.add_sample(args.prefix,args.bam,[caller_output[2]],sbatch_ID[1],"CNVnator",output)
        caller_vcf=outputVCF
    else:
        tracker=tracking_module.update_status(args.prefix,"CNVnator",output)
        tracker=tracking_module.update_status(args.prefix,"FindTranslocations",output)
        
    #combine them
    if not args.prefix in tracker["FindSV"]["combine"]: 
        outputVCF,combine_ID=combine_module(args,config,output,scripts,programDirectory,caller_vcf,sbatch_ID)
        tracking_module.add_sample(args.prefix,caller_vcf,[outputVCF],combine_ID,"combine",output)
        combine_vcf=outputVCF
        
    else:
        tracker=tracking_module.update_status(args.prefix,"combine",output)
    
    #annotate the vcf  
    if not args.prefix in tracker["FindSV"]["annotation"]:
        outputVCF,annotation_ID=annotation(args,config,output,scripts,programDirectory,combine_vcf,combine_ID)   
        tracking_module.add_sample(args.prefix,combine_vcf,[outputVCF],annotation_ID,"annotation",output)
    else:
        tracker=tracking_module.update_status(args.prefix,"annotation",output)


parser = argparse.ArgumentParser("FindSV core module",add_help=False)
parser.add_argument('--bam', type=str,help="analyse the bam file using FindSV")
parser.add_argument("--folder", type=str,help="analyse every bam file within a folder using FindSV")
parser.add_argument('--output', type=str,default=None,help="the output is stored in this folder")
parser.add_argument("--no_tracking",action="store_true",help="run all input samples, even if they already have been analysed")
parser.add_argument("--update_tracker",action='store_true',help="update the tracker of one of a selected output folder(default output if none is chosen)") 
parser.add_argument("--config",type=str, default=None,help="the location of the config file(default= the same folder as the FindSV-core script")
parser.add_argument("--test",action="store_true",help="Check so that all the required files are accessible")
parser.add_argument("--install",action="store_true",help="Install the FindSV pipeline")
parser.add_argument("--restart",action="store_true",help="restart module: perform the selected restart on the specified folder(default output if none is chosen)")
args, unknown = parser.parse_known_args()

programDirectory = os.path.dirname(os.path.abspath(__file__))
#test to see if all components are setup
if args.test:
    print("Testing the pipeline components")
    if args.config:
        config_path=args.config
    else:
        config_path=os.path.join(programDirectory,"config.txt")
    caller_error,annotation_error=test_modules.main(config_path)
    print("-----results-----")
    if caller_error:
        print("ERROR, the callers are not properly setup")
    else:
        print("the callers are ok!")
    if annotation_error:
        print("ERROR, the annotation tools are not properly setup")
    else:
        print("the annotation tools are ok!")
    if caller_error or annotation_error:
        print("all the errors must be fixed before running FindSV. In order to get the best posible results, the warnings should be fixed as well")
#install the pipeline
elif args.install:
    
    parser = argparse.ArgumentParser("FindSV core module:Install module")
    parser.add_argument("--auto",action="store_true",help="install all required software automaticaly")
    parser.add_argument("--manual",action="store_true",help="the config file is generated, the user have to set each option manually")
    parser.add_argument("--conda",action="store_true",help="Install conda modules, the user needs to install the other software manually as well as to set the path correctly in the config file")
    parser.add_argument("--UPPMAX",action="store_true",help="set the pipeline to run on UPPMAX, install all the required software")
    args, unknown = parser.parse_known_args()
    if not os.path.exists(os.path.join(programDirectory,"config.txt")):
        setup.generate_config(programDirectory)
        if args.UPPMAX:
            setup.UPPMAX(programDirectory)
        elif args.conda:
            setup.conda(programDirectory)
        elif args.auto:
            setup.auto(programDirectory)
    else:
        print("warning: a config file is already installed, delete it or move before generating another one")
#update the status of analysed files   
elif args.update_tracker:
    parser = argparse.ArgumentParser("FindSV core module:tracker module")
    parser.add_argument("--update_tracker",type=str,nargs="*",help="update the tracker of one of a selected output folder(default output if none is chosen)")
    parser.add_argument("--config",type=str, default=None,help="the location of the config file(default= the same folder as the FindSV-core script")
    args, unknown = parser.parse_known_args()
    if not args.update_tracker:
        config=readconfig(args.config,args);
        args.update_tracker= [config["FindSV"]["general"]["output"]]
    for tracker in args.update_tracker:
        tracking_module.update_tracker(tracker)
 
#analyse one single bam file   
elif args.bam:
    caller_error=False;annotation_error=False
    try:
        if args.config:
            config_path=args.config
        else:
            config_path=os.path.join(programDirectory,"config.txt")
        caller_error,annotation_error=test_modules.main(config_path)
    except:
        pass
        
    if caller_error or annotation_error:
        print("FindSV is not correctly setup, all errors must be solved before running")
    else:
        if os.path.exists(os.path.join(programDirectory,"config.txt")):
            main(args)
        else:
            print("use the install module to generate a config file before running the pipeline")
#analyse all bamfiles within a folder(recursive searching)
elif args.folder:
    caller_error=False;annotation_error=False
    try:
        if args.config:
            config_path=args.config
        else:
            config_path=os.path.join(programDirectory,"config.txt")
        caller_error,annotation_error=test_modules.main(config_path)
    except:
        pass
        
    if caller_error or annotation_error:
        print("FindSV is not correctly setup, all errors must be solved before running")
    else:
        if os.path.exists(os.path.join(programDirectory,"config.txt")):
            for root, dirnames, filenames in os.walk(args.folder):
                for filename in fnmatch.filter(filenames, '*.bam'):
                    bam_file=os.path.join(root, filename)
                    args.bam=bam_file
                    main(args)
        else:
            print("use the install module to generate a config file before running the pipeline")

#the restart module
elif args.restart:
    parser = argparse.ArgumentParser("FindSV core module:restart module")
    parser.add_argument("--failed",action="store_true",help="restart all failed samples")
    parser.add_argument("--cancelled",action="store_true",help="restart all cancelled samples")
    parser.add_argument("--combine",action="store_true",help="restart all samples within this tracker to the combine step")
    parser.add_argument("--annotation",action="store_true",help="reruns the annotation step on all samples")
    parser.add_argument("--full",action="store_true",help="restarts the analysis from scratch")
    parser.add_argument("--restart",type=str,nargs="*",help="restart module: perform the selected restart on the specified folder(default output if none is chosen)")
    parser.add_argument("--config",type=str, default=None,help="the location of the config file(default= the same folder as the FindSV-core script")
    args, unknown = parser.parse_known_args()
    if not args.restart:
        config=readconfig(args.config,args);
        args.update_tracker= [config["FindSV"]["general"]["output"]]
    for tracker in args.restart:
        tracking_module.restart(tracker, args)
    
else:
    parser.print_help()
