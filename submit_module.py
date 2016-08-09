import os
import FindSV_modules
import tracking_module
import subprocess
import re

def run_callers(tracker,args,output,config,account):
    scripts=FindSV_modules.main()
    programDirectory = os.path.dirname(os.path.abspath(__file__))
    #run the callers
    if not args.prefix in tracker["FindSV"]["CNVnator"] and not args.prefix in tracker["FindSV"]["FindTranslocations"]:    
        outputVCF,sbatch_ID=calling(args,config,output,scripts,programDirectory,account)
        caller_output=outputVCF.split()

        tracker = tracking_module.add_sample(args.prefix,args.bam,[caller_output[0],caller_output[1]],sbatch_ID[0],"FindTranslocations",output,account,tracker)
        tracker = tracking_module.add_sample(args.prefix,args.bam,[caller_output[2]],sbatch_ID[1],"CNVnator",output,account,tracker)
        caller_vcf=outputVCF
    else:
        tracker=tracking_module.update_status(args.prefix,"CNVnator",output,tracker)
        tracker=tracking_module.update_status(args.prefix,"FindTranslocations",output,tracker)
        vcf_output=tracker["FindSV"]["CNVnator"][args.prefix]["output"]+tracker["FindSV"]["FindTranslocations"][args.prefix]["output"]
        caller_vcf=" ".join(vcf_output)
        sbatch_ID=[tracker["FindSV"]["CNVnator"][args.prefix]["sbatch"],tracker["FindSV"]["FindTranslocations"][args.prefix]["sbatch"]]
    return(tracker,caller_vcf,sbatch_ID)

def run_combine(tracker,args,output,config,account,caller_vcf,sbatch_ID):
    #fetch the scripts
    scripts=FindSV_modules.main()
    programDirectory = os.path.dirname(os.path.abspath(__file__))
    if not args.prefix in tracker["FindSV"]["combine"]: 
        outputVCF,combine_ID=combine_module(args,config,output,scripts,programDirectory,caller_vcf,sbatch_ID,account)
        tracker=tracking_module.add_sample(args.prefix,caller_vcf,[outputVCF],combine_ID,"combine",output,account,tracker)
        combine_vcf=outputVCF
        
    else:
        tracker=tracking_module.update_status(args.prefix,"combine",output,tracker)
        combine_vcf=tracker["FindSV"]["combine"][args.prefix]["output"]
        combine_ID=tracker["FindSV"]["combine"][args.prefix]["sbatch"]
        combine_vcf=combine_vcf[0]

    return(tracker,combine_vcf,combine_ID)        

def run_annotation(tracker,args,output,config,account,combine_vcf,combine_ID):
    #fetch the scripts
    scripts=FindSV_modules.main()
    programDirectory = os.path.dirname(os.path.abspath(__file__))
    if not args.prefix in tracker["FindSV"]["annotation"]:
        outputVCF,annotation_ID=annotation(args,config,output,scripts,programDirectory,combine_vcf,combine_ID,account)   
        tracker=tracking_module.add_sample(args.prefix,combine_vcf,[outputVCF],annotation_ID,"annotation",output,account,tracker)
    else:
        tracker=tracking_module.update_status(args.prefix,"annotation",output,tracker)
    return(tracker)
    
#the module used to perform the variant calling
def calling(args,config,output,scripts,programDirectory,account):
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
            FT=scripts["FindSV"]["header"].format(account=account,time="30:00:00",name=job_name,filename=process_files)
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
            CNVNator=scripts["FindSV"]["header"].format(account=account,time="30:00:00",name=job_name,filename=process_files)
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
def combine_module(args,config,output,scripts,programDirectory,input_vcf,sbatch_ID,account):
    general_config=config["FindSV"]["general"]
    annotation_config=config["FindSV"]["annotation"]
    #combine module; combine all the caller modules into one VCF
    prefix=args.prefix
    output_prefix=os.path.join(output,prefix)
    job_name="combine_{}".format(prefix)
    process_files=os.path.join(output,"slurm/combine/",job_name)
    merge_VCF_path=annotation_config["DB"]["DB_script_path"]
    contig_sort=os.path.join(programDirectory,"internal_scripts","contigSort.py")
    combine=scripts["FindSV"]["header"].format(account=account,time="3:00:00",name=job_name,filename=process_files)
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
def annotation(args,config,output,scripts,programDirectory,outputVCF,combine_ID,account):
    general_config=config["FindSV"]["general"]
    prefix=args.prefix
    output_prefix=os.path.join(output,prefix)

    #annotation module; filter and annotate the samples
    annotation_config=config["FindSV"]["annotation"]
    job_name="annotation_{}".format(prefix)
    process_files=os.path.join(output,"slurm/annotation/",job_name)
    annotation = scripts["FindSV"]["header"].format(account=account,time="10:00:00",name=job_name,filename=process_files)
    annotation += scripts["FindSV"]["annotation"]["header"].format(combine_script_id=combine_ID)
    
    #if uppmax modules are chosen, the vep module is loaded
    if not general_config["UPPMAX"] == "":
        annotation +=scripts["FindSV"]["UPPMAX"].format(modules="bioinfo-tools vep")
    #otherwise if the vep conda environment is installed, we will use it
    elif not config["FindSV"]["conda"]["vep"] == "":
        annotation +=scripts["FindSV"]["conda"].format(environment="VEP_FINDSV")
        
    #add vep annotation
    cache_dir=""
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

 #create a cleaned vcf
    inputVCF=outputVCF
    outputVCF=output_prefix+"_cleaned.vcf"
    clean_VCF_path=os.path.join(programDirectory,"internal_scripts","cleanVCF.py")
    annotation += scripts["FindSV"]["annotation"]["cleaning"].format(output=output_prefix,VCFTOOLS_path=clean_VCF_path,input_vcf=inputVCF,output_vcf=outputVCF)
    
    #merge the breakpoints
    inputVCF=outputVCF
    outputVCF=output_prefix+"_merged.vcf"
    contig_sort=os.path.join(programDirectory,"internal_scripts","contigSort.py")
    annotation +=scripts["FindSV"]["conda"].format(environment="numpy_FINDSV")
    annotation += scripts["FindSV"]["annotation"]["merge"].format(merge_vcf_path=annotation_config["DB"]["DB_script_path"],input_vcf=inputVCF,output_vcf=outputVCF)
    
    #sort according to the contig order of the reference
    inputVCF=outputVCF
    outputVCF=output_prefix+"_contigSort.vcf"
    if not general_config["UPPMAX"] == "":
        annotation +=scripts["FindSV"]["UPPMAX"].format(modules="bioinfo-tools samtools")
    #if we are not on uppmax and the samtools conda module is installed
    elif not config["FindSV"]["conda"]["samtools"] == "":
        annotation +=scripts["FindSV"]["conda"].format(environment="samtools_FINDSV")
    
    annotation += scripts["FindSV"]["annotation"]["sort"].format(input_vcf=inputVCF,output_vcf=outputVCF,contig_sort_path=contig_sort,bam_path=args.bam)


    
    inputVCF=outputVCF
    outputVCF=output_prefix+"_annotator.vcf"
    #use annotator to add omim information
    annotator=os.path.join(programDirectory,"internal_scripts","the_annotator.py")
    annotator_db=os.path.join(programDirectory,"gene_keys")
    annotation += scripts["FindSV"]["annotation"]["annotator"].format(annotator_path=annotator,folder_path=annotator_db,input_vcf=inputVCF,output_vcf=outputVCF)
    
    inputVCF=outputVCF
    #add genmod annotation
    if not annotation_config["GENMOD"]["GENMOD_rank_model_path"] == "":
        #use the genmod conda module if the user wishes to do so
        if not config["FindSV"]["conda"]["samtools"] == "":
            annotation +=scripts["FindSV"]["conda"].format(environment="GENMOD_FINDSV")
        inputVCF=outputVCF
        outputVCF=output_prefix+"_genmod.vcf"
        genmod_sort=os.path.join(programDirectory,"internal_scripts","genmod_stable_sort.py")
        annotation += scripts["FindSV"]["annotation"]["GENMOD"].format(genmod_score_path=annotation_config["GENMOD"]["GENMOD_rank_model_path"],output=output_prefix,input_vcf=inputVCF,output_vcf=outputVCF)

    #add frequency database annotation
    if not annotation_config["DB"]["DB_script_path"] == "" and not annotation_config["DB"]["DB_path"] == "":
        inputVCF=outputVCF
        outputVCF=output_prefix+"_finished.vcf"
        annotation +=scripts["FindSV"]["conda"].format(environment="numpy_FINDSV")
        annotation += scripts["FindSV"]["annotation"]["DB"].format(query_script=annotation_config["DB"]["DB_script_path"],output=output_prefix,db_folder_path=annotation_config["DB"]["DB_path"],input_vcf=inputVCF,output_vcf=outputVCF)
    
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
