import yaml
import os
import subprocess
import shlex
import submit_module

#generate the tracker file
def generate_tracker(directory):

    combine={"combine":{}}
    annotation={"annotation":{}}
    tracker={ "FindSV":{"CNVnator":{},"FindTranslocations":{},"combine":{},"annotation":{} } }
    
    
    f = open(os.path.join(directory,"tracker.yml"), 'w')
    track=[tracker]
    for entry in track:
        f.write(yaml.dump(entry).strip())

#add a sample to the tracker file
def add_sample(ID,input,output,sbatch_id,process,directory,account,tracker):
    tracker["FindSV"][process][ID]={"input":input,"output":output,"sbatch":sbatch_id,"status":"SUBMITED","account":account}
    return(tracker)

#update the status of a sample
def update_status(ID,process,directory,tracker):
    sbatch=tracker["FindSV"][process][ID]["sbatch"]
    
    SLURM_EXIT_CODES = {"PENDING": "PENDING","RUNNING": "RUNNING","RESIZING": "RUNNING","SUSPENDED": "RUNNING","COMPLETED": "COMPLETED","CANCELLED": "CANCELLED","FAILED": "FAILED",
    "TIMEOUT": "FAILED","PREEMPTED": "FAILED","BOOT_FAIL": "FAILED","NODE_FAIL": "FAILED"}

    #try connect to sbatch using sacct, if we cannot fetch the status, there is probably some temporal error
    try:
        check_cl = "sacct -n -j {0} -o STATE".format(sbatch)
        job_status = subprocess.check_output(shlex.split(check_cl))
    except:
        tracker["FindSV"][process][ID]["status"] = "CONNECTION_ERROR"
        f = open(os.path.join(directory,"tracker.yml"), 'w')
        return(tracker)

    if not job_status:
        return(tracker)
    else:
        #if a sample failed or timed out, print it to the stoud
        status=job_status.split()[0].strip("+")
        if SLURM_EXIT_CODES[status] == "FAILED":
            print "sample {0} FAILED:{}".format(status)
        elif SLURM_EXIT_CODES[status] == "TIMEOUT":
            print "sample {0} ERROR:{}".format(status)
        try:
            tracker["FindSV"][process][ID]["status"] = SLURM_EXIT_CODES[status]
        except (IndexError, KeyError, TypeError) as e:
            tracker["FindSV"][process][ID]["status"] = "UNKNOWN"
        return(tracker)

#update the status of each sample within the tracker
def update_tracker(directory):
    with open(os.path.join(directory,"tracker.yml"), 'r') as stream:
        tracker=yaml.load(stream)
        for process in tracker["FindSV"]:
            for sample in tracker["FindSV"][process]:
                tracker=update_status(sample,process,directory,tracker)
    f = open(os.path.join(directory,"tracker.yml"), 'w')
    f.write(yaml.dump(tracker).strip())
    return(tracker)

#restart all steps of all samples within the project               
def full_restart(prefix,output,tracker,args,config):

    print("restarting sample:{}".format(prefix))
    account=tracker["FindSV"]["FindTranslocations"][prefix]["account"]
    args.bam=tracker["FindSV"]["FindTranslocations"][prefix]["input"]
    
    #delete the process
    for process in tracker["FindSV"]:
        del tracker["FindSV"][process][prefix]
    #run the sample
    tracker,caller_vcf,sbatch_ID = submit_module.run_callers(tracker,args,output,config,account)
    #combine them
    tracker,combine_vcf,combine_ID = submit_module.run_combine(tracker,args,output,config,account,caller_vcf,sbatch_ID)
    print(combine_vcf)
    #annotate the vcf  
    tracker=submit_module.run_annotation(tracker,args,output,config,account,combine_vcf,combine_ID)
    return(tracker)

def combine_restart(prefix,output,tracker,args,config):
    print("restarting sample:{}".format(prefix))
    account=tracker["FindSV"]["FindTranslocations"][prefix]["account"]
    args.bam=tracker["FindSV"]["FindTranslocations"][prefix]["input"]

    #clear the combine and annotation status
    for process in tracker["FindSV"]:
        if process == "combine" or process == "annotation":
            del tracker["FindSV"][process][prefix]

    #run the sample
    tracker,caller_vcf,sbatch_ID = submit_module.run_callers(tracker,args,output,config,account)
    #combine them
    tracker,combine_vcf,combine_ID = submit_module.run_combine(tracker,args,output,config,account,caller_vcf,sbatch_ID)
    
    #annotate the vcf  
    tracker=submit_module.run_annotation(tracker,args,output,config,account,combine_vcf,combine_ID)
    return(tracker)


def annotation_restart(prefix,output,tracker,args,config):
    print("restarting sample:{}".format(prefix))
    account=tracker["FindSV"]["FindTranslocations"][prefix]["account"]
    args.bam=tracker["FindSV"]["FindTranslocations"][prefix]["input"]

    #clear the combine and annotation status
    for process in tracker["FindSV"]:
        if process == "annotation":
            del tracker["FindSV"][process][prefix]

    #run the sample
    tracker,caller_vcf,sbatch_ID = submit_module.run_callers(tracker,args,output,config,account)
    #combine them
    tracker,combine_vcf,combine_ID = submit_module.run_combine(tracker,args,output,config,account,caller_vcf,sbatch_ID)
    
    #annotate the vcf  
    tracker=submit_module.run_annotation(tracker,args,output,config,account,combine_vcf,combine_ID)
    return(tracker)

def status(prefix,output,tracker,args,config):
    keyword="FAILED"
    if(args.cancelled):
        keyword="CANCELLED"

    if tracker["FindSV"]["FindTranslocations"][prefix]["status"] == keyword or tracker["FindSV"]["CNVnator"][prefix]["status"] == keyword:
        tracker=full_restart(prefix,output,tracker,args,config)
    elif tracker["FindSV"]["combine"][prefix]["status"] == keyword:
        tracker=combine_restart(prefix,output,tracker,args,config)
    elif tracker["FindSV"]["annotation"][prefix]["status"] == keyword:
        tracker=annotation_restart(prefix,output,tracker,args,config)
    return(tracker)


#this function is used to restart samples based on their status or a selected step of the pipeline
def restart(directory,args,config):
    #update the tracker before loading it
    update_tracker(directory)
    with open(os.path.join(directory,"tracker.yml"), 'r') as stream:
        tracker=yaml.load(stream)
    #restart all samples within the project
    if args.full:
        prefix_list = list(tracker["FindSV"]["FindTranslocations"].keys() )
        for prefix in prefix_list:
            args.prefix=prefix
            tracker=full_restart(prefix,directory,tracker,args,config)

    #restart the combine step of all the samples     
    elif args.combine:
        prefix_list = list(tracker["FindSV"]["combine"].keys() )
        for prefix in prefix_list:
            args.prefix=prefix
            tracker=combine_restart(prefix,directory,tracker,args,config)

    #redo the annotation on all samples
    elif args.annotation:
        prefix_list = list(tracker["FindSV"]["annotation"].keys() )
        for prefix in prefix_list:
            args.prefix=prefix
            tracker=annotation_restart(prefix,directory,tracker,args,config)

    #restart only the failed or cancelled samples
    elif args.cancelled or args.failed:
        prefix_list = list(tracker["FindSV"]["FindTranslocations"].keys() )
        for prefix in prefix_list:
            args.prefix=prefix
            tracker=status(prefix,directory,tracker,args,config)

    f = open(os.path.join(directory,"tracker.yml"), 'w')
    f.write(yaml.dump(tracker).strip())
    return(tracker)
            
