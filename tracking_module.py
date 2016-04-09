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
def add_sample(ID,input,output,sbatch_id,process,directory,account):
    with open(os.path.join(directory,"tracker.yml"), 'r') as stream:
        tracker=yaml.load(stream)
    tracker["FindSV"][process][ID]={"input":input,"output":output,"sbatch":sbatch_id,"status":"SUBMITED","account":account}
    
    f = open(os.path.join(directory,"tracker.yml"), 'w')
    track=[tracker]
    for entry in track:
        f.write(yaml.dump(entry).strip())

#update the status of a sample
def update_status(ID,process,directory):
    with open(os.path.join(directory,"tracker.yml"), 'r') as stream:
        tracker=yaml.load(stream)
    sbatch=tracker["FindSV"][process][ID]["sbatch"]
    
    SLURM_EXIT_CODES = {"PENDING": "PENDING","RUNNING": "RUNNING","RESIZING": "RUNNING","SUSPENDED": "RUNNING","COMPLETED": "COMPLETED","CANCELLED": "CANCELLED","FAILED": "FAILED",
    "TIMEOUT": "TIMEOUT","PREEMPTED": "FAILED","BOOT_FAIL": "FAILED","NODE_FAIL": "FAILED"}

    #try connect to sbatch using sacct, if we cannot fetch the status, there is probably some temporal error
    try:
        check_cl = "sacct -n -j {0} -o STATE".format(sbatch)
        job_status = subprocess.check_output(shlex.split(check_cl))
    except:
        tracker["FindSV"][process][ID]["status"] = "CONNECTION_ERROR"
        f = open(os.path.join(directory,"tracker.yml"), 'w')
        track=[tracker]
        for entry in track:
            f.write(yaml.dump(entry).strip())
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
        f = open(os.path.join(directory,"tracker.yml"), 'w')
        track=[tracker]
        for entry in track:
            f.write(yaml.dump(entry).strip())
        return(tracker)

#update the status of each sample within the tracker
def update_tracker(directory):
    with open(os.path.join(directory,"tracker.yml"), 'r') as stream:
        tracker=yaml.load(stream)
        for process in tracker["FindSV"]:
            for sample in tracker["FindSV"][process]:
                update_status(sample,process,directory)

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
    
    #annotate the vcf  
    tracker=submit_module.run_annotation(tracker,args,output,config,account,combine_vcf,combine_ID)
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
            
    elif args.combine:
        print("not implemented")
        pass
    #redo the annotation on all samples
    elif args.annotation:
        print("not implemented")
        pass
    #restart only the failed or cancelled samples
    elif args.cancelled or failed:
        print("not implemented")
        if args.cancelled:
            pass
        if args.failed:
            pass
    return(tracker)
            
