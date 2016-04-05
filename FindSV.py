import argparse
import sys
import yaml
import os
import re
import fnmatch
import tracking_module
import setup
import test_modules
import submit_module



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
    
    account=config["FindSV"]["general"]["account"]
    if args.account:
        account= args.account
        
    #run the callers
    tracker,caller_vcf,sbatch_ID = submit_module.run_callers(tracker,args,output,config,account)
    #combine them
    tracker,combine_vcf,combine_ID = submit_module.run_combine(tracker,args,output,config,account,caller_vcf,sbatch_ID)
    
    #annotate the vcf  
    tracker=submit_module.run_annotation(tracker,args,output,config,account,combine_vcf,combine_ID)


parser = argparse.ArgumentParser("FindSV core module",add_help=False)
parser.add_argument('--bam', type=str,help="analyse the bam file using FindSV")
parser.add_argument("--account", type=str,help="The slurm account that will be used to analyse a sample or folder, overrides the account in the config file")
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
        args.restart= [config["FindSV"]["general"]["output"]]

    config=readconfig(args.config,args);
    
    for tracker in args.restart:
        tracking_module.restart(tracker, args,config)
    
else:
    parser.print_help()
