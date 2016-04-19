import yaml
import os
import subprocess


#all possible command line questions to set up the config file
def questions(question):
    answer={}
    answer["account"]="enter a project account(or leave blank for no account)"
    answer["output"]="enter set standard output path(or leave blank)"
    answer["TMPDIR"]="set the path to the tmp/scratch directory, or leave blank to perform computations directly on the fileystem of the bam file"
    answer["UPPMAX"]="type anything to use the uppmax module system, leave blank otherwise"
    answer["tracking"]="type anything to use the FindSV tracking system to remember which files has been analysed, leave blank to skip tracking"
    answer["FT_path"]="enter the path to the FindTranslocations executable"
    answer["minimum_supporting_pairs"]="enter the minimum number of pairs used to call a variant(Findtranslocations)"
    answer["CNVnator_path"]="set the cnvnator path"
    answer["CNVnator2vcf_path"]="set the path to the cnvnator2VCF.pl script"
    answer["ROOTSYS"]="set the path to the root folder(only needed if root is not added to path already)"
    answer["bin_size"]=""
    answer["reference_dir"]="set the path containing the reference files"
    answer["VEP.pl_path"]=""
    answer["cache_dir"]="change the vep directory, or use the default(visit the vep documentation for more info)"
    answer["port"]=""
    answer["DB_script_path"]="set the path to the internal db query script(query_db.py)"
    answer["DB_path"]="add the path to a frequency database folder, or leave blank"
    answer["overlap_parameter"]=""
    answer["GENMOD_rank_model_path"]="add the path to a genmod ini file, or leave blank"

    print(answer[question])
    return(raw_input())

#install the conde environments
def conda_environments(config,programDirectory):
    print("installing conda environments")    
    config["FindSV"]["conda"]["genmod"]="True"
    config["FindSV"]["conda"]["samtools"]="True"
    config["FindSV"]["conda"]["vep"]="True"
    config["FindSV"]["annotation"]["DB"]["DB_script_path"]=os.path.join(programDirectory,"SVDB/SVDB.py")
    command=["{} {} {}".format(os.path.join(programDirectory,"internal_scripts/CONDA/create_conda_env.sh"), os.path.join(programDirectory,"internal_scripts/CONDA/"), programDirectory)]
    tmp=subprocess.check_output(command,shell = True)
    
    return(config)
    
#install FT
def FT_install(config,programDirectory,UPPMAX):
    print("installing FindTranslocations")
    config["FindSV"]["calling"]["FT"]["FT_path"]=os.path.join(programDirectory,"FindTranslocations/bin/FindTranslocations")
    
    command=["{} {}".format(os.path.join(programDirectory,"internal_scripts/install_FT.sh"),programDirectory)]
    tmp=subprocess.check_output(command,shell = True)
    
    return(config)

def cnvnator_install(config,programDirectory,args):
    config["FindSV"]["calling"]["CNVnator"]["CNVnator_path"]=os.path.join(programDirectory,"CNVnator_v0.3.1/src/cnvnator")
    config["FindSV"]["calling"]["CNVnator"]["CNVnator2vcf_path"]=os.path.join(programDirectory,"CNVnator_v0.3.1/cnvnator2VCF.pl")

    if not args.no_root:
        config["FindSV"]["calling"]["CNVnator"]["ROOTSYS"]=os.path.join(programDirectory,"root")

        if not args.compile_root:
            env = os.environ.copy()
            env["ROOTSYS"]= config["FindSV"]["calling"]["CNVnator"]["ROOTSYS"]
            env["ROOTLIB"]= os.path.join(config["FindSV"]["calling"]["CNVnator"]["ROOTSYS"],"lib")
            command=["{} {}".format(os.path.join(programDirectory,"internal_scripts/install_cnvnator.sh"),programDirectory)]
            tmp=subprocess.check_output(command,shell = True,env=env)
        else:
            command=["{} {}".format(os.path.join(programDirectory,"internal_scripts/install_cnvnator_compile_root.sh"),programDirectory)]
            tmp=subprocess.check_output(command,shell = True)
    else:
        command=["{} {}".format(os.path.join(programDirectory,"internal_scripts/install_cnvnator_no_root.sh"),programDirectory)]
        tmp=subprocess.check_output(command,shell = True)

        
    return(config)
    
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
    

#generate an empty config file
def generate_config(programDirectory):
    config=[]
    general={"account":"","output":"","TMPDIR":"","UPPMAX":"","tracking":"enabled"}
    calling={"FT":{"FT_path":"","minimum_supporting_pairs":"6"},"CNVnator":     {"CNVnator_path":"cnvnator","CNVnator2vcf_path":"cnvnatorVCF.pl","ROOTSYS":"","bin_size":"1000","reference_dir":""}}
    conda={"genmod":"","samtools":"","vep":""}
    filter={"VEP":{"VEP.pl_path":"","cache_dir":"","port":"3337"},"DB":{"DB_script_path":"","DB_path":"","overlap_parameter":"0.7"},"GENMOD":{"GENMOD_rank_model_path":""}}

    config=[{"FindSV":{"general":general,"calling":calling,"annotation":filter,"conda":conda}}]

    f = open(os.path.join(programDirectory,"config.txt"), 'w')
    for entry in config:
        f.write(yaml.dump(entry).strip())
        
#UPPMAX configuration of the config file
def UPPMAX(programDirectory):
    config=readconfig(0,os.path.join(programDirectory,"config.txt"));
    config["FindSV"]["general"]["UPPMAX"]="True"
    
    #install FT and the conda environments
    config=FT_install(config,programDirectory,True)
    config=conda_environments(config,programDirectory)
    
    print("leave blank to exit setup, enter anything else to continue setting up path to the annotation tools and references")
    #the user may choose the set up everything manually(probably faster for the experienced user)
    if raw_input() == "":
        print("installation finished, please set up the path to cnvnators reference directory")
        print("also add the frequency DB path, and genmod ini file path, or these annotation tools will be skipped")
    #the user may also choose to go through the last options step by step
    else:
        config["FindSV"]["calling"]["CNVnator"]["reference_dir"]=questions("reference_dir")
        config["FindSV"]["annotation"]["GENMOD"]["GENMOD_rank_model_path"]=questions("GENMOD_rank_model_path")
        config["FindSV"]["annotation"]["DB"]["DB_path"]=questions("DB_path")
        config["FindSV"]["annotation"]["VEP"]["cache_dir"]=questions("cache_dir")
        config["FindSV"]["general"]["account"]=questions("account")
        config["FindSV"]["general"]["output"]=questions("output")
        print("installation complete!")
    f = open(os.path.join(programDirectory,"config.txt"), 'w')
    f.write(yaml.dump(config).strip())

#setup the config file so that conda environments are used, do not install any other tool and let the user set these tools manually
def conda(programDirectory):
    config=readconfig(0,os.path.join(programDirectory,"config.txt"));
    #install the conda environments
    config=conda_environments(config,programDirectory)
    print("leave blank to exit setup, enter anything else to continue setting up path to the annotation tools and references")
    #the user may choose the set up everything manually(probably faster for the experienced user)
    if raw_input() == "":
        print("installation finished, to properly run FindSV, the path to cnvnator,cnvnator2vcf.pl, and FindTranslocations must be added to the config file")
        print("also add the frequency DB path, and genmod ini file path, or these annotation tools will be skipped")
    #the user may also choose to go through the last options step by step
    else:
        config["FindSV"]["calling"]["FT"]["FT_path"]=questions("FT_path")
        config["FindSV"]["annotation"]["DB"]["DB_script_path"]=questions("DB_script_path")
        config["FindSV"]["calling"]["CNVnator"]["CNVnator_path"]=questions("CNVnator_path")
        config["FindSV"]["calling"]["CNVnator"]["CNVnator2vcf_path"]=questions("CNVnator2vcf_path")
        config["FindSV"]["calling"]["CNVnator"]["reference_dir"]=questions("reference_dir")
        config["FindSV"]["calling"]["CNVnator"]["ROOTSYS"]=questions("ROOTSYS")
        config["FindSV"]["annotation"]["GENMOD"]["GENMOD_rank_model_path"]=questions("GENMOD_rank_model_path")
        config["FindSV"]["annotation"]["DB"]["DB_path"]=questions("DB_path")
        config["FindSV"]["annotation"]["VEP"]["cache_dir"]=questions("cache_dir")
        config["FindSV"]["general"]["account"]=questions("account")
        config["FindSV"]["general"]["output"]=questions("output")
        print("installation complete!")
    
    f = open(os.path.join(programDirectory,"config.txt"), 'w')
    f.write(yaml.dump(config).strip())
    
#automatic install on any system
def auto(programDirectory,args):
    config=readconfig(0,os.path.join(programDirectory,"config.txt"));
    #install conda environments
    config=conda_environments(config,programDirectory)
    #install ft
    config=FT_install(config,programDirectory,True)
    #install cnvnator
    config=cnvnator_install(config,programDirectory,args)
    print("leave blank to exit setup, enter anything else to continue setting up path to the annotation tools and references")
    if raw_input() == "":
        print("installation finished!")
        print("Add the frequency DB path, and genmod ini file path, or these annotation tools will be skipped")
    else:
        config["FindSV"]["calling"]["CNVnator"]["reference_dir"]=questions("reference_dir")
        config["FindSV"]["annotation"]["GENMOD"]["GENMOD_rank_model_path"]=questions("GENMOD_rank_model_path")
        config["FindSV"]["annotation"]["DB"]["DB_path"]=questions("DB_path")
        config["FindSV"]["annotation"]["VEP"]["cache_dir"]=questions("cache_dir")
        config["FindSV"]["general"]["account"]=questions("account")
        config["FindSV"]["general"]["output"]=questions("output")
        
        print("installation complete!")
    f = open(os.path.join(programDirectory,"config.txt"), 'w')
    f.write(yaml.dump(config).strip())
