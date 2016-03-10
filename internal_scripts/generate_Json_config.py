import yaml

config=[]
general={"account":"","output":"","TMPDIR":"","UPPMAX":""}
calling={"FT":{"FT_path":"","minimum_supporting_pairs":"6"},"CNVnator":{"CNVnator_path":"cnvnator","CNVnator2vcf_path":"CNVnator2vcf.pl","ROOTSYS":"","bin_size":"1000","reference_dir":""}}
conda={"genmod":"","samtools":"","vep":""}
filter={"VEP":{"VEP.pl_path":"","cache_dir":"","port":"3337"},"DB":{"DB_script_path":"","DB_path":"","overlap_parameter":"0.7"},"GENMOD":{"GENMOD_rank_model_path":""}}

config=[{"FindSV":{"general":general,"calling":calling,"annotation":filter,"conda":conda}}]
for entry in config:
    print(yaml.dump(entry).strip())
