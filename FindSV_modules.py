import yaml


def main():
    #source activate for conda modules
    conda="\nsource activate {environment}\n"
    
    #if uppmax module system is activted, the available uppmax modules will be used
    uppmax="\nmodule load {modules}\n"
    
    #Slurm header
    header="""#! /bin/bash -l
#SBATCH -A {account}
#SBATCH -o {filename}.out
#SBATCH -e {filename}.err
#SBATCH -J {name}.job
#SBATCH -p core
#SBATCH -t {time}"""

    #the FT script
    FT="""
{FT_path} --sv  -bam {bam_path} -pairs {minimum_suporting} -o {output}
rm {output}.tab"""
    
    #the cnvnator script
    ROOTPATH ="""
START=$(pwd)
cd {rootdir}/bin
source thisroot.sh
cd $START

"""

    CNVnator="""
{CNVnator_path} -root {output}.root -tree {bam_path}
{CNVnator_path} -root {output}.root -his {bin_size} -d {reference_dir}
{CNVnator_path} -root {output}.root -stat {bin_size} >> {output}.cnvnator.log
{CNVnator_path} -root {output}.root -partition {bin_size}
{CNVnator_path} -root {output}.root -call {bin_size} > {output}.cnvnator.out
{CNVnator2vcf_path} {output}.cnvnator.out  >  {output}_CNVnator.vcf
rm {output}.root
rm {output}.cnvnator.out

"""
    calling={"FT":FT,"CNVnator":CNVnator}
    
    #The combine script
    combine="""
python {merge_vcf_path} --vcf {input_vcf} > {output_vcf}.unsorted
python {contig_sort_path} --vcf {output_vcf}.unsorted --bam {bam_path} > {output_vcf}
rm {output_vcf}.unsorted"""
    combine={"combine":combine}
    
    #The annotation header
    annotation_header="""
#SBATCH -d afterok:{combine_script_id}

"""
    #the DB section
    DB="python {query_script} --variations {input_vcf} --db {db_folder_path} > {output_vcf}\n"
    #the vep section
    VEP="perl {vep_path} --cache --force_overwrite --poly b -i {input_vcf} -o {output_vcf} --buffer_size 5 --port {port} --vcf --per_gene --format vcf  {cache_dir} -q\n"
    UPPMAX_VEP="variant_effect_predictor.pl --cache --force_overwrite --poly b -i {input_vcf}  -o {output_vcf} --buffer_size 5 --port {port} --vcf --per_gene --format vcf  {cache_dir} -q\n"
    #the genmod section
    GENMOD="genmod score -c {genmod_score_path} {input_vcf}  > {output_vcf}\n"
    
    cleaning="python {VCFTOOLS_path} --vcf {input_vcf} > {output_vcf} \n"
    filter={"header":annotation_header,"VEP":VEP,"UPPMAX_VEP":UPPMAX_VEP,"DB":DB,"GENMOD":GENMOD,"cleaning":cleaning}
    
    
    scripts={"FindSV":{"calling":calling,"annotation":filter,"combine":combine,"header":header,"UPPMAX":uppmax,"afterok":"\n#SBATCH -d afterok:{slurm_IDs}\n","ROOTSYS":ROOTPATH,"conda":conda}}
    return(scripts)
