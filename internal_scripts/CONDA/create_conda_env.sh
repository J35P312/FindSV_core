cd $1
conda config --add channels r
conda config --add channels bioconda

conda env create --force --name GENMOD_FINDSV -f GENMOD.yml
conda env create --force --name samtools_FINDSV -f SAMTOOLS.yml
conda env create --force --name VEP_FINDSV -f VEP_FINDSV.yml
