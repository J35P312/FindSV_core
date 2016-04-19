cd $1
conda config --add channels r
conda config --add channels bioconda

conda env create --force --name CYTHON_FINDSV -f CYTHON.yml
conda env create --force --name GENMOD_FINDSV -f GENMOD.yml
conda env create --force --name samtools_FINDSV -f SAMTOOLS.yml
conda env create --force --name VEP_FINDSV -f VEP_FINDSV.yml

echo $2
cd $2
git clone https://github.com/J35P312/SVDB.git
cd SVDB
python setup.py build_ext --inplace
