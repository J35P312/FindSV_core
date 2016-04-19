cd $1
git clone https://github.com/J35P312/SVDB.git
cd SVDB
source activate CYTHON_FINDSV
python setup.py build_ext --inplace
