cd $1
#install root and ad it to path
wget https://root.cern.ch/download/root_v5.34.34.Linux-ubuntu14-x86_64-gcc4.8.tar.gz
tar -xvf root_v5.34.34.Linux-ubuntu14-x86_64-gcc4.8.tar.gz
#cd root/bin/
#source thisroot.sh
#cd $1
#download and install cnvnator
wget https://github.com/abyzovlab/CNVnator/releases/download/v0.3.1/CNVnator_v0.3.1.zip
unzip CNVnator_v0.3.1.zip
cd CNVnator_v0.3.1
cd src/samtools
make 
cd ..
make
