# [TEclass](https://www.compgen.uni-muenster.de/tools/teclass/index.hbi?)
Repeat sequence classification

## Requirements (这两个软件下载和编译可能有问题，在bin里面有编译好的文件在linux中可以直接使用)
* [lvq_pak](http://www.cis.hut.fi/research/som-research/nnrc-programs.shtml)
  * [download](http://www.cis.hut.fi/research/lvq_pak/lvq_pak-3.1.tar)
* [librf](http://mtv.ece.ucsb.edu/benlee/librf.html)
  * [download](https://github.com/tearshark/librf/archive/refs/tags/2.9.10.tar.gz)

## Database
* [classifiers](https://www.compgen.uni-muenster.de/tools/teclass/download/classifiers.tar.gz)
  
Third-party
-----------
```
cd /Work #软件安装路径
mkdir TEclass
cd TEclass
git clone https://github.com/zxgsy520/TEclass.git
mv TEclass temp
mkdir v2.1.4
cd v2.1.4
conda env create --prefix=/Work/TEclass/v2.1.4 -f environment.yml
#or(或者如下)
#conda create --prefix=/Work/TEclass/v2.1.4 -c bioconda blast-legacy=2.2.26 glimmer=3.02 -c conda-forge libsvm=325
cd bin
chmod 755 ../../temp/*
cp ../../temp/* .
ln -s TEclassBuild.pl TEclassBuild
ln -s TEclassTest.pl TEclassTest
ln -s /Work/TEclass_db/classifiers/   #classifiers 的模型数据库
```
