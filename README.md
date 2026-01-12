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
#conda create --prefix=/Work/TEclass/v2.1.4 -c bioconda -c conda-forge blast-legacy=2.2.26 glimmer=3.02 libsvm=325
#or(或者如下)
#mamba create --prefix=/Work/TEclass/v2.1.4 -c bioconda -c conda-forge blast-legacy=2.2.26 glimmer=3.02 libsvm=325
cd bin
chmod 755 ../../temp/*.pl
chmod 755 ../../temp/bin/*
mv ../../temp/*.pl .
mv ../../temp/*.pm .
mv ../../temp/bin/* .
ln -s TEclassBuild.pl TEclassBuild
ln -s TEclassTest.pl TEclassTest
ln -s /Work/TEclass_db/classifiers/   #classifiers 的模型数据库
```

## References
-----------
[1] Zhang L, Shi Y, Gong W, et al. The tetraploid Camellia oleifera genome provides insights into evolution, agronomic traits, and genetic architecture of oil Camellia plants. Cell Rep. 2024;43(11):114902. doi:[10.1016/j.celrep.2024.114902](https://www.cell.com/cell-reports/fulltext/S2211-1247(24)01253-1?uuid=uuid%3Ae8f16f2d-e8b5-4510-a644-c4a910ae9884)
