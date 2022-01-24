# [TEclass](https://www.compgen.uni-muenster.de/tools/teclass/index.hbi?)
Repeat sequence classification

## Requirements
* [lvq_pak](http://www.cis.hut.fi/research/som-research/nnrc-programs.shtml)
  * [download](http://www.cis.hut.fi/research/lvq_pak/lvq_pak-3.1.tar)
* [librf](http://mtv.ece.ucsb.edu/benlee/librf.html)
  * [download](https://github.com/tearshark/librf/archive/refs/tags/2.9.10.tar.gz)
  
Third-party
-----------
```
cd /Work
mkdir TEclass
cd TEclass
git clone https://github.com/zxgsy520/TEclass.git
mv TEclass v2.1.4
cd v2.1.4
conda env create --prefix=/Work/TEclass/v2.1.4 -f environment.yml

```
