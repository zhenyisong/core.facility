# quality control pipeline (python) implemented using industry standard

## Project aim

This project is the initial step to provide a deep service 
for the NGS platform users at GuoZhong to pre-analyze their data  and double-check them
from the outside sequencing services. This tool is suggested to be used
as an independent source to assess the sequencing results, which include, 
but not limited to,

- [x] check ribosome rRNA contamination
- [x] check the strandedness of the library
- [x] check the read length
- [x] count the total read number
- [x] check the unique mapping reads percentage
- [x] check the mycoplasma contamination


## Installing
you should follow [installation manual in this folder](install.md) and 
the step and requirement to install this
python QC pipeline

you should first create and run a specific conda envieroment.
This is a brief step below.

```console
source activate <your_env_name>
```
and then

```console
pip install --no-cache-dir --index-url https://test.pypi.org/simple/ completeQC
```
## Usage

use the following command to read the help menu.

```console
quality_control_industry.py -h
```

or use the following command to perform debug procedure

```console
quality_control_industry.py  -g 'mm10' -l 'PE' -s 'RF'  -u True
```

## Uninstall

```console
pip uninstall completeQC
```

## Authors
Dr. Yisong Zhen hosts his personal web source for cardiovascular community at [CardioSignal](http://www.cardiosignal.org/).
Yisong has his [Linkin account](https://www.linkedin.com/in/yisongzhen/)

## Acknowledgments
Jingzhou Chen, [Deputy Director of GuoZhong](http://www.sklcvd.org/WebShowPage/PsnlIntroDetail.aspx?indexID=6&smallID=30&id=71)