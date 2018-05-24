# quality control pipeline (python) implemented using industry standard

## Project aim

This project is the initial step to provide a deep service 
for the NGS platform users at GuoZhong to pre-analyze their data from
the outside sequencing services. This tool is suggested to be used
as an independent source to assess the sequencing results, which include, 
but not limited to,
[] check ribosome rRNA contamination
[] check the strandedness of the library
[] check the read length
[] count the total read number
[] check the unique mapping reads percentage
[] check the mycoplasma contamination


## Installing
you should follow installation manual in this folder and 
the step and requirement to install this
python QC pipeline

```
pip install --no-cache-dir --index-url https://test.pypi.org/simple/ completeQC
```
## Usage

use this command to read the help menu.

```
quality_control_industry.py -h
```

or use this command to perform test set checking

```
quality_control_industry.py  -g 'mm10' -l 'PE' -s 'RF'  -u True
```

## Authors
Dr. Yisong Zhen at [CardioSignal](http://www.cardiosignal.org/)
Yisong has his [Linkin account](https://www.linkedin.com/in/yisongzhen/)
## Acknowledgments
Jingzhou Chen, [Deputy Director of GuoZhong](http://www.sklcvd.org/WebShowPage/PsnlIntroDetail.aspx?indexID=6&smallID=30&id=71)