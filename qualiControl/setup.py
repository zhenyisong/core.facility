# @author Yisong
# @since  2018-05-04
# @update 2018-05-08
# https://groups.google.com/a/continuum.io/forum/#!topic/anaconda/Ef7sV_Y1wY4
# Anaconda: Permanently include external packages (like in PYTHONPATH)
#---
from setuptools import setup

setup(
    name         = 'completeQualityControl',
    version      = '0.11dev',
    description  = 'the QC pipeline designed for GuoZhong',
    url          = 'https://github.com/zhenyisong/core.facility/',
    author       = 'Yisong Zhen',
    author_email = 'zhenyisong@cardiosignal.org',
    packages     = ['code','test'],
    license      = 'MIT',
    scripts      = ['code/quality_control_industry.py'],
    install_requires = ['plumbum','numpy','psutil'],
    long_description = open('README.md').read()
)