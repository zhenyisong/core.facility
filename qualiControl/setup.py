# @author Yisong
# @since 2018-05-04
# https://groups.google.com/a/continuum.io/forum/#!topic/anaconda/Ef7sV_Y1wY4
# Anaconda: Permanently include external packages (like in PYTHONPATH)
#---
from setuptools import setup

setup(
    name         = 'completeQualityControl',
    version      = '0.1dev',
    description  = 'the QC pipeline designed for GuoZhong',
    url          = '',
    author       = 'Yisong Zhen',
    author_email = 'zhenyisong@cardiosignal.org',
    packages     = ['qualiControl'],
    license      = 'Creative Commons Attribution-Noncommercial-Share Alike license',
    long_description = open('README.txt').read()
)