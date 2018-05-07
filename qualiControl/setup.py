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
    url          = 'https://github.com/zhenyisong/core.facility/',
    author       = 'Yisong Zhen',
    author_email = 'zhenyisong@cardiosignal.org',
    packages     = ['qualiControl'],
    license      = 'MIT',
    long_description = open('README.md').read()
)