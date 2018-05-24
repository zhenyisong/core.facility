# @author Yisong
# @since  2018-05-04
# @update 2018-05-08
# https://groups.google.com/a/continuum.io/forum/#!topic/anaconda/Ef7sV_Y1wY4
# Anaconda: Permanently include external packages (like in PYTHONPATH)
#---
from setuptools import setup

setup(
    name         = 'completeQC',
    version      = '0.12',
    description  = 'the QC pipeline designed for GuoZhong',
    url          = 'https://github.com/zhenyisong/core.facility/',
    author       = 'Yisong Zhen',
    author_email = 'zhenyisong@cardiosignal.org',
    packages     = ['completeQC'],
    license      = 'MIT',
    platforms    = 'Linux',
    scripts      = ['script/quality_control_industry.py'],
    python_requires      = '>=3.6, < 4',
    install_requires     = [ 'plumbum','numpy','psutil','setuptools'],
    include_package_data = True,
    long_description     = open('README.md').read()
)