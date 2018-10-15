### My first Python Packge from mRNA complete quanlity control
[Python package upload howto](https://packaging.python.org/guides/using-testpypi/)
1. modify the local source code
2. update the version number in setup.py

   ```console
   cd sourcecode/core.facility/completeQC
   vi setup.py
   ```
3. generate the local source code package

   ```console
   python setup.py sdist
   ```
4. remove the old package installed locally

    ```console
    pip uninstall completeQC
    ```
5. install the updated pacckage. and test if it is OK.

    ```console
    pip install dist/completeQC-0.13.tar.gz
    ```
6. publish the updated package.

    ```console
    twine upload --repository testpypi dist/*
    ```