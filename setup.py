#!/usr/bin/env python
import os
from setuptools import setup, find_packages

NAME = "BIOx"
VERSION = "0.1"

DESCRIPTION = "A Python library which PIPA uses to perform bioinformatics analysis."
LONG_DESCRIPTION = open(os.path.join(os.path.dirname(__file__), 'README.rst')).read()
AUTHOR = "Bioinformatics Laboratory at University of Ljubljana"
AUTHOR_EMAIL = "marko.toplak@fri.uni-lj"
URL = "https://bitbucket.org/markotoplak/biox"
LICENSE = "Proprietary software"

if __name__ == '__main__':
    setup(
        name=NAME,
        version=VERSION,
        description=DESCRIPTION,
        long_description=LONG_DESCRIPTION,
        author=AUTHOR,
        author_email=AUTHOR_EMAIL,
        url=URL,
        license=LICENSE,
        packages=find_packages(),
        package_data={},
        classifiers=[
            'Development Status :: 5 - Production/Stable',
            'Environment :: Console',
            'Intended Audience :: Developers',
            'License :: Other/Proprietary License',
            'Operating System :: OS Independent',
            'Programming Language :: Python',
        ],
        include_package_data=True,
        zip_safe=False,
        install_requires=[],
    )
