# coding: utf-8

import sys
from setuptools import setup, find_packages

NAME = "act"
VERSION = "0.0.1"

# To install the library, run the following
#
# python setup.py install
#
# prerequisite: setuptools
# http://pypi.python.org/pypi/setuptools

setup(
    name=NAME,
    version=VERSION,
    description="Allele Calling Service",
    author_email="mhalagan@nmdp.org",
    url="",
    keywords=["Swagger", "Allele Calling Service"],
    packages=find_packages(exclude=['contrib', 'docs', 'tests']),
    package_data={'': ['swagger/swagger.yaml'],'gfe_typing': ['data/*.structure']},
    install_requires=['py2neo', 'pandas', 'connexion', 'nose', 'biopython'],
    include_package_data=True,
    long_description="""\
    The Allele Calling  service provides an API for converting raw sequence data to GFE. It provides both a RESTful API and a simple user interface for converting raw sequence data to GFE results. Sequences can be submitted one at a time or as a fasta file. This service uses &lt;a href&#x3D;\&quot;https://github.com/nmdp-bioinformatics/service-feature\&quot;&gt;nmdp-bioinformatics/service-feature&lt;/a&gt; for encoding the raw sequence data and &lt;a href&#x3D;\&quot;https://github.com/nmdp-bioinformatics/HSA\&quot;&gt;nmdp-bioinformatics/HSA&lt;/a&gt; for aligning the raw sequence data. The code is open source, and available on &lt;a href&#x3D;\&quot;https://github.com/nmdp-bioinformatics/service-act\&quot;&gt;GitHub&lt;/a&gt;.&lt;br&gt;&lt;br&gt;Go to &lt;a href&#x3D;\&quot;http://service-act.readthedocs.io\&quot;&gt;service-act.readthedocs.io&lt;/a&gt; for more information
    """
)

