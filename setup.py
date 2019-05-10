#!/usr/bin/env python

from setuptools import setup, Extension

version = open('VERSION').read().strip()
download_url = "https://github.com/vishnubob/ssw/archive/v%s.tar.gz" % version
long_description = \
"""ssw is a fast implementation of the Smith-Waterman algorithm, which
uses the Single-Instruction Multiple-Data (SIMD) instructions to parallelize
the algorithm at the CPU level.  This repository wraps the SSW library into an
easy to install, high-level python interface with no external library dependancies.

The SSW library is written by Mengyao Zhao and Wan-Ping Lee, and this python
interface is maintained by Giles Hall.
"""

libssw_ext = {
    "sources": ["src/ssw/ssw.c"],
    "include_dirs": ["src/ssw"],
}

config = {
    "name": "ssw", 
    "version": version,
    "description": "Smith-Waterman Sequence Aligner",
    "author": "Giles Hall",
    "author_email": "giles@polymerase.org",
    "package_dir": {"ssw": "src"},
    "packages": ["ssw"],
    "classifiers": [
        "Development Status :: 5 - Production/Stable",
        "Programming Language :: Python",
        "Programming Language :: Python :: 2",
        "Programming Language :: Python :: 2.7",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.3",
        "Programming Language :: Python :: 3.4",
        "Programming Language :: Python :: 3.5",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "License :: OSI Approved :: MIT License",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    "keywords": [
        "bioinformatics",
        "sequence alignment",
        "smith-waterman",
        "genomics",
        "proteomics"
    ],
    "install_requires": [
        "six",
    ],
    "platforms": "any",
    "ext_modules": [Extension("_libssw", **libssw_ext)],
    "zip_safe": False,
    "download_url": download_url,
    "url": "https://github.com/vishnubob/ssw",
    "long_description": long_description,
    "long_description_content_type": "text/x-rst",
}

if __name__ == "__main__":
    setup(**config)
