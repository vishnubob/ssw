#!/usr/bin/env python

from setuptools import setup, Extension

libssw_ext = {"sources": ["src/ssw/ssw.c"], "include_dirs": ["src/ssw"]}

config = {
    "name": "ssw", 
    "version": "0.1",
    "description": "Smith-Waterman Sequence Aligner",
    "author": "Giles Hall",
    "author_email": "giles@polymerase.org",
    "package_dir": {"ssw": "src"},
    "packages": ["ssw"],
    "classifiers": [
        "Development Status :: 4 - Beta",
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
    "ext_modules": [Extension("_libssw", **libssw_ext)],
    "zip_safe": False,
}

if __name__ == "__main__":
    setup(**config)
