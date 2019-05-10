#!/usr/bin/env python

from setuptools import setup, Extension

version = open('VERSION').read().strip()
download_url = "https://github.com/vishnubob/ssw/archive/v%s.tar.gz" % version

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
        "Development Status :: 4 - Beta",
        "Programming Language :: Python",
        "Programming Language :: Python :: 2",
        "Programming Language :: Python :: 2.6",
        "Programming Language :: Python :: 2.7",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.3",
        "Programming Language :: Python :: 3.4",
        "Programming Language :: Python :: 3.5",
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
}

if __name__ == "__main__":
    setup(**config)
