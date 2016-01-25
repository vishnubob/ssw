#!/usr/bin/env python

from setuptools import setup, Extension

libssw_ext = {"sources": ["src/ssw/ssw.c"], "include_dirs": ["src/ssw"]}

config = {
    "name": "ssw", 
    "version": "0.1",
    "description": "Complete Striped Smith Waterman Library",
    "author": "Mengyao Zhao et al.",
    "author_email": "zhaomengyao@gmail.com",
    "package_dir": {"ssw": "src"},
    "packages": ["ssw"],
    "ext_modules": [Extension("libssw", **libssw_ext)],
}

if __name__ == "__main__":
    setup(**config)
