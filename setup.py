#!/usr/bin/env python
# coding: utf-8

from setuptools import setup, find_packages


with open("README.md", "r") as fh:
    long_description = fh.read()

setup_args = dict(
    name="abipy_panel",
    description="Panels/Widgets for abipy",
    long_description=long_description,
    long_description_content_type="text/markdown",
    version="v2021.05.0",
    author="Adam Fekete",
    author_email="adam@fekete.co.uk",
    url="https://github.com/fekad/abipy-panel",
    license="BSD",
    platforms="Linux, Mac OS X, Windows",
    keywords=["Jupyter", "Widgets", "IPython"],
    classifiers=[
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        "Framework :: Jupyter",
    ],
    packages=find_packages(),
    python_requires=">=3.7",
    install_requires=[
        "abipy" "ipywidgets>=7.0.0",
    ],
    extras_require={
        "docs": ["mkdocs", "mkdocs-material"],
    },
)

if __name__ == "__main__":
    setup(**setup_args)