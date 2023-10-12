# easyfinemap


[![pypi](https://img.shields.io/pypi/v/easyfinemap.svg)](https://pypi.org/project/easyfinemap/)
[![python](https://img.shields.io/pypi/pyversions/easyfinemap.svg)](https://pypi.org/project/easyfinemap/)
[![Build Status](https://github.com/Jianhua-Wang/easyfinemap/actions/workflows/dev.yml/badge.svg)](https://github.com/Jianhua-Wang/easyfinemap/actions/workflows/dev.yml)
[![codecov](https://codecov.io/gh/Jianhua-Wang/easyfinemap/branch/main/graphs/badge.svg)](https://codecov.io/github/Jianhua-Wang/easyfinemap)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![PyPI download month](https://img.shields.io/pypi/dm/easyfinemap.svg)](https://pypi.org/project/easyfinemap/)
[![Build Status](https://github.com/Jianhua-Wang/easyfinemap/actions/workflows/python-package-conda.yml/badge.svg)](https://github.com/Jianhua-Wang/easyfinemap/actions/workflows/python-package-conda.yml)


A flexible framework for GWAS fine-mapping


* Documentation: <https://Jianhua-Wang.github.io/easyfinemap>
* GitHub: <https://github.com/Jianhua-Wang/easyfinemap>
* PyPI: <https://pypi.org/project/easyfinemap/>
* License: MIT


## Features
* Formatting of summary statistics using smunger.
* Fast extraction of summary statistics for fine-mapping using tabix/bgzip.
* Checking and formatting of LD reference.
* Representation of variants using Unique SNP ID.
* Support for identifying independent loci using three methods: distance, LD clumping, and conditional analysis.
* Fine-mapping without the need for LD reference.
* Support for four LD-based fine-mapping tools and function-based fine-mapping.
* Fine-mapping combined with conditional analysis.

<!-- ## Finemapping approaches
* LD-free
    * aBF
* LD-based
    * FINEMAP
    * CAVIARBF
    * PAINTOR -->