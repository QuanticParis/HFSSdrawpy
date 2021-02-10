[![Python version](https://img.shields.io/badge/python-3.6%2B-blue)](https://www.python.org/downloads/release/python-360/)
[![License](https://img.shields.io/badge/license-MIT-yellow)](https://github.com/QuanticParis/HFSSdrawpy/blob/refactoring-2021/LICENSE)

# HFSSdrawpy

This package allows to draw circuits in HFSS using Python. It also directly 
generates GDS files.

## Installation
Install the package from pip using:
```shell
pip install HFSSdrawpy
```

## Development
### Requirements
[![Python version](https://img.shields.io/badge/python-3.6%2B-blue)](https://www.python.org/downloads/release/python-360/)

The project was written using **Python 3.6+**, you need to have a 
compatible Python version (>= 3.6) installed on your computer.

This project uses [Poetry](https://python-poetry.org/), please follow the 
instructions [here](https://python-poetry.org/docs/#installation) to install
it.

### Setup
Clone the repository and dive into it:
```shell
git clone https://github.com/QuanticParis/HFSSdrawpy.git
cd HFSSdrawpy
```

To setup the project run in the project directory:
```shell
poetry install
```

:point_up: On Windows you need to manually install `gdspy` from the 
pre-compiled binaries:
1. Head to [`gdspy` release page](https://github.com/heitzmann/gdspy/releases).
2. Download the release corresponding to your python and windows version 
   (latest: `gdspy-1.5.2-cp38-cp38-win_amd64.whl`).
3. Head to the `HFSSdrawpy` repository directory and run the following to 
   install `gdspy` from the pre-compiled binaries in your project Poetry 
   virtualenv:
   ```shell
   # to install the gdspy-1.5.2-cp38-cp38-win_amd64.whl gdspy release
   poetry run pip install path/to/downloaded/gdspy/release/gdspy-1.5.2-cp38-cp38-win_amd64.whl
   ```
