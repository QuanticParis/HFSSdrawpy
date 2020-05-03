# HFSSdrawpy

This package allows to draw circuits in HFSS using Python. It also directly generates GDS files.

[Source Code](https://github.com/leghtas/HFSSdrawpy/)

## Installation

Installing using pip will soon be available. Meanwhile download the [source code](https://github.com/leghtas/HFSSdrawpy/).

### Linux / MacOS

Go to the HFSSdrawpy directory and run: 

`python setup.py install`

To easily get the latest changes from git, run instead:

`python setup.py develop`

### Windows

You need to manually install gdspy from the pre-compiled binaries, [here](https://github.com/heitzmann/gdspy/releases). Download the one corresponding to your python and windows version (latest: gdspy-1.5.2-cp38-cp38-win_amd64.whl). At the downloaded file location run:

`pip install gdspy-1.5.2-cp38-cp38-win_amd64.whl`

with the relevant name.

Then proceed as in the Linux / MacOs install.

