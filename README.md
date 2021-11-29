# HFSSdrawpy

This package allows to draw circuits in HFSS using Python. It also directly generates GDS files.

[Source Code](https://github.com/QuanticParis/HFSSdrawpy/)

## Installation

Installing using pip will soon be available. Meanwhile download the [source code](https://github.com/QuanticParis/HFSSdrawpy/).

### Linux / MacOS

Go to the HFSSdrawpy directory and run:

```
python setup.py develop
```

This kind of install enables the latest changes from Github conveniently. In the near future
when the code will be more stable we will organize a release. From then, the install wil be done by running

```
pip install HFSSdrawpy
```

### Windows

You need to manually install gdspy from the pre-compiled binaries, [here](https://github.com/heitzmann/gdspy/releases). Download the one corresponding to your python and windows version (latest: gdspy-1.5.2-cp38-cp38-win_amd64.whl). At the downloaded file location run:

```
pip install gdspy-1.5.2-cp38-cp38-win_amd64.whl
```

with the relevant name.

Then proceed as in the Linux / MacOs install.

### Documentation

Install Sphinx with

```
pip install Sphinx
```

and build the documentation.

```
cd doc
make html
```
Enjoy!

`./doc/build/html/index.html`