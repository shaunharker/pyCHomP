# pyCHomP

CHomP (Computational Homology Project) with Python bindings

## Installation

```bash
pip install pychomp
```

## Troubleshooting / Installing from source

To reinstall:

```bash
pip uninstall pychomp
pip --no-cache-dir install pychomp -v -v -v
```

To build from sources:

```bash
git clone https://github.com/shaunharker/pyCHomP.git
cd pyCHomP
git submodule update --init --recursive
pip install . --ignore-installed --no-cache-dir -v -v -v
```
