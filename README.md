# HTAN h5ad Validator

This is a bare bones h5ad validator for HTAN.

To further develop this code, please use Python 3.6 or above.

Next, it is recommended that you create a virtual environment:

```commandline
python -m venv .venv
```

To activate the virtual environment:

```commandline
source .venv/bin/activate
```

Then, install all dependencies:

```commandline
pip install -r requirements.txt
```

To run the validator from the command line, run:

```commandline
python validate.py example/HTAN_h5ad_exemplar_2025_03_03.h5ad
```

To run the unit tests, run:

```commandline
pytest tests
```
