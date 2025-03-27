# HTAN h5ad Validator

This is an (AnnData 0.10) h5ad validator for HTAN Phase 2 single cell/single nuclei RNA-sequencing data. It:
1. runs [cellxgene-schema validate](https://github.com/chanzuckerberg/single-cell-curation) to check if the h5ad file conforms to CELLxGENE's schema; and
2. Runs additional scripts to validate HTAN-specific requirements.

## Installation

1. It is recommended that you create and a virtual environment:

```commandline
python -m venv /path/to/new/virtual/environment
```
Then, activate the virtual environment:

```commandline
source /path/to/new/virtual/environment/bin/activate
```

2. Install the HTAN h5ad package using pip:

```commandline
pip install HTAN-h5ad-validate
```

3. To check your installation, run:

```commandline
HTAN-h5ad-validate example/HTAN_h5ad_exemplar_2025_03_03.h5ad output_file.txt
```
You should see the following messages on your screen in green:

```
File:  example/HTAN_h5ad_exemplar_2025_03_03.h5ad
Running cellxgene-schema
Cellxgene run successful. Please check the output file to see if warnings exist.
Running HTAN-specific validation
Validation Passed!
```
You should also see a new file called "output_file.txt" with the following messages:

```
cellxgene-schema output: Starting validation...
WARNING: Dataframe 'var' only has 12641 rows. Features SHOULD NOT be filtered from expression matrix.
Validation complete in 0:00:03.421082 with status is_valid=True
```

Please note this is a test h5ad. Your h5ad should not have filtered features!

## Use
To validate an h5ad file, run the following from the command line:

```commandline
HTAN-h5ad-validate {path/to/your/h5ad} {output_filename}.txt
```
replacing {path/to/your/h5ad} and {output_filename} with appropriate text strings. (Please see Installation instructions for an example.)

Information regarding whether the h5ad file passed or failed validation will be printed to the terminal window. Warnings and error messages will be printed to {output_filename}.txt.

If the h5ad file fails validation, please resolve the errors noted and retry validation.

## Information for Developers
To further develop this code, please use Python >= 3.8.

1. It is recommended that you create a virtual environment:

```commandline
python -m venv /path/to/new/virtual/environment
```
Then, activate the virtual environment:

```commandline
source /path/to/new/virtual/environment/bin/activate
```

2. Clone this github repository.

3. Run pip install locally with optional developer dependencies:

```commandline
pip install -e ".[dev]"
```

4. To run the validator from the command line, run:

- For successful validation:
```commandline
HTAN-h5ad-validate example/HTAN_h5ad_exemplar_2025_03_03.h5ad output_file.txt
```

- For errors in both CELLxGENE validation and HTAN validations:
```commandline
HTAN-h5ad-validate example/HTAN_h5ad_error_exemplar.h5ad output_file.txt
```

- For successful CELLxGENE, but HTAN validation errors, run:
```commandline
HTAN-h5ad-validate example/HTAN_h5ad_exemplar_HTANonly_error.h5ad output_file.txt
```

5. To run the unit tests, run:

```commandline
pytest tests
```
Please see the README.md in the tests folder for information about the tests.

### Validate.py

Provided a valid h5ad file path and an output filename, Validate.py will: 
1. Create a Validator object with:
    - self.error_list: a list of errors from HTAN specific validations.
    - self.pass_code: a two item list where 0 indicates pass and 1 indicates fail.
        - [0,0] = pass both cellxgene-schema and HTAN-specific validation.
        - [1,0] = fail cellxgene-schema, pass HTAN-specific validation.
        - [0,1] = pass cellxgene-schema, fail HTAN-specific validation.
        - [1,1] = fail both cellxgene-schema and HTAN-specific validation.
    - Validation functions (see descriptions below):
        - self.check_cell_x_gene(h5ad_path, output_file)
        - self.check_donor_ids(adata.obs)
        - self.check_sample_ids(adata.obs)
        - self.check_cell_enrichment(adata.obs)
        - self.check_intron_inclusion(adata.obs)
2. Evaluate self.pass_code to produce a final message (Validation Passed or HTAN validation failure).
3. Print all function warnings and errors to the provided output filename.

### Specific Validation Functions (htan/validator.py)

#### CELLxGENE Validation
- self.check_cell_x_gene(h5ad_path, output_file)
    - runs [cellxgene-schema](https://github.com/chanzuckerberg/single-cell-curation).
    - writes warnings and error messages to output_file.
    - writes cellxgene-schema pass or fail message to screen.
    - if any errors, pass_code[0] is set to 1

#### HTAN-specific Validation
- self.check_donor_ids(adata.obs)
    - checks for presence of obs.donor_id.
    - checks all unique values in obs.donor_id match r"^(HTA20[0-9])(?:_0000)?(?:_\d+)?(?:_EXT\d+)?"
    - if any errors, pass_code[1] is set to 1 and error strings added to error_list.
 - self.check_sample_ids(adata.obs)
    - checks for presence of obs.sample_id.
    - checks all unique values in obs.sample_id match r"^(HTA20[0-9])(?:_0000)?(?:_\d+)?(?:_EXT\d+)?_(B|D)\d{1,50}$"
    - if any errors, pass_code[1] is set to 1 and error strings added to error_list.
- self.check_cell_enrichment(adata.obs)
    - checks for presence of obs.cell_enrichment.
    - checks that all unique values in obs.cell_enrichment match r"^CL:(00000000|[0-9]{7}[+-])$"
    - strips + or - character from the cell_enrichment term to check if the CL term is valid.
        - valid CL terms taken from file htan/CL_codes_human.tsv
        - CL:00000000 (no enrichment) added to htan/CL_codes_human.tsv before term checked.
    - if any errors, pass_code[1] is set to 1 and error strings added to error_list.
- self.check_intron_inclusion(adata.obs)
    - checks for presence of obs.intron_inclusion.
    - verifies all unique values in obs.intron_inclusion are in ['yes','no']
    - if any errors, pass_code[1] is set to 1 and error strings added to error_list.
