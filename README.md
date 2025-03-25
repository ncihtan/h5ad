# HTAN h5ad Validator

This is an (AnnData 0.10) h5ad validator for HTAN Phase 2 single cell/single nuclei RNA-sequencing data. It:
1. runs [cellxgene-schema validate](https://github.com/chanzuckerberg/single-cell-curation) to check if the h5ad file conforms to CELLxGENE's schema; and
2. Runs additional scripts to validate HTAN-specific requirements.

## Installation

(Pip install capability and directions pending)

To check your installation, run:

```commandline
python validate.py example/HTAN_h5ad_exemplar_2025_03_03.h5ad output_file.txt
```

## Use
To validate an h5ad file, run the following from the command line:

```commandline
python validate.py {path/to/your/h5ad} {output_filename}.txt
```
replacing {path/to/your/h5ad} and {output_filename} with appropriate text strings. (Please see Installation instructions for an example.)

Information regarding whether the h5ad file passed or failed validation will be printed to the terminal window. Warnings and error messages will be printed to {output_filename}.txt.

If the h5ad file fails validation, please resolve the errors noted and retry validation.

## Information for Developers
To further develop this code, please use Python >= 3.8.

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
python validate.py example/HTAN_h5ad_exemplar_2025_03_03.h5ad output_file.txt
```

To run the unit tests, run:

```commandline
pytest tests
```
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
