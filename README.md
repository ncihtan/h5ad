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
python validate.py example/HTAN_h5ad_exemplar_2025_03_03.h5ad output_file.txt
```

To run the unit tests, run:

```commandline
pytest tests
```

## Validate.py

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
2. Evaluate the self.pass_code to produce final message (Validation Passed) or failure messsage.
3. Print all function warnings and errors to the provided output file.

## Specific Validation Functions

### CELLxGENE Validation
- self.check_cell_x_gene(h5ad_path, output_file)
    - runs cellxgene-schema.
    - writes warnings and error messages to the provided output_file.
    - writes cellxgene-schema pass or fail message to screen.
    - if any errors, pass_code[0] is set to 1

### HTAN-specific Validation
- self.check_donor_ids(adata.obs)
    - checks for presence of obs.donor_id.
    - checks all unique values in obs.donor_id for match to r"^(HTA20[0-9])(?:_0000)?(?:_\d+)?(?:_EXT\d+)?"
    - if any errors, pass_code[1] is set to 1 and error strings added to error_list.
 - self.check_sample_ids(adata.obs)
    - checks for presence of obs.sample_id.
    - checks all unique values in obs.sample_id for match to r"^(HTA20[0-9])(?:_0000)?(?:_\d+)?(?:_EXT\d+)?_(B|D)\d{1,50}$"
    - if any errors, pass_code[1] is set to 1 and error strings added to error_list.
- self.check_cell_enrichment(adata.obs)
    - checks for presence of obs.cell_enrichment.
    - checks that all unique values in obs.cell_enrichment for match to r"^CL:(00000000|[0-9]{7}[+-])$"
    - strips + or - character from the cell_enrichment term to check if the CL term is valid.
        - valid CL terms taken from file htan/CL_codes_human.tsv
        - CL:00000000 (no enrichment) added to htan/CL_codes_human.tsv before term checked.
    - if any errors, pass_code[1] is set to 1 and error strings added to error_list.
- self.check_intron_inclusion(adata.obs)
    - checks for presence of obs.intron_inclusion.
    - verifies all unique values in obs.intron_inclusion are in ['yes','no']
    - if any errors, pass_code[1] is set to 1 and error strings added to error_list.
