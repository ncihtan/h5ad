# tests

The following tests exist in this folder:

1. test_0_validator.py

    Simulation of incorrect values for donor_id, sample_id, cell_enrichment and intron_inclusion. Should also produce cellxgene-schema errors which are captured in test_out.txt.

    - creates a dummy adata object and saves it in h5ad format as "test.h5ad"
        - adata has errors in ensemble_id, donor_id, sample_id, cell_enrichment and intron_inclusion.
        - one error each for ensemble_id, donor_id and sample_id.
        - multiple cell_enrichment errors:
            - invalid format (CL term, missing + or -)
            - invalid CL term
            - invalid format + not CL term
        - multiple intron_inclusion errors:
            - blank 
            - true (instead of 'yes' or 'no')
            - false (instead of 'yes' or 'no')
    - runs htan/validator.py with arguments adata, test.h5ad and tests/test_out.txt.
    - asserts:
        - 8 total error messages;
        - 1 specific error message each for donor_id and sample_id;
        - 3 specific error messages each for cell_enrichment and intron_inclusion; and
        - final pass code is [1,1]. (failed both cellxgene-schema.py and HTAN-specific validations.)
    - does not assert specific cellxgene-schema errors, but these can be found in test_out.txt


2. test_1_validator.py

    Simulation of missing obs.donor_id, obs.sample_id, obs.cell_enrichment and obs.intron_inclusion. Should also produce cellxgene-schema errors which are captured in test_out.txt.

    - creates a dummy adata object and saves it in h5ad format as "test.h5ad"
        - adata has error in ensemble_id
        - adata missing obs.donor_id, obs.sample_id, obs.cell_enrichment and obs.intron_inclusion.
    - runs htan/validator.py with arguments adata, test.h5ad and tests/test_out.txt.
    - asserts:
        - 4 total error messages;
        - 4 specific error messages for missing obs.donor_id, obs.sample_id, obs.cell_enrichment and obs.intron_inclusion; and
        - final pass code is [1,1]. (failed both cellxgene-schema.py and HTAN-specific validations.)
    - does not assert specific cellxgene-schema errors, but these can be found in test_out.txt

3. test_2_validator.py

    Simulation of incorrect path for h5ad file and/or corrupt h5ad.

    - sets incorrect_path to nonsense path
    - sets non_h5ad to tests/not_real_h5ad.h5ad (a text file saved with .h5ad extension).
    - runs validator.py with incorrect path, asserts error exit code.
    - runs validator.py with non_h5ad, assert error exit code.