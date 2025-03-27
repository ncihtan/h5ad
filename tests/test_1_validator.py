import anndata
import numpy as np
from h5ad_validate.validator import Validator


def test_validator1():
    """Test Gene Validator with missing HTAN specific obs"""
    # Create a matrix for 5 cells x 5 genes
    x = np.random.poisson(1, size=(5, 5))
    adata = anndata.AnnData(x)

    # Add Fake Cell Bar Codes
    adata.obs_names = [f"cell_bar_code{i:d}" for i in range(adata.n_obs)]

    # Ensembl IDs
    human_ensembl_ids = [
        "ENSG00000139618",
        "ENSG00000141510",
        "ENSG00000157764",
        "ENSG00000157761",
        "EGFR",
    ]

    adata.var_names = human_ensembl_ids

    # store adata for cellxgene check
    test_name = "tests/test.h5ad"
    adata.write_h5ad(test_name)

    validator = Validator(adata, test_name, "tests/test_out.txt")

    error_list = set(validator.error_list)
    donor_error = "donor_id was not found in obs"
    sample_error = "sample_id was not found in obs"
    cell_enrich_error = "cell_enrichment was not found in obs"
    intron_inclusion_error = "intron_inclusion was not found in obs"

    pass_code = validator.pass_code

    assert len(error_list) == 4
    assert donor_error in error_list
    assert sample_error in error_list
    assert cell_enrich_error in error_list
    assert intron_inclusion_error in error_list
    assert pass_code == [1, 1]
