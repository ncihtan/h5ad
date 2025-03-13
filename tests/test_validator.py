import anndata
import numpy as np
from htan.validator import Validator


def test_validator0():
    """Test Gene Validator"""
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
    donor_ids = [
        "HTA200_1",
        "HTA200_1",
        "HTA200_1",
        "HTA200_1",
        "UUID_12321",
    ]
    sample_ids = [
        "HTA200_1_B1",
        "HTA200_1_B1",
        "HTA200_1_B1",
        "HTA200_1_D1",
        "UUID_12321",
    ]
    adata.var_names = human_ensembl_ids
    adata.obs["donor_id"] = donor_ids
    adata.obs["sample_id"] = sample_ids

    #store adata -- this will be replaced with permanent example
    test_name = "tests/test.h5ad"
    adata.write_h5ad(test_name)


    validator = Validator(adata, test_name, "tests/test_out.txt")

    error_list = set(validator.error_list)
    donor_error = "Invalid donor_id: UUID_12321"
    sample_error = "Invalid sample_id: UUID_12321"
    cell_enrich_error = "cell_enrichment was not found in obs"
    intron_inclusion_error = "intron_inclusion was not found in obs"

    assert len(error_list) == 4
    assert donor_error in error_list
    assert sample_error in error_list
    assert cell_enrich_error in error_list
    assert intron_inclusion_error in error_list
