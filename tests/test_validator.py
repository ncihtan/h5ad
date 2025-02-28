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
        "HTA200_1_1",
        "HTA200_1_1",
        "HTA200_1_1",
        "HTA200_1_1",
        "UUID_12321",
    ]
    adata.var_names = human_ensembl_ids
    adata.obs["donor_id"] = donor_ids
    adata.obs["sample_id"] = sample_ids
    validator = Validator(adata)

    error_list = validator.error_list
    assert len(error_list) == 3
    assert error_list[0] == "Invalid gene identifier:  EGFR"
    assert error_list[1] == "Invalid donor_id:  UUID_12321"
    assert error_list[2] == "Invalid sample_id:  UUID_12321"
