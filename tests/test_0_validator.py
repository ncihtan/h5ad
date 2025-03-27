import anndata
import numpy as np
from h5ad_validate.validator import Validator


def test_validator0():
    """Test Gene Validator with incorrect ensembl_ids, donor_ids,
    samples_ids, cell_enrichment, intron_inclusion."""
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
    cell_enrichment = [
        "CL:0000182+",
        "CL:9999999+",
        "CL:00000000",
        "CL:0000182",
        "t-cells"
    ]
    intron_inclusion = [
        "",
        "no",
        "yes",
        "true",
        "false"
    ]

    adata.var_names = human_ensembl_ids
    adata.obs["donor_id"] = donor_ids
    adata.obs["sample_id"] = sample_ids
    adata.obs["cell_enrichment"] = cell_enrichment
    adata.obs["intron_inclusion"] = intron_inclusion

    # store adata for cellxgene check
    test_name = "tests/test.h5ad"
    adata.write_h5ad(test_name)

    validator = Validator(adata, test_name, "tests/test_out.txt")

    error_list = set(validator.error_list)
    donor_error = "Invalid donor_id: UUID_12321"
    sample_error = "Invalid sample_id: UUID_12321"
    cell_enrich_error1 = ("Invalid cell_enrichment term CL:9999999+. "
                          "CL_term is not in CL_codes_human.tsv")
    cell_enrich_error2 = ("Invalid cell_enrichment term CL:0000182. "
                          "obs.cell_enrichment must be "
                          "CL term followed by a '+' or '-' "
                          "sign or CL:00000000 if no enrichment.")
    cell_enrich_error3 = ("Invalid cell_enrichment term t-cells. "
                          "obs.cell_enrichment must be "
                          "CL term followed by a '+' or '-' "
                          "sign or CL:00000000 if no enrichment.")
    intron_inclusion_error1 = ("Invalid intron_inclusion term: true. "
                               "Must be 'yes' or 'no'.")
    intron_inclusion_error2 = ("Invalid intron_inclusion term: . "
                               "Must be 'yes' or 'no'.")
    intron_inclusion_error3 = ("Invalid intron_inclusion term: false. "
                               "Must be 'yes' or 'no'.")

    pass_code = validator.pass_code

    assert len(error_list) == 8
    assert donor_error in error_list
    assert sample_error in error_list
    assert cell_enrich_error1 in error_list
    assert cell_enrich_error2 in error_list
    assert cell_enrich_error3 in error_list
    assert intron_inclusion_error1 in error_list
    assert intron_inclusion_error2 in error_list
    assert intron_inclusion_error3 in error_list
    assert pass_code == [1, 1]
