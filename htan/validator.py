import anndata
import re


class Validator:
    """HTAN h5ad Validator."""

    def __init__(self, adata):
        """Constructor with anndata data structure."""
        self.adata = adata
        self.error_list = []
        #gene_list = self.adata.var.index.tolist()
        #self._check_genes(gene_list)
        self.check_donor_ids(adata.obs)
        self.check_sample_ids(adata.obs)
        self.check_cell_enrichment(adata.obs)
        self.check_intron_inclusion(adata.obs)

    def check_donor_ids(self, obs):
        """Check Donor IDs."""
        pattern = r"^(HTA20[0-9])(?:_0000)?(?:_\d+)?(?:_EXT\d+)?"
        if "donor_id" in obs:
            donor_list = list(obs.donor_id.unique())
            for donor_id in donor_list:
                # TODO:  Add RegEx Matching
                if re.match(pattern, donor_id):   
                    pass
                else:
                    self.error_list.append("Invalid donor_id:  " + donor_id)
        else:
            self.error_list.append("donor_id was not found in obs")

    def check_sample_ids(self, obs):
        """Check Sample IDs."""
        pattern = r"^(HTA20[0-9])(?:_0000)?(?:_\d+)?(?:_EXT\d+)?_(B|D)\d{1,50}$"
        if "sample_id" in obs:
            sample_list = list(obs.sample_id.unique())  
            for sample_id in sample_list:
                # TODO:  Add RegEx Matching
                if re.match(pattern, sample_id):
                    pass
                else:
                    self.error_list.append("Invalid sample_id:  " + sample_id)
        else:
            self.error_list.append("sample_id was not found in obs")

    def check_cell_enrichment(self, obs):
        """Check Cell Enrichment."""
        # POSSIBLE TO DO: add step to check for valid CL term 
        pattern = r"^CL:(00000000|[0-9]{7}[+-])$"
        if "cell_enrichment" in obs:
            cell_enrichment_list = list(obs.cell_enrichment.unique())
            for cell_enrich_term in cell_enrichment_list:
                if re.match(pattern, cell_enrich_term):
                    pass
                else:
                    self.error_list.append("Invalid cell_enrichment term " + cell_enrich_term)
        else:
            self.error_list.append("cell_enrichment was not found in obs")

    def check_intron_inclusion(self, obs):
        """Check intron inclusion"""
        valid_values = ['yes','no']
        intron_inclusion_list = list(obs.intron_inclusion.unique().astype(str))
        if "intron_inclusion" in obs:
            for intron_include_term in intron_inclusion_list:
                if intron_include_term in valid_values:
                    pass
                else:
                    self.error_list.append("Invalid intron_inclusion term: " + intron_include_term + ". Must be 'yes' or 'no'.")
        else:
            self.error_list.append("intron_inclusion was not found in obs")

