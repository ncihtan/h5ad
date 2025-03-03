import anndata


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

    #def _check_genes(self, gene_list):
    #    """Check Genes."""
    #    for gene in gene_list:
    #        if not gene.startswith("ENSG"):
    #            self.error_list.append("Invalid gene identifier:  " + gene)

    def check_donor_ids(self, obs):
        """Check Donor IDs."""
        if "donor_id" in obs:
            donor_list = obs.donor_id.to_list()
            for donor_id in donor_list:
                # TODO:  Add RegEx Matching
                if donor_id.startswith("HTA"):
                    pass
                else:
                    self.error_list.append("Invalid donor_id:  " + donor_id)
        else:
            self.error_list.append("donor_id was not found in obs")

    def check_sample_ids(self, obs):
        """Check Sample IDs."""
        if "sample_id" in obs:
            sample_list = obs.sample_id.to_list()
            for sample_id in sample_list:
                # TODO:  Add RegEx Matching
                if sample_id.startswith("HTA"):
                    pass
                else:
                    self.error_list.append("Invalid sample_id:  " + sample_id)
        else:
            self.error_list.append("sample_id was not found in obs")
