"""
HTAN h5ad Validator.
"""
from htan.validator import Validator
import anndata
import click


@click.command()
@click.argument("h5ad_path")
def validate(h5ad_path):
    """HTAN h5ad Validator."""
    # TODO:  We should first run the built in CellXGene Validator, e.g via os()
    click.echo(click.style("HTAN h5ad File Validator", bg="blue", fg="white"))
    click.echo(click.style("File:  " + h5ad_path, fg="green"))
    adata = anndata.read_h5ad(h5ad_path)
    validator = Validator(adata)
    error_list = validator.error_list
    if len(error_list) == 0:
        print("Validation Passed!")
    else:
        print("Validation Failed")
        for error in error_list:
            print(error)


if __name__ == "__main__":
    validate()
