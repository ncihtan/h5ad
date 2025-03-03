"""
HTAN h5ad Validator.
"""
from htan.validator import Validator
import anndata
import click
import os


@click.command()
@click.argument("h5ad_path")
def validate(h5ad_path):
    """HTAN h5ad Validator."""
    # TODO:  We should first run the built in CellXGene Validator, e.g via os()
    click.echo(click.style("HTAN h5ad File Validator", bg="blue", fg="white"))
    click.echo(click.style("File:  " + h5ad_path, fg="green"))
    
    click.echo(click.style("Running cellxgene-schema", fg="green"))
    click.echo(click.style("The following output is from cellxgene-schema", fg="green"))

    os.environ['H5AD_PATH'] = h5ad_path
    os.system("cellxgene-schema validate $H5AD_PATH")  
    click.echo(click.style("Cellxgene run complete. Please note any errors or warnings in the output above this line.", 
                           fg="green"))
    
    try: 
        adata = anndata.read_h5ad(h5ad_path)
    except:
        click.echo(click.style("Unable to open File:  " + h5ad_path, fg="red"))
        return

    
    #### TO DO: need to not report validation passed if cellxgene Schema has error
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
