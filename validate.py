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

    try: 
        adata = anndata.read_h5ad(h5ad_path)
    except:
        click.echo(click.style("Unable to open File:  " + h5ad_path, fg="red"))
        click.echo(click.style("Please check file path and file integrity.", fg="red"))
        return
    
    os.environ['H5AD_PATH'] = h5ad_path
    cellxgene = os.WEXITSTATUS(os.system("cellxgene-schema validate $H5AD_PATH"))  
    if cellxgene !=0:
        error_list = ['cellxgene-schema error']
        click.echo(click.style("Cellxgene run has errors. Please note errors or warnings in the output above this line.", 
                           fg="red"))
    else: 
        error_list = []
        click.echo(click.style("Cellxgene run successful. Please note any warnings in the output above this line.", 
                           fg="green"))

    
    click.echo(click.style("Running HTAN specific checks.", fg="green"))
    validator = Validator(adata)
    error_list.extend(validator.error_list)
    if len(error_list) == 0:
        print("Validation Passed!")
    else:
        print("Validation Failed")
        ## to do -- create error log
        for error in error_list:
            print(error)


if __name__ == "__main__":
    validate()
