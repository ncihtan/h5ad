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
    click.echo(click.style("HTAN h5ad File Validator", bg="blue", fg="white"))
    click.echo(click.style("File:  " + h5ad_path, fg="green"))
   
    try: 
        adata = anndata.read_h5ad(h5ad_path)
    except:
        click.echo(click.style("Unable to open File:  " + h5ad_path, fg="red"))
        click.echo(click.style("Please check file path and file integrity.", fg="red"))
        return
    
    validator = Validator(adata, h5ad_path)
    error_list = validator.error_list
    pass_code = validator.pass_code
    if len(error_list) == 0:
        click.echo(click.style("Validation Passed!", 
                           fg="green"))
    else:
        click.echo(click.style("Validation Failed.", 
                           fg="red"))
        ## to do -- create error log
        for error in error_list:
            print(error)


if __name__ == "__main__":
    validate()
