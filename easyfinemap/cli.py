"""Console script for easy_finemap."""

import typer

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
app = typer.Typer(context_settings=CONTEXT_SETTINGS, add_completion=False)


@app.command()
def main():
    """Main entrypoint."""
    typer.echo("easy_finemap")
    typer.echo("=" * len("easy_finemap"))
    typer.echo("user-friendly pipeline for GWAS fine-mapping")


if __name__ == "__main__":
    app()
