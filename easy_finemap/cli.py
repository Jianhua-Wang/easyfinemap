"""Console script for easy_finemap."""

import typer


def main():
    """Main entrypoint."""
    typer.echo("easy_finemap")
    typer.echo("=" * len("easy_finemap"))
    typer.echo("user-friendly pipeline for GWAS fine-mapping")


if __name__ == "__main__":
    typer.run(main)
