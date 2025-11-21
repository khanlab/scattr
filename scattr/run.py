#!/usr/bin/env python3
import os
from pathlib import Path

from snakebids import bidsapp, plugins

try:
    from scattr.workflow.lib import (
        utils as utils,  # Works when run as a package
    )
except ImportError:
    from workflow.lib import utils as utils  # Works when run directly

if "__file__" not in globals():
    __file__ = "../scattr/run.py"

app = bidsapp.app(
    [
        plugins.SnakemakeBidsApp(Path(__file__).resolve().parent),
        plugins.BidsValidator(),
        plugins.CliConfig("parse_args"),
        plugins.ComponentEdit("pybids_inputs"),
    ]
)

# Set the conda prefix directory
conda_prefix = Path(utils.get_download_dir()) / "conda"

# Set the environment variable SNAKEMAKE_CONDA_PREFIX if not already set
if "SNAKEMAKE_CONDA_PREFIX" not in os.environ:
    os.environ["SNAKEMAKE_CONDA_PREFIX"] = str(conda_prefix)

def get_parser():
    """Exposes parser for sphinx doc generation, cwd is the docs dir"""
    return app.build_parser().parser


if __name__ == "__main__":
    app.run()