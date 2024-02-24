#!/usr/bin/env python
from pathlib import Path

from snakebids import bidsapp, plugins

app = bidsapp.app(
    [
        plugins.SnakemakeBidsApp(Path(__file__).resolve().parent),
        plugins.Version(distribution="scattr"),
    ]
)


def get_parser():
    """Exposes parser for sphinx doc generation, cwd is the docs dir"""
    return app.build_parser().parser


def main():
    app.run()


if __name__ == "__main__":
    main()
