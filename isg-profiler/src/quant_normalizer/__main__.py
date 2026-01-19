# SPDX-License-Identifier: GPL-3.0-only
# SPDX-FileCopyrightText: Copyright 2026 Luca Nishimura & Jumpei Ito

"""Entry point for isg-profiler command-line python part."""

import sys


def main():
    """Run the ISG-Profiler python CLI.

    This invokes the main function from cli.py.
    """
    try:
        from .cli import main as run_pipeline
    except ImportError as e:
        print("Error: Failed to import 'cli'.", file=sys.stderr)
        print(f"Current sys.path: {sys.path}", file=sys.stderr)
        raise e

    run_pipeline()


if __name__ == "__main__":
    main()
