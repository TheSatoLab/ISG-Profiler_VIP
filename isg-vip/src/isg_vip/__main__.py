# SPDX-License-Identifier: GPL-3.0-only
# SPDX-FileCopyrightText: Copyright 2026 Hiroaki Unno & Jumpei Ito

"""Entry point for isg-vip command-line tool.

This module provides a thin wrapper around the existing prediction script,
allowing the tool to be run as a command after installation with pip.
"""

import sys


def main():
    """Run the ISG-VIP prediction pipeline.

    This invokes the main function from pipelines.py.
    """
    try:
        from .pipelines import main as run_pipeline
    except ImportError as e:
        print("Error: Failed to import 'cli'.", file=sys.stderr)
        print(f"Current sys.path: {sys.path}", file=sys.stderr)
        raise e

    run_pipeline()


if __name__ == "__main__":
    main()
