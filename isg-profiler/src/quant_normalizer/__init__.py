# SPDX-License-Identifier: GPL-3.0-only
# SPDX-FileCopyrightText: Copyright 2026 Luca Nishimura & Jumpei Ito

"""Quant Normalizer: Calculate normalization count form salmon quant.sf files."""

import importlib.metadata

try:
    __version__ = importlib.metadata.version("isg_profiler")
except importlib.metadata.PackageNotFoundError:
    __version__ = "unknown"
