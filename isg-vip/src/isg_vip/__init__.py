# SPDX-License-Identifier: GPL-3.0-only
# SPDX-FileCopyrightText: Copyright 2026 Hiroaki Unno & Jumpei Ito

"""ISG-VIP: ISG-based Viral Infection Predictor.

Machine learning-based viral infection prediction using ISG (Interferon-Stimulated Gene)
expression profiles with ensemble methods.
"""

import importlib.metadata

try:
    __version__ = importlib.metadata.version("isg_vip")
except importlib.metadata.PackageNotFoundError:
    __version__ = "unknown"
