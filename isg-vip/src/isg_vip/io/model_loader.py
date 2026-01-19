# SPDX-License-Identifier: GPL-3.0-only
# SPDX-FileCopyrightText: Copyright 2026 Hiroaki Unno & Jumpei Ito

import logging
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Any, ClassVar, Dict, Tuple

import joblib
import numpy as np

logger = logging.getLogger(__name__)


class ModelType(str, Enum):
    LGB = "lgb"
    LR = "lr"
    META = "meta"


@dataclass(frozen=True)
class ISGModelArtifacts:
    """
    An immutable container holding loaded ISG model artifacts.
    """

    encoders: Dict[ModelType, Tuple[Any, ...]]
    final_models: Dict[ModelType, Tuple[Any, ...]]
    normalizers: Tuple[Any, ...]
    thresholds: Dict[ModelType, Any]

    _N_FOLDS: ClassVar[int] = 5

    @classmethod
    def from_directory(cls, model_dir: Path) -> "ISGModelArtifacts":
        temp_encoders = {t: [] for t in ModelType}
        temp_models = {t: [] for t in ModelType}
        temp_normalizers = []
        temp_thresholds = {}

        missing_files = []

        for fold in range(cls._N_FOLDS):
            # Normalizer
            norm_path = model_dir / f"normalizer_{fold}.pkl"
            if not norm_path.exists():
                missing_files.append(str(norm_path))
            else:
                temp_normalizers.append(joblib.load(norm_path))

            # Encoder & Final Model
            for m_type in ModelType:
                enc_path = model_dir / f"encoder_{m_type.value}_{fold}.pkl"
                mod_path = model_dir / f"final_model_{m_type.value}_{fold}.pkl"

                if not enc_path.exists():
                    missing_files.append(str(enc_path))
                else:
                    temp_encoders[m_type].append(joblib.load(enc_path))

                if not mod_path.exists():
                    missing_files.append(str(mod_path))
                else:
                    temp_models[m_type].append(joblib.load(mod_path))

        # read thresholds
        for m_type in ModelType:
            th_path = model_dir / f"thresholds_{m_type.value}.npy"
            if not th_path.exists():
                missing_files.append(str(th_path))
            else:
                temp_thresholds[m_type] = np.load(th_path)

        if missing_files:
            raise FileNotFoundError(
                f"Missing required model files ({len(missing_files)}):\n" + "\n".join(missing_files)
            )

        logger.info("Success to load models.")
        return cls(
            encoders={k: tuple(v) for k, v in temp_encoders.items()},
            final_models={k: tuple(v) for k, v in temp_models.items()},
            normalizers=tuple(temp_normalizers),
            thresholds=temp_thresholds,
        )

    def _validate_fold(self, fold: int):
        if not (0 <= fold < self._N_FOLDS):
            raise IndexError(f"Fold index out of range: {fold} (0 ~ {self._N_FOLDS - 1})")

    def get_encoder(self, model_type: ModelType, fold: int) -> Any:
        self._validate_fold(fold)
        return self.encoders[model_type][fold]

    def get_model(self, model_type: ModelType, fold: int) -> Any:
        self._validate_fold(fold)
        return self.final_models[model_type][fold]

    def get_normalizer(self, fold: int) -> Any:
        self._validate_fold(fold)
        return self.normalizers[fold]

    def get_threshold(self, model_type: ModelType) -> Any:
        return self.thresholds[model_type]
