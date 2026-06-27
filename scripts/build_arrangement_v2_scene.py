#!/usr/bin/env python3
"""Build the arrangement V2 protein-placement experiment scene."""

from __future__ import annotations

import os
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
os.environ["ARRANGEMENT_CONFIG_MODULE"] = "arrangement_v2_config"

from build_arrangement_v1_scene import main  # noqa: E402


if __name__ == "__main__":
    main()
