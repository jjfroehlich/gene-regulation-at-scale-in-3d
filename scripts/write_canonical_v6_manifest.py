#!/usr/bin/env python3
"""Write the resolved canonical V6 scene manifest."""

from __future__ import annotations

import json

from canonical_v6_config import MANIFEST_PATH, write_manifest


def main() -> None:
    manifest = write_manifest()
    print(json.dumps({"manifest": str(MANIFEST_PATH), "title": manifest["title"]}, indent=2))


if __name__ == "__main__":
    main()

