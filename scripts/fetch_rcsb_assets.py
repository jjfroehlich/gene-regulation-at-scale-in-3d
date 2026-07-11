#!/usr/bin/env python3
"""Download RCSB mmCIF assets used by canonical or retained experiments."""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import urllib.error
import urllib.request
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "config" / "scene_manifest.json"
RCSB_DIR = ROOT / "assets" / "rcsb"
MCP_DIR = ROOT / "assets" / "blender_mcp"
REPORT_PATH = ROOT / "outputs" / "canonical" / "asset_download_report.json"
MCP_ADDON_URL = "https://raw.githubusercontent.com/ahujasid/blender-mcp/main/addon.py"


def load_manifest() -> dict:
    return json.loads(MANIFEST.read_text(encoding="utf-8"))


def download(url: str, path: Path, force: bool = False) -> dict:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists() and path.stat().st_size > 0 and not force:
        data = path.read_bytes()
        return {
            "url": url,
            "path": str(path),
            "bytes": len(data),
            "sha256": hashlib.sha256(data).hexdigest(),
            "status": "cached",
        }
    req = urllib.request.Request(url, headers={"User-Agent": "gene-expression-scene-builder/1.0"})
    with urllib.request.urlopen(req, timeout=120) as response:
        data = response.read()
    path.write_bytes(data)
    return {
        "url": url,
        "path": str(path),
        "bytes": len(data),
        "sha256": hashlib.sha256(data).hexdigest(),
        "status": "downloaded",
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--force", action="store_true", help="Redownload existing files.")
    parser.add_argument("--skip-mcp", action="store_true", help="Skip Blender MCP add-on download.")
    parser.add_argument(
        "--pdb-id",
        action="append",
        default=[],
        help="Download only the specified PDB ID; repeat for multiple calibrators.",
    )
    args = parser.parse_args()

    manifest = load_manifest()
    pdb_ids = {pdb_id.upper() for pdb_id in args.pdb_id}
    if not pdb_ids:
        pdb_ids = {asset["pdb_id"].upper() for asset in manifest["pdb_assets"]}
        for group in manifest.get("nucleic_acid_calibrators", {}).values():
            for asset in group:
                pdb_ids.add(asset["pdb_id"].upper())
    pdb_ids = sorted(pdb_ids)
    records = []
    errors = []
    for pdb_id in pdb_ids:
        url = f"https://files.rcsb.org/download/{pdb_id}.cif"
        path = RCSB_DIR / f"{pdb_id}.cif"
        try:
            records.append(download(url, path, force=args.force))
        except (OSError, urllib.error.URLError, urllib.error.HTTPError) as exc:
            errors.append({"pdb_id": pdb_id, "url": url, "error": str(exc)})

    if not args.skip_mcp and not args.pdb_id:
        try:
            records.append(download(MCP_ADDON_URL, MCP_DIR / "addon.py", force=args.force))
        except (OSError, urllib.error.URLError, urllib.error.HTTPError) as exc:
            errors.append({"asset": "blender_mcp_addon", "url": MCP_ADDON_URL, "error": str(exc)})

    report = {"records": records, "errors": errors}
    REPORT_PATH.parent.mkdir(parents=True, exist_ok=True)
    REPORT_PATH.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))
    return 1 if errors else 0


if __name__ == "__main__":
    raise SystemExit(main())
