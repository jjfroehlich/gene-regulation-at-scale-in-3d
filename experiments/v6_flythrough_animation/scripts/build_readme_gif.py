#!/usr/bin/env python3
"""Create a small README GIF fallback from existing tracked V6 stills."""

from __future__ import annotations

import argparse
from pathlib import Path

from PIL import Image, ImageDraw


ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUTPUT = ROOT / "docs" / "images" / "v6-flythrough-preview.gif"
DEFAULT_INPUTS = [
    ROOT / "docs" / "images" / "v6-overview.jpg",
    ROOT / "docs" / "images" / "v6-p53-dna.jpg",
    ROOT / "docs" / "images" / "v6-transcription-end.jpg",
    ROOT / "docs" / "images" / "v6-nucleosome-loop.jpg",
    ROOT / "docs" / "images" / "v6-translation.jpg",
    ROOT / "docs" / "images" / "v6-actin.jpg",
]
CAPTIONS = [
    "V6 flythrough: annotated canonical overview",
    "p53 tetramer bound to DNA",
    "RNA polymerase II at the gene end and nascent RNA 3′ end",
    "Nucleosome loop on ACTB DNA",
    "Ribosome and tRNA translation machinery",
    "ACTB protein endpoint",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--width", type=int, default=720)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    frames = []
    for path, caption in zip(DEFAULT_INPUTS, CAPTIONS):
        image = Image.open(path).convert("RGB")
        height = round(image.height * args.width / image.width)
        image = image.resize((args.width, height), Image.Resampling.LANCZOS)
        canvas = Image.new("RGB", (args.width, height + 54), "white")
        canvas.paste(image, (0, 0))
        draw = ImageDraw.Draw(canvas)
        draw.rectangle((0, height, args.width, height + 54), fill=(245, 247, 248))
        draw.text((18, height + 18), caption, fill=(18, 22, 26))
        frames.append(canvas)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    frames[0].save(
        args.output,
        save_all=True,
        append_images=frames[1:],
        duration=900,
        loop=0,
        optimize=True,
    )
    print(args.output)


if __name__ == "__main__":
    main()
