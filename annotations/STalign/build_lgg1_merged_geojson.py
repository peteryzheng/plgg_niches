#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path


def detect_workdir() -> str:
    home = os.path.expanduser("~")
    if home in ["/Users/youyun", "/Users/youyunzheng"]:
        return home + "/Documents/HMS/PhD/beroukhimlab/dfci_mount/"
    if home == "/home/yz762":
        return "/mnt/storage/dept/medonc/beroukhim/"
    if home == "/PHShome/yz762":
        return "/data/beroukhim1/"
    return "/data/beroukhim1/"


def load_feature_collection(path: Path) -> dict:
    with open(path, "r", encoding="utf-8") as handle:
        obj = json.load(handle)

    if not isinstance(obj, dict):
        raise ValueError(f"{path} is not a JSON object")
    if obj.get("type") != "FeatureCollection":
        raise ValueError(f"{path} is not a FeatureCollection")
    features = obj.get("features")
    if not isinstance(features, list):
        raise ValueError(f"{path} has no list-valued 'features'")
    return obj


def is_classified(feature: dict) -> bool:
    props = feature.get("properties") if isinstance(feature, dict) else None
    cls = props.get("classification") if isinstance(props, dict) else None
    return cls not in (None, "", {}, [])


def classified_counts(features: list[dict]) -> tuple[int, int]:
    with_cls = sum(1 for f in features if is_classified(f))
    without_cls = len(features) - with_cls
    return with_cls, without_cls


def main() -> None:
    workdir = Path(detect_workdir())
    geojson_dir = workdir / "youyun/plgg/data/Xenium_annotations/geojsons"

    parser = argparse.ArgumentParser(
        description="Merge LGG1 FeatureCollections for STalign notebooks."
    )
    parser.add_argument(
        "--input-a",
        type=Path,
        default=geojson_dir / "230918_Xenium_CytAssist_LGG1.geojson",
        help="First input FeatureCollection",
    )
    parser.add_argument(
        "--input-b",
        type=Path,
        default=geojson_dir / "230918_Xenium_CytAssist_LGG1_2_sample1_fc.geojson",
        help="Second input FeatureCollection",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=geojson_dir / "230918_Xenium_CytAssist_LGG1_merged.geojson",
        help="Output merged FeatureCollection",
    )
    parser.add_argument(
        "--drop-unclassified",
        action="store_true",
        help="Drop features lacking properties.classification after merging.",
    )
    args = parser.parse_args()

    fc_a = load_feature_collection(args.input_a)
    fc_b = load_feature_collection(args.input_b)
    a_features = fc_a["features"]
    b_features = fc_b["features"]

    merged = {
        "type": "FeatureCollection",
        "features": a_features + b_features,
    }

    merged_pre = merged["features"]
    if args.drop_unclassified:
        merged["features"] = [f for f in merged_pre if is_classified(f)]

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w", encoding="utf-8") as handle:
        json.dump(merged, handle, indent=2)
        handle.write("\n")

    a_with, a_without = classified_counts(a_features)
    b_with, b_without = classified_counts(b_features)
    pre_with, pre_without = classified_counts(merged_pre)
    post_with, post_without = classified_counts(merged["features"])

    print(
        f"input_a={args.input_a} total={len(a_features)} "
        f"classified={a_with} unclassified={a_without}"
    )
    print(
        f"input_b={args.input_b} total={len(b_features)} "
        f"classified={b_with} unclassified={b_without}"
    )
    print(
        f"merged_pre total={len(merged_pre)} "
        f"classified={pre_with} unclassified={pre_without}"
    )
    print(
        f"merged_post total={len(merged['features'])} "
        f"classified={post_with} unclassified={post_without}"
    )
    print(f"output={args.output} drop_unclassified={args.drop_unclassified}")


if __name__ == "__main__":
    main()
