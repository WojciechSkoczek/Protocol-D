import argparse
import yaml
import os
import sys

parser = argparse.ArgumentParser()
parser.add_argument("--config", required=True)
args = parser.parse_args()

with open(args.config, "r", encoding="utf-8") as f:
    cfg = yaml.safe_load(f)

missing = []

for ds in cfg.get("datasets", []):
    if not ds.get("enabled", True):
        continue
    path = ds["input"]["path"]
    if not os.path.exists(path):
        missing.append(path)

if missing:
    print("Missing raw files:")
    for m in missing:
        print("  ", m)
    sys.exit(1)

print("All raw files present.")