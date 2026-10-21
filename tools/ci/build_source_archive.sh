#!/usr/bin/env bash
# Build the installable source ZIP from public inputs, excluding ignored output.
# Usage: bash tools/ci/build_source_archive.sh [_dist/source.zip]
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "$SCRIPT_DIR/../.." && pwd)"
if [[ $# -gt 1 ]]; then
    echo "usage: bash tools/ci/build_source_archive.sh [source.zip]" >&2
    exit 2
fi

SOURCE_ZIP="${1:-_dist/source.zip}"
cd "$REPO_ROOT"
python3 - "$SOURCE_ZIP" <<'PY'
import subprocess
import sys
import zipfile
from pathlib import Path

archive = Path(sys.argv[1]).resolve()
roots = ["README.md", "RELEASE_NOTES.md", "LICENSE.md", "LINKING_EXCEPTION.md", "rdistnn.pkg", "stata.toc"]
public = subprocess.check_output([
    "git", "ls-files", "--cached", "--others", "--exclude-standard", "-z",
    "--", "src/", "docs/", "examples/",
]).decode().split("\0")
deleted = set(subprocess.check_output([
    "git", "ls-files", "--deleted", "-z", "--", "src/", "docs/", "examples/",
]).decode().split("\0"))
files = sorted(set(roots + [name for name in public if name and name not in deleted]))
for name in files:
    if not Path(name).is_file():
        raise SystemExit("Source input is missing: " + name)
    if Path(name).resolve() == archive:
        raise SystemExit("Archive output must not replace a source input.")
archive.parent.mkdir(parents=True, exist_ok=True)
with zipfile.ZipFile(archive, "w", compression=zipfile.ZIP_DEFLATED) as source:
    for name in files:
        source.write(name, arcname=name)
print(str(archive) + ": " + str(len(files)) + " source files")
PY
bash "$SCRIPT_DIR/archive_check.sh" "$SOURCE_ZIP"
