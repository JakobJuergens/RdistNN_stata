#!/usr/bin/env bash
# Build the named source download and a checksum verifiable after downloading.
# Usage: bash tools/ci/build_release_assets.sh [_dist]
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "$SCRIPT_DIR/../.." && pwd)"
if [[ $# -gt 1 ]]; then
    echo "usage: bash tools/ci/build_release_assets.sh [release-directory]" >&2
    exit 2
fi

RELEASE_DIR="${1:-$REPO_ROOT/_dist}"
mkdir -p -- "$RELEASE_DIR"
RELEASE_DIR="$(cd -- "$RELEASE_DIR" && pwd)"
bash "$SCRIPT_DIR/build_source_archive.sh" "$RELEASE_DIR/rdistnn-source.zip"
(
    cd -- "$RELEASE_DIR"
    sha256sum rdistnn-source.zip > SHA256SUMS-source.txt
    sha256sum --check SHA256SUMS-source.txt
)
