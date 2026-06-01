#!/usr/bin/env bash
# archive_check.sh -- validate an RdistNN_stata source archive.
#
# Usage:
#   bash tools/ci/archive_check.sh archives/source.zip

set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "$SCRIPT_DIR/../.." && pwd)"

if [[ $# -ne 1 ]]; then
	echo "usage: bash tools/ci/archive_check.sh <source.zip>" >&2
	exit 2
fi

ZIP_PATH="$1"
if [[ ! -f "$ZIP_PATH" && -f "$REPO_ROOT/$ZIP_PATH" ]]; then
	ZIP_PATH="$REPO_ROOT/$ZIP_PATH"
fi
if [[ ! -f "$ZIP_PATH" ]]; then
	echo "MISSING  archive not found: $1" >&2
	exit 1
fi

cd "$REPO_ROOT"

LIST_FILE="$(mktemp)"
trap 'rm -f "$LIST_FILE"' EXIT

if command -v unzip >/dev/null 2>&1; then
	unzip -Z1 "$ZIP_PATH" >"$LIST_FILE"
elif command -v python3 >/dev/null 2>&1; then
	python3 - "$ZIP_PATH" >"$LIST_FILE" <<'PY'
import sys
import zipfile

with zipfile.ZipFile(sys.argv[1]) as archive:
    for name in archive.namelist():
        print(name)
PY
elif command -v python >/dev/null 2>&1; then
	python - "$ZIP_PATH" >"$LIST_FILE" <<'PY'
import sys
import zipfile

with zipfile.ZipFile(sys.argv[1]) as archive:
    for name in archive.namelist():
        print(name)
PY
elif command -v powershell.exe >/dev/null 2>&1; then
	ZIP_PATH_FOR_PS="$ZIP_PATH" powershell.exe -NoProfile -Command \
		'Add-Type -AssemblyName System.IO.Compression.FileSystem; [IO.Compression.ZipFile]::OpenRead((Resolve-Path -LiteralPath $env:ZIP_PATH_FOR_PS).Path).Entries | ForEach-Object { $_.FullName }' \
		>"$LIST_FILE"
elif command -v powershell >/dev/null 2>&1; then
	ZIP_PATH_FOR_PS="$ZIP_PATH" powershell -NoProfile -Command \
		'Add-Type -AssemblyName System.IO.Compression.FileSystem; [IO.Compression.ZipFile]::OpenRead((Resolve-Path -LiteralPath $env:ZIP_PATH_FOR_PS).Path).Entries | ForEach-Object { $_.FullName }' \
		>"$LIST_FILE"
else
	tar -tf "$ZIP_PATH" >"$LIST_FILE"
fi

mapfile -t entries < <(sed 's#^\./##' "$LIST_FILE")
ERRORS=0

if [[ ${#entries[@]} -eq 0 ]]; then
	echo "archive-check FAILED: archive contains no readable entries."
	exit 1
fi

contains_entry() {
	local needle="$1"
	printf '%s\n' "${entries[@]}" | grep -Fx -- "$needle" >/dev/null
}

echo "=== Checking archive required files ==="
required=("README.md" "stata.toc" "rdistnn.pkg")
while IFS= read -r file; do
	required+=("$file")
done < <(awk '/^f[[:space:]]+/ {print $2}' rdistnn.pkg)
mapfile -t required < <(printf '%s\n' "${required[@]}" | awk '!seen[$0]++')

for file in "${required[@]}"; do
	if contains_entry "$file"; then
		echo "  OK       $file"
	else
		echo "  MISSING  $file"
		ERRORS=$((ERRORS + 1))
	fi
done

echo ""
echo "=== Checking archive forbidden paths ==="
forbidden=(".git/" ".github/" "archives/" "GH_logs/" "_dist/" "tests.log")
for pattern in "${forbidden[@]}"; do
	if printf '%s\n' "${entries[@]}" | grep -F -- "$pattern" >/dev/null; then
		echo "  FOUND    forbidden path pattern: $pattern"
		ERRORS=$((ERRORS + 1))
	else
		echo "  OK       absent: $pattern"
	fi
done

echo ""
echo "=== Checking archive duplicate entries ==="
duplicates="$(printf '%s\n' "${entries[@]}" | sort | uniq -d || true)"
if [[ -n "$duplicates" ]]; then
	echo "$duplicates" | sed 's/^/  DUP      /'
	ERRORS=$((ERRORS + 1))
else
	echo "  OK       no duplicate entries"
fi

echo ""
if [[ $ERRORS -gt 0 ]]; then
	echo "archive-check FAILED with $ERRORS error(s)."
	exit 1
fi

echo "archive-check passed."
