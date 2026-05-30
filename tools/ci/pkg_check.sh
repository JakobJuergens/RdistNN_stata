#!/usr/bin/env bash
# pkg_check.sh — validate RdistNN Stata package integrity
#
# Checks:
#   1. Every file listed under "f " in rdistnn.pkg exists on disk.
#   2. Every public .ado command in src/ has a matching .sthlp help file.
#   3. Every public .ado/.sthlp surface is listed in rdistnn.pkg.
#   4. The installable plugin is listed in rdistnn.pkg.
#   5. Active ado/test entry files target Stata 17+.
#   6. stata.toc references rdistnn.pkg.
#
# Exit code: 0 = all checks passed, 1 = one or more failures.
# Run from the repository root.

set -euo pipefail

ERRORS=0

# ── 1. rdistnn.pkg file manifest ────────────────────────────────────────────
echo "=== Checking rdistnn.pkg file entries ==="
while IFS= read -r line; do
	if [[ "$line" =~ ^f[[:space:]]+(.+)$ ]]; then
		filepath="${BASH_REMATCH[1]}"
		if [[ -f "$filepath" ]]; then
			echo "  OK       $filepath"
		else
			echo "  MISSING  $filepath"
			ERRORS=$((ERRORS + 1))
		fi
	fi
done <rdistnn.pkg

# ── 2. .ado / .sthlp pairing ────────────────────────────────────────────────
echo ""
echo "=== Checking .ado / .sthlp pairing in src/ ==="
for ado in src/*.ado; do
	base=$(basename "$ado" .ado)
	if [[ "$base" == _* ]]; then
		echo "  SKIP     $ado  (internal helper)"
		continue
	fi
	sthlp="src/${base}.sthlp"
	if [[ -f "$sthlp" ]]; then
		echo "  OK       $base  (.ado + .sthlp)"
	else
		echo "  MISSING  $sthlp"
		ERRORS=$((ERRORS + 1))
	fi
done

# ── 3. public command manifest entries ──────────────────────────────────────
echo ""
echo "=== Checking public command manifest entries ==="
for ado in src/*.ado; do
	base=$(basename "$ado" .ado)
	if [[ "$base" == _* ]]; then
		echo "  SKIP     $ado  (internal helper)"
		continue
	fi
	sthlp="src/${base}.sthlp"
	if grep -qE "^f[[:space:]]+${ado}([[:space:]]|$)" rdistnn.pkg; then
		echo "  OK       rdistnn.pkg includes $ado"
	else
		echo "  MISSING  rdistnn.pkg has no 'f $ado' entry"
		ERRORS=$((ERRORS + 1))
	fi
	if grep -qE "^f[[:space:]]+${sthlp}([[:space:]]|$)" rdistnn.pkg; then
		echo "  OK       rdistnn.pkg includes $sthlp"
	else
		echo "  MISSING  rdistnn.pkg has no 'f $sthlp' entry"
		ERRORS=$((ERRORS + 1))
	fi
done

# ── 4. plugin manifest entry ────────────────────────────────────────────────
echo ""
echo "=== Checking plugin manifest entry ==="
if grep -qE "^f[[:space:]]+src/rdistnn[.]plugin([[:space:]]|$)" rdistnn.pkg; then
	echo "  OK       rdistnn.pkg includes src/rdistnn.plugin"
else
	echo "  MISSING  rdistnn.pkg has no 'f src/rdistnn.plugin' entry"
	ERRORS=$((ERRORS + 1))
fi

# ── 5. Stata version floor ──────────────────────────────────────────────────
echo ""
echo "=== Checking Stata version floor ==="
VERSION_ERRORS=0
for f in src/*.ado tests/*.do; do
	if grep -nE "^[[:space:]]*version[[:space:]]+16([[:space:]]|$)" "$f" >/dev/null; then
		echo "  OLD      $f declares version 16; use version 17 or later"
		grep -nE "^[[:space:]]*version[[:space:]]+16([[:space:]]|$)" "$f" \
			| sed "s/^/           /"
		VERSION_ERRORS=$((VERSION_ERRORS + 1))
	fi
done

if [[ $VERSION_ERRORS -eq 0 ]]; then
	echo "  OK       active ado/test files declare Stata 17+"
else
	ERRORS=$((ERRORS + VERSION_ERRORS))
fi

# ── 6. stata.toc has a package entry for rdistnn ────────────────────────────
# In stata.toc format, packages are listed as "p <name>" (no .pkg extension).
echo ""
echo "=== Checking stata.toc ==="
if [[ -f stata.toc ]] && grep -qE "^p[[:space:]]+rdistnn([[:space:]]|$)" stata.toc; then
	echo "  OK       stata.toc has entry: p rdistnn"
else
	echo "  MISSING  stata.toc has no 'p rdistnn' package entry"
	ERRORS=$((ERRORS + 1))
fi

# ── Result ───────────────────────────────────────────────────────────────────
echo ""
if [[ $ERRORS -gt 0 ]]; then
	echo "pkg-check FAILED with $ERRORS error(s)."
	exit 1
else
	echo "pkg-check passed."
fi
