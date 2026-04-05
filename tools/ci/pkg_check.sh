#!/usr/bin/env bash
# pkg_check.sh — validate RdistNN Stata package integrity
#
# Checks:
#   1. Every file listed under "f " in rdistnn.pkg exists on disk.
#   2. Every public .ado command in src/ has a matching .sthlp help file.
#   3. stata.toc references rdistnn.pkg.
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
done < rdistnn.pkg

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

# ── 3. stata.toc has a package entry for rdistnn ────────────────────────────
# In stata.toc format, packages are listed as "p <name>" (no .pkg extension).
echo ""
echo "=== Checking stata.toc ==="
if grep -qE "^p[[:space:]]+rdistnn([[:space:]]|$)" stata.toc; then
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
