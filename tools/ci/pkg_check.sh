#!/usr/bin/env bash
# pkg_check.sh — validate RdistNN Stata package integrity
#
# Checks:
#   1. Every file listed under "f " in rdistnn.pkg exists on disk.
#   2. Every public .ado command in src/ has a matching .sthlp help file.
#   3. Every public .ado/.sthlp surface is listed in rdistnn.pkg.
#   4. The installable plugin is listed in rdistnn.pkg.
#   5. Active ado/test entry files target Stata 17+.
#   6. cate_dnn help retains the nuisance-inference contract.
#   7. DNN kernel-order documentation retains the exact support direction.
#   8. Experimental PLQR and CATE paper-diagnostic contracts remain explicit.
#   9. Generated archives and CI logs are not tracked.
#  10. stata.toc references rdistnn.pkg.
#
# Exit code: 0 = all checks passed, 1 = one or more failures.
# Run from the repository root.

set -euo pipefail

ERRORS=0

# ── 1. rdistnn.pkg file manifest ────────────────────────────────────────────
echo "=== Checking rdistnn.pkg file entries ==="
while IFS= read -r line; do
	if [[ "$line" =~ ^f[[:space:]]+(.+)$ ]]; then
		filepath="${BASH_REMATCH[1]%$'\r'}"
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

# ── 6. cate_dnn nuisance-inference documentation contract ──────────────────
echo ""
echo "=== Checking cate_dnn inference documentation ==="
DOC_CONTRACT_ERRORS=0
for phrase in \
	"{title:Inference contract}" \
	"Orthogonality controls first-order population nuisance drift" \
	"same-sample nuisance estimates" \
	"Valid full-sample reuse additionally requires" \
	"nuisance-rate" \
	"full-state stability" \
	"cannot diagnose"; do
	if grep -qF "$phrase" src/cate_dnn.sthlp; then
		echo "  OK       cate_dnn.sthlp contains: $phrase"
	else
		echo "  MISSING  cate_dnn.sthlp must contain: $phrase"
		DOC_CONTRACT_ERRORS=$((DOC_CONTRACT_ERRORS + 1))
	fi
done

if [[ $DOC_CONTRACT_ERRORS -gt 0 ]]; then
	ERRORS=$((ERRORS + DOC_CONTRACT_ERRORS))
fi

# ── 7. exact DNN support-direction documentation ───────────────────────────
echo ""
echo "=== Checking DNN kernel-order documentation ==="
DNN_DOC_ERRORS=0
DNN_DOCS=(
	README.md
	src/cate_dnn.sthlp
	src/np_dnn.sthlp
	src/plqr_dnn.sthlp
	src/weights_dnn.sthlp
)
for file in "${DNN_DOCS[@]}"; do
	for phrase in "n - s + k" "n - s + 1" "opposite to the usual bandwidth convention"; do
		if grep -qFi "$phrase" "$file"; then
			echo "  OK       $file contains: $phrase"
		else
			echo "  MISSING  $file must contain: $phrase"
			DNN_DOC_ERRORS=$((DNN_DOC_ERRORS + 1))
		fi
	done
done

DNN_DOC_TEXT=""
for file in "${DNN_DOCS[@]}"; do
	DNN_DOC_TEXT+=" $(tr -d '\r`' <"$file" | tr '\n' ' ')"
done
if printf '%s\n' "$DNN_DOC_TEXT" | grep -Eqi \
	'larger[[:space:]]+s[[:space:]]+(includes|spreads|assigns|puts)[^.]{0,80}more'; then
	echo "  INVALID  documentation says larger s spreads weight over more observations"
	DNN_DOC_ERRORS=$((DNN_DOC_ERRORS + 1))
fi
if printf '%s\n' "$DNN_DOC_TEXT" | grep -Eqi \
	'larger[[:space:]]+s[[:space:]]+(increases|expands|widens)[^.]{0,40}support'; then
	echo "  INVALID  documentation says larger s expands positive-weight support"
	DNN_DOC_ERRORS=$((DNN_DOC_ERRORS + 1))
fi

if [[ $DNN_DOC_ERRORS -eq 0 ]]; then
	echo "  OK       larger s is documented as tighter localization"
else
	ERRORS=$((ERRORS + DNN_DOC_ERRORS))
fi

# ── 8. experimental PLQR and CATE paper-diagnostic contracts ───────────────
echo ""
echo "=== Checking estimator-surface documentation contracts ==="
SURFACE_ERRORS=0
for phrase in \
	"Experimental status" \
	"not a coefficient estimator" \
	"No local root solving was performed" \
	"thetaeval()" \
	"Compatibility aliases"; do
	if grep -qFi "$phrase" src/plqr_dnn.sthlp src/plqr_dnn.ado; then
		echo "  OK       PLQR surface contains: $phrase"
	else
		echo "  MISSING  PLQR surface must contain: $phrase"
		SURFACE_ERRORS=$((SURFACE_ERRORS + 1))
	fi
done
if grep -q "ereturn display" src/plqr_dnn.ado; then
	echo "  INVALID  experimental PLQR display must not use coefficient-style ereturn display"
	SURFACE_ERRORS=$((SURFACE_ERRORS + 1))
fi
for phrase in physical_radius_90 treated_ess control_ess top5_score_share; do
	if grep -qF "$phrase" src/cate_dnn.ado src/cate_dnn_diagnostics.ado; then
		echo "  OK       CATE diagnostics contain: $phrase"
	else
		echo "  MISSING  CATE diagnostics must contain: $phrase"
		SURFACE_ERRORS=$((SURFACE_ERRORS + 1))
	fi
done
ERRORS=$((ERRORS + SURFACE_ERRORS))

# ── 9. generated-artifact tracking hygiene ─────────────────────────────────
echo ""
echo "=== Checking generated-artifact tracking hygiene ==="
TRACKED_GENERATED="$(git ls-files -- GH_logs archives/source.zip)"
if [[ -n "$TRACKED_GENERATED" ]]; then
	echo "  TRACKED  generated files must be removed:"
	printf '%s\n' "$TRACKED_GENERATED" | sed 's/^/           /'
	ERRORS=$((ERRORS + 1))
else
	echo "  OK       GH_logs and archives/source.zip are not tracked"
fi
for pattern in "GH_logs/" "archives/source.zip"; do
	if grep -qFx "$pattern" .gitignore; then
		echo "  OK       .gitignore contains: $pattern"
	else
		echo "  MISSING  .gitignore must contain: $pattern"
		ERRORS=$((ERRORS + 1))
	fi
done

# ── 10. stata.toc has a package entry for rdistnn ───────────────────────────
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
