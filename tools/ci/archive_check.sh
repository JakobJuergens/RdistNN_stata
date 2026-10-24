#!/usr/bin/env bash
# Validate the source ZIP's own installation inventory and public payload.
# Usage: bash tools/ci/archive_check.sh /tmp/rdistnn-source.zip
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

# Python3 is already required by the source builder and package hygiene gate.
# Read the archive directly: no extraction or dependence on checkout contents.
python3 - "$ZIP_PATH" "$SCRIPT_DIR" <<'PY'
from collections import Counter
from pathlib import PurePosixPath
import stat
import sys
import zipfile

sys.path.insert(0, sys.argv[2])
from package_hygiene import installation_manifest, relative_path_issue

errors = []
with zipfile.ZipFile(sys.argv[1]) as archive:
    entries = archive.infolist()
    names = [entry.filename for entry in entries]
    for name, count in Counter(names).items():
        if count > 1:
            errors.append('duplicate archive entry: ' + name)
    files = {entry.filename for entry in entries if not entry.is_dir()}
    for required in ('README.md', 'rdistnn.pkg', 'stata.toc'):
        if required not in files:
            errors.append('missing archive input: ' + required)
    manifest_files = []
    if 'rdistnn.pkg' in files:
        installed = {name for name in files
                     if PurePosixPath(name).parent == PurePosixPath('src')
                     and PurePosixPath(name).suffix in ('.ado', '.sthlp')}
        manifest_files, manifest_errors = installation_manifest(
            archive.read('rdistnn.pkg').decode('utf-8'), files.__contains__, installed,
        )
        errors.extend(manifest_errors)
    # The builder includes public source/docs/examples plus root installation
    # inputs. Other checkout-only roots and generated output are not payload.
    root_files = {'README.md', 'rdistnn.pkg', 'stata.toc'} | {
        name for name in manifest_files if len(PurePosixPath(name).parts) == 1
    }
    generated_suffixes = {'.log', '.o', '.obj', '.so', '.dll', '.dylib', '.mlib', '.mo', '.pyc', '.pyo'}
    for entry in entries:
        name = entry.filename.rstrip('/') if entry.is_dir() else entry.filename
        issue = relative_path_issue(name)
        if issue:
            errors.append(f'archive entry {issue}: {entry.filename}')
            continue
        path = PurePosixPath(name)
        if path.parts[0] not in ('src', 'docs', 'examples') and name not in root_files:
            errors.append('nonpublic archive input: ' + entry.filename)
        if (any(part.startswith('.') or part == '__pycache__' for part in path.parts)
                or path.suffix.lower() in generated_suffixes
                or name == 'src/plugin/rdistnn.plugin'):
            errors.append('generated/private archive input: ' + entry.filename)
        kind = stat.S_IFMT(entry.external_attr >> 16)
        if kind not in (0, stat.S_IFREG, stat.S_IFDIR):
            errors.append('archive input must be an ordinary file/directory: ' + entry.filename)
    corrupt = archive.testzip()
    if corrupt:
        errors.append('corrupt archive input: ' + corrupt)
if errors:
    raise SystemExit('archive-check FAILED:\n' + '\n'.join('- ' + error for error in errors))
print(f'archive-check passed: {len(files)} public files, {len(manifest_files)} installation inputs.')
PY
