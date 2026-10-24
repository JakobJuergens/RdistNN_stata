#!/usr/bin/env python3
"""Check installation completeness, private Mata names, and optional R mirrors."""
import argparse
import re
from collections import Counter
from pathlib import Path, PurePosixPath

SHARED_HEADERS = ('compensated_arithmetic.hpp', 'distance_arithmetic.hpp')
PACKAGE_VERSION = r'[0-9]+\.[0-9]+\.[0-9]+(?:\.[0-9]+)?'
MATA_FUNCTION = re.compile(
    r'^[ \t]*(?:(?:real|string|complex|transmorphic)\s+(?:scalar|vector|rowvector|colvector|matrix)|void)'
    r'\s+([A-Za-z_]\w*)\s*\(', re.M,
)


def relative_path_issue(name):
    """Return the shared path restriction violated by an installation input."""
    path = PurePosixPath(name)
    if path.is_absolute() or '..' in path.parts or '\\' in name or ':' in name:
        return 'must be a local relative path'
    if str(path) != name:
        return 'must use a canonical relative path'
    return None


def installation_manifest(text, file_exists, required=()):
    """Check a manifest against filesystem or archive inventory without I/O."""
    files = re.findall(r'^[fF]\s+(\S+)', text, re.M)
    errors = []
    for name, count in Counter(str(PurePosixPath(x)) for x in files).items():
        if count > 1:
            errors.append(f'duplicate installation entry: {name}')
    for name in files:
        issue = relative_path_issue(name)
        if issue:
            errors.append(f'installation entry {issue}: {name}')
        elif not file_exists(name):
            errors.append(f'missing installation input: {name}')
    for name in sorted(set(required).difference(files)):
        errors.append(f'unregistered installation input: {name}')
    return files, errors


def package_issues(root):
    manifest = root / 'rdistnn.pkg'
    if not manifest.is_file():
        return ['missing rdistnn.pkg']
    required = {p.relative_to(root).as_posix() for pattern in ('*.ado', '*.sthlp')
                for p in (root / 'src').glob(pattern)}
    _, errors = installation_manifest(
        manifest.read_text(encoding='utf-8'), lambda name: (root / name).is_file(), required,
    )
    function_owners = {}
    for source in sorted((root / 'src').glob('*.ado')):
        text = source.read_text(encoding='utf-8')
        versions = re.findall(r'^\s*version\s+([0-9]+(?:[.][0-9]+)?)\b', text, re.M)
        if not versions or any(float(x) < 17 for x in versions):
            errors.append(f'{source.name}: declare version 17 or later')
        for name in MATA_FUNCTION.findall(text):
            if not name.startswith('rdistnn_'):
                errors.append(f'{source.name}: private Mata function needs rdistnn_ prefix: {name}')
            if name in function_owners:
                errors.append(f'duplicate Mata function {name}: {function_owners[name]}, {source.name}')
            else:
                function_owners[name] = source.name
    return errors


def package_version_issues(root):
    """Check current package claims, preserving file revisions and release pins."""
    manifest = root / 'rdistnn.pkg'
    if not manifest.is_file():
        return ['missing rdistnn.pkg package version']
    declarations = re.findall(r'^d Version:([^\n]*)$',
                              manifest.read_text(encoding='utf-8'), re.M)
    if len(declarations) != 1 or not re.fullmatch(PACKAGE_VERSION, declarations[0].strip()):
        return ['rdistnn.pkg: declare exactly one three- or four-component package Version']
    version = declarations[0].strip()
    owners = (
        ('README.md', r'^Version \*\*(' + PACKAGE_VERSION + r')\*\*'),
        ('RELEASE_NOTES.md', r'^##\s+(' + PACKAGE_VERSION + r')(?=\s|$)'),
        ('src/rdistnn.sthlp', r'\bPackage \{bf:(' + PACKAGE_VERSION + r')\}'),
    )
    errors = []
    for name, pattern in owners:
        path = root / name
        if not path.is_file():
            errors.append(f'{name}: missing current package version owner')
            continue
        text = path.read_text(encoding='utf-8')
        matches = re.findall(pattern, text, re.M)
        # Later release-note headings record immutable historical versions.
        if name == 'RELEASE_NOTES.md':
            headings = re.findall(r'^##\s+([^\n]+)', text, re.M)
            first = re.match(r'(' + PACKAGE_VERSION + r')(?=\s|$)', headings[0]) if headings else None
            matches = [first[1]] if first else []
        if len(matches) != 1:
            errors.append(f'{name}: declare one current package version ({version})')
        elif matches[0] != version:
            errors.append(f'{name}: current package version {matches[0]} differs from rdistnn.pkg {version}')

    claim = re.compile(r'\bpackage(?:\s+version)?\s+(?:\*\*|\{bf:)('
                       + PACKAGE_VERSION + r')(?:\*\*|\})', re.I)
    development = re.compile(r'\bcurrent\s+\*\*(' + PACKAGE_VERSION + r') development\*\*')
    # Public articles/help describe the current interface unless a literal is
    # explicitly qualified as production or historical. Headers, language floors,
    # native/schema versions, fixed test provenance and release URLs are separate.
    for path in sorted([*(root / 'src').glob('*.sthlp'), *(root / 'docs').glob('*.md')]):
        name = path.relative_to(root).as_posix()
        if name == 'src/rdistnn.sthlp':
            continue
        for line in path.read_text(encoding='utf-8').splitlines():
            if re.search(r'\b(production|historical|released|previous|older|pinned)\b', line, re.I):
                continue
            for match in [*claim.finditer(line), *development.finditer(line)]:
                if match[1] != version:
                    errors.append(f'{name}: current package version {match[1]} differs from rdistnn.pkg {version}')
    return errors


def shared_header_issues(root, r_root):
    errors = []
    for name in SHARED_HEADERS:
        stata = root / 'src/plugin' / name
        r = r_root / 'inst/include' / name
        if not stata.is_file() or not r.is_file():
            errors.append(f'missing mirrored header: {name}')
        elif stata.read_bytes() != r.read_bytes():
            errors.append(f'mirrored header differs: {name}')
    return errors


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo-root', type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument('--r-root', type=Path, help='also require identical shared headers in this R checkout')
    args = parser.parse_args()
    errors = package_issues(args.repo_root)
    errors.extend(package_version_issues(args.repo_root))
    if args.r_root is not None:
        errors.extend(shared_header_issues(args.repo_root, args.r_root))
    if errors:
        parser.exit(1, 'Package hygiene failed:\n' + '\n'.join(f'- {x}' for x in errors) + '\n')
    print('Package version, installation inventory and private Mata namespaces pass.' +
          (' Shared R/C++ headers match.' if args.r_root is not None else ''))


if __name__ == '__main__':
    main()
