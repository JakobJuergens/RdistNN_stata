#!/usr/bin/env python3
"""Check installation completeness, private Mata names, and optional R mirrors."""
import argparse
import re
from collections import Counter
from pathlib import Path, PurePosixPath

SHARED_HEADERS = ('compensated_arithmetic.hpp', 'rank_weight_tables.hpp')
MATA_FUNCTION = re.compile(
    r'^[ \t]*(?:(?:real|string|complex|transmorphic)\s+(?:scalar|vector|rowvector|colvector|matrix)|void)'
    r'\s+([A-Za-z_]\w*)\s*\(', re.M,
)


def package_issues(root):
    manifest = root / 'rdistnn.pkg'
    if not manifest.is_file():
        return ['missing rdistnn.pkg']
    files = re.findall(r'^[fF]\s+(\S+)', manifest.read_text(encoding='utf-8'), re.M)
    errors = []
    for name, count in Counter(str(PurePosixPath(x)) for x in files).items():
        if count > 1:
            errors.append(f'duplicate installation entry: {name}')
    for name in files:
        path = PurePosixPath(name)
        if path.is_absolute() or '..' in path.parts or '\\' in name or ':' in name:
            errors.append(f'installation entry must be a local relative path: {name}')
        elif str(path) != name:
            errors.append(f'installation entry must use a canonical relative path: {name}')
        elif not (root / name).is_file():
            errors.append(f'missing installation input: {name}')
    required = {p.relative_to(root).as_posix() for pattern in ('*.ado', '*.sthlp')
                for p in (root / 'src').glob(pattern)}
    for name in sorted(required.difference(files)):
        errors.append(f'unregistered installation input: {name}')
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
    if args.r_root is not None:
        errors.extend(shared_header_issues(args.repo_root, args.r_root))
    if errors:
        parser.exit(1, 'Package hygiene failed:\n' + '\n'.join(f'- {x}' for x in errors) + '\n')
    print('Installation inventory and private Mata namespaces pass.' +
          (' Shared R/C++ headers match.' if args.r_root is not None else ''))


if __name__ == '__main__':
    main()
