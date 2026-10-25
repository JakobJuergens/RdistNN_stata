#!/usr/bin/env python3
"""Check installation completeness, private Mata names, and optional R mirrors."""
import argparse
import datetime
import re
from collections import Counter
from pathlib import Path, PurePosixPath

SHARED_HEADERS = ('compensated_arithmetic.hpp', 'distance_arithmetic.hpp')
PACKAGE_VERSION = r'[0-9]+\.[0-9]+\.[0-9]+(?:\.[0-9]+)?'
MATA_FUNCTION = re.compile(
    r'^[ \t]*(?:(?:real|string|complex|transmorphic)\s+(?:scalar|vector|rowvector|colvector|matrix)|void)'
    r'\s+([A-Za-z_]\w*)\s*\(', re.M,
)
STAMP_MARKER = re.compile(r'^(?:\{\* )?\*![ \t]*version\b')
STAMP = re.compile(r'\*! version (' + PACKAGE_VERSION
                   + r')[ \t]+([0-9]{4}-[0-9]{2}-[0-9]{2})(?:[ \t]+.*)?')


def package_version(text):
    """Read the single package version owner from a checkout or archive."""
    declarations = re.findall(r'^d Version:([^\n]*)$', text, re.M)
    if len(declarations) != 1 or not re.fullmatch(PACKAGE_VERSION, declarations[0].strip()):
        raise ValueError('rdistnn.pkg: declare exactly one three- or four-component package Version')
    return declarations[0].strip()


def source_version_stamp_issues(name, text, version):
    """Check one installed ado/help stamp without interpreting language/schema versions."""
    lines = text.splitlines()
    stamps = [(i, line) for i, line in enumerate(lines) if STAMP_MARKER.match(line)]
    if len(stamps) != 1:
        return [f'{name}: declare exactly one package version stamp ({version})']
    index, line = stamps[0]
    help_file = name.endswith('.sthlp')
    if index != int(help_file) or (help_file and lines[0] != '{smcl}'):
        return [f'{name}: place the package version stamp at the top, after SMCL for help']
    if help_file:
        if not line.startswith('{* ') or not line.endswith('}{...}'):
            return [f'{name}: malformed SMCL package version stamp']
        line = line[3:-len('}{...}')]
    match = STAMP.fullmatch(line)
    if match is None:
        return [f'{name}: package version stamp must contain a version and ISO date']
    errors = []
    if match[1] != version:
        errors.append(f'{name}: package version stamp {match[1]} differs from rdistnn.pkg {version}')
    try:
        datetime.date.fromisoformat(match[2])
    except ValueError:
        errors.append(f'{name}: package version stamp must contain a valid ISO date')
    return errors


def source_version_paths(root):
    return sorted([*(root / 'src').glob('*.ado'), *(root / 'src').glob('*.sthlp')])


def error_rendering_issues(root):
    """Enforce custom rendering boundaries while permitting native validation/signaling."""
    errors = []
    ado_error = re.compile(r'\b(?:di|dis|disp|displ|displa|display)\s+as\s+err(?:o(?:r)?)?\b')
    for path in sorted((root / 'src').glob('*.ado')):
        lines = path.read_text(encoding='utf-8').splitlines()
        for index, line in enumerate(lines):
            if line.lstrip().startswith(('*', '//')):
                continue
            # Ignore literal examples in owned messages; this is a source-policy
            # check, not a replacement for the Stata interpreter/runtime tests.
            code = re.sub(r"""`".*?"'""", '""', line)
            code = re.sub(r'"(?:[^"\\]|\\.)*"', '""', code).split('//', 1)[0]
            if re.search(r'\berrprintf\s*\(', code):
                errors.append(f'{path.name}:{index + 1}: render custom Mata errors through rdistnn_fail()')
            if not ado_error.search(code):
                continue
            renderer = path.name == '_rdistnn_messages.ado' and 'display as error _asis' in code
            bootstrap = (path.name == '_rdistnn_mata.ado' and 'display as error _asis' in code
                         and "`plugin_rc'" in line and index + 1 < len(lines)
                         and lines[index + 1].strip() == "exit `plugin_rc'")
            if not (renderer or bootstrap):
                errors.append(f'{path.name}:{index + 1}: render custom ado errors through _rdistnn_messages')
    return errors


def sync_source_version_stamps(root, version, date):
    """Prepare all header edits before writing; preserve executable text and other notes."""
    changes = []
    for path in source_version_paths(root):
        name = path.relative_to(root).as_posix()
        original = path.read_bytes()
        lines = original.decode('utf-8').splitlines(keepends=True)
        help_file = path.suffix == '.sthlp'
        index = int(help_file)
        stamps = [(i, line) for i, line in enumerate(lines) if STAMP_MARKER.match(line)]
        if len(stamps) > 1 or (stamps and stamps[0][0] != index):
            raise ValueError(f'{name}: ambiguous package version stamps; resolve them before syncing')
        if help_file and (not lines or lines[0].rstrip('\r\n') != '{smcl}'):
            raise ValueError(f'{name}: help must begin with SMCL before syncing its version stamp')
        suffix = ''
        if stamps:
            old = stamps[0][1].rstrip('\r\n')
            if help_file:
                if not old.startswith('{* ') or not old.endswith('}{...}'):
                    raise ValueError(f'{name}: malformed SMCL version stamp')
                old = old[3:-len('}{...}')]
            parts = old.split(maxsplit=4)
            if len(parts) < 4 or parts[:2] != ['*!', 'version']:
                raise ValueError(f'{name}: malformed version stamp')
            if len(parts) == 5:
                suffix = '  ' + parts[4]
        newline = '\r\n' if lines and lines[0].endswith('\r\n') else '\n'
        stamp = f'*! version {version}  {date.isoformat()}{suffix}'
        if help_file:
            stamp = '{* ' + stamp + '}{...}'
        if stamps:
            lines[index] = stamp + newline
        else:
            lines.insert(index, stamp + newline)
        updated = ''.join(lines).encode('utf-8')
        issues = source_version_stamp_issues(name, updated.decode('utf-8'), version)
        if issues:
            raise ValueError('\n'.join(issues))
        if updated != original:
            changes.append((path, updated))
    for path, updated in changes:
        path.write_bytes(updated)
    return len(changes)


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
    """Check package claims and source stamps while preserving historical release pins."""
    manifest = root / 'rdistnn.pkg'
    if not manifest.is_file():
        return ['missing rdistnn.pkg package version']
    try:
        version = package_version(manifest.read_text(encoding='utf-8'))
    except ValueError as error:
        return [str(error)]
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
    # explicitly qualified as production or historical. Language floors,
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
    for path in source_version_paths(root):
        errors.extend(source_version_stamp_issues(
            path.relative_to(root).as_posix(), path.read_text(encoding='utf-8'), version,
        ))
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
    parser.add_argument('--sync-version-stamps', action='store_true',
                        help='synchronize ado/help comment stamps to rdistnn.pkg before checking')
    parser.add_argument('--stamp-date', type=datetime.date.fromisoformat,
                        help='ISO date for synchronized stamps (default: today)')
    args = parser.parse_args()
    if args.stamp_date and not args.sync_version_stamps:
        parser.error('--stamp-date requires --sync-version-stamps')
    if args.sync_version_stamps:
        try:
            version = package_version((args.repo_root / 'rdistnn.pkg').read_text(encoding='utf-8'))
            changed = sync_source_version_stamps(
                args.repo_root, version, args.stamp_date or datetime.date.today(),
            )
        except ValueError as error:
            parser.exit(1, str(error) + '\n')
        print(f'Synchronized {changed} ado/help version stamps to {version}.')
    errors = package_issues(args.repo_root)
    errors.extend(package_version_issues(args.repo_root))
    errors.extend(error_rendering_issues(args.repo_root))
    if args.r_root is not None:
        errors.extend(shared_header_issues(args.repo_root, args.r_root))
    if errors:
        parser.exit(1, 'Package hygiene failed:\n' + '\n'.join(f'- {x}' for x in errors) + '\n')
    print('Package version, installation inventory and private Mata namespaces pass.' +
          (' Shared R/C++ headers match.' if args.r_root is not None else ''))


if __name__ == '__main__':
    main()
