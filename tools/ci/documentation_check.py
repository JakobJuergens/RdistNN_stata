#!/usr/bin/env python3
"""Check public documentation links, installed help and optional R topic coverage."""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path
from urllib.parse import unquote, urlsplit


def headings(path: Path) -> set[str]:
    anchors: set[str] = set()
    counts: dict[str, int] = {}
    fenced = False
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.lstrip().startswith("```"):
            fenced = not fenced
        if fenced or not re.match(r"^#{1,6}\s", line):
            continue
        title = re.sub(r"^#+\s+|\s+#+\s*$", "", line).strip().lower()
        title = re.sub(r"[^\w\-\s]", "", title)
        slug = re.sub(r"\s", "-", title)
        count = counts.get(slug, 0)
        counts[slug] = count + 1
        anchors.add(slug if count == 0 else f"{slug}-{count}")
    return anchors


def check(root: Path, r_root: Path | None) -> list[str]:
    errors: list[str] = []
    markdown = [p for p in root.glob("*.md") if p.name not in {"AGENTS.md", "CLAUDE.md"}]
    for folder in ("docs", "examples", "tests", "src/plugin"):
        markdown.extend((root / folder).rglob("*.md"))
    links = 0
    for path in markdown:
        content = path.read_text(encoding="utf-8")
        for target in re.findall(r"\[[^\]]*\]\(([^\s)]+)(?:\s+\"[^\"]*\")?\)", content):
            parsed = urlsplit(target)
            if parsed.scheme or parsed.netloc:
                continue
            resolved = (path.parent / unquote(parsed.path)).resolve() if parsed.path else path
            if not resolved.is_relative_to(root):
                errors.append(f"{path.relative_to(root)}: nonportable local link {target}")
            elif not resolved.exists():
                errors.append(f"{path.relative_to(root)}: missing link {target}")
            elif parsed.fragment and resolved.suffix == ".md" and unquote(parsed.fragment) not in headings(resolved):
                errors.append(f"{path.relative_to(root)}: missing anchor {target}")
            links += 1

    manifest = (root / "rdistnn.pkg").read_text(encoding="utf-8")
    installed = set(re.findall(r"^[fF]\s+(\S+)", manifest, re.M))
    commands = sorted(p.stem for p in (root / "src").glob("*.ado") if not p.stem.startswith("_"))
    index = (root / "docs/index.md").read_text(encoding="utf-8")
    overview_path = root / "src/rdistnn.sthlp"
    overview = overview_path.read_text(encoding="utf-8") if overview_path.exists() else ""
    for name in ["rdistnn", *commands]:
        help_path = root / f"src/{name}.sthlp"
        if not help_path.exists():
            errors.append(f"Missing installed help for {name}")
        if f"src/{name}.sthlp" not in installed:
            errors.append(f"Help not in package manifest: {name}")
        if name != "rdistnn" and (name not in index or f"{{help {name}}}" not in overview):
            errors.append(f"Public command not routed by index/overview: {name}")
    helps = sorted((root / "src").glob("*.sthlp"))
    for path in helps:
        content = path.read_text(encoding="utf-8")
        if content.count("{") != content.count("}"):
            errors.append(f"Unbalanced SMCL braces: {path.name}")
        for topic, anchor in re.findall(r'\{viewerjumpto\s+"[^"]*"\s+"([^"#]+)##([^"#]+)"\}', content):
            destination = root / f"src/{topic}.sthlp"
            if not destination.exists() or f"{{marker {anchor}}}" not in destination.read_text(encoding="utf-8"):
                errors.append(f"Broken help jump: {path.name} -> {topic}##{anchor}")
        for topic in re.findall(r"\{helpb?\s+([a-z_]+)", content):
            if (topic == "rdistnn" or "_dnn" in topic) and not (root / f"src/{topic}.sthlp").exists():
                errors.append(f"Broken local help link: {path.name} -> {topic}")

    examples = sorted((root / "examples").glob("*.do"))
    for path in examples:
        content = path.read_text(encoding="utf-8")
        if not re.search(r"(?m)^set seed \d+", content):
            errors.append(f"Example has no explicit seed: {path.name}")
    print(f"Checked {len(markdown)} Markdown files, {links} local links, {len(helps)} helps, {len(commands)} commands, {len(examples)} examples.")
    if r_root:
        crosswalk = (root / "docs/r-documentation-crosswalk.md").read_text(encoding="utf-8")
        articles = sorted((r_root / "vignettes").glob("*.Rmd"))
        groups = re.findall(r"^  - title: (.+)$", (r_root / "_pkgdown.yml").read_text(encoding="utf-8").split("\nreference:\n", 1)[1], re.M)
        for topic in [p.name for p in articles] + groups:
            if topic not in crosswalk:
                errors.append(f"R documentation topic absent from crosswalk: {topic}")
        print(f"Checked R coverage against {len(articles)} articles and {len(groups)} reference families at {r_root}.")
    return errors


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--r-root", type=Path, help="Optional live sibling RdistNN source for topic-inventory checks")
    args = parser.parse_args()
    errors = check(args.root.resolve(), args.r_root.resolve() if args.r_root else None)
    for error in errors:
        print(f"ERROR: {error}", file=sys.stderr)
    if errors:
        return 1
    print("Documentation checks passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
