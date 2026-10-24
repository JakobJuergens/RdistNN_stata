#!/usr/bin/env python3
"""Inspect actual Stata SVG paths against exported numeric plot-data CSVs.

Standard-library only; does not run Stata or mutate source/exported evidence.
Run after tests/test-reporting-export-contracts.do; writes an ignored JSON receipt
in --evidence-root (default: build/premerge).
"""

import argparse
import csv
import json
import re
from pathlib import Path
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[2] / "build/premerge"
NS = "{http://www.w3.org/1999/xlink}href"
NUM = r"[-+]?(?:\d*\.\d+|\d+\.?)(?:[eE][-+]?\d+)?"


def points(path):
    d = path.attrib.get("d", "")
    if re.search(r"[ACQSTHV]", d, re.I):
        return []
    values = [float(v) for v in re.findall(NUM, d)]
    if len(values) % 2:
        return []
    answer = []
    for point in zip(values[::2], values[1::2]):
        if not answer or point != answer[-1]:
            answer.append(point)
    if len(answer) > 1 and answer[-1] == answer[0]:
        answer.pop()
    return answer


def affine(xs, ys):
    xbar, ybar = sum(xs) / len(xs), sum(ys) / len(ys)
    denominator = sum((x - xbar) ** 2 for x in xs)
    assert denominator > 0, "degenerate independent affine axis"
    slope = sum((x - xbar) * (y - ybar) for x, y in zip(xs, ys)) / denominator
    return ybar - slope * xbar, slope


def area(poly):
    return abs(sum(a[0] * b[1] - a[1] * b[0]
                   for a, b in zip(poly, poly[1:] + poly[:1]))) / 2


def svg_read(filename):
    tree = ET.parse(ROOT / filename).getroot()
    paths = [x for x in tree.iter() if x.tag.endswith("}path")]
    return tree, paths


def filled(tree, path):
    ident = path.attrib.get("id")
    styles = [path.attrib.get("style", "")]
    if ident:
        styles += [x.attrib.get("style", "") for x in tree.iter()
                   if x.attrib.get(NS) == "#" + ident]
    return any(re.search(r"(?:^|;)fill:(?!none)[^;]+", style)
               and not re.search(r"(?:^|;)fill-opacity:0(?:;|$)", style)
               for style in styles)


def check_band(stem, outcome=False):
    tree, paths = svg_read(stem + ".svg")
    with (ROOT / (stem + "-data.csv")).open(newline="") as stream:
        data = sorted(csv.DictReader(stream), key=lambda r: float(r["x0_1"]))
    x = [float(r["x0_1"]) for r in data]
    mean = "weighted_mean_y" if outcome else "estimate"
    lower = "weighted_y_q10" if outcome else "ci_lower"
    upper = "weighted_y_q90" if outcome else "ci_upper"
    y = [float(r[mean]) for r in data]
    lines = [points(p) for p in paths if "Z" not in p.attrib.get("d", "").upper()
             and len(points(p)) == len(data)]
    assert len(lines) == 1, (stem, "expected one estimate line", len(lines))
    line = lines[0]
    xa, xb = affine(x, [p[0] for p in line])
    ya, yb = affine(y, [p[1] for p in line])
    assert xb > 0 and yb < 0, (stem, "axis direction")
    residual = max(max(abs(px - (xa + xb * dx)), abs(py - (ya + yb * dy)))
                   for dx, dy, (px, py) in zip(x, y, line))
    assert residual < .5, (stem, "mean-line geometric error", residual)
    bands = [p for p in paths if "Z" in p.attrib.get("d", "").upper()
             and len(points(p)) == 2 * len(data) and filled(tree, p)]
    assert len(bands) == 1, (stem, "expected one genuine filled band", len(bands))
    polygon = points(bands[0])
    expected = [(xa + xb * float(row["x0_1"]), ya + yb * float(row[bound]))
                for bound in (lower, upper) for row in data]
    maximum_error = max(min(max(abs(a - c), abs(b - d)) for c, d in polygon)
                        for a, b in expected)
    # Stata rounds the different SVG layers separately to hundredths; the
    # ~0.17-unit layer displacement is far below one rendered pixel here.
    assert maximum_error < .5, (stem, "bounds do not match plot data", maximum_error)
    polygon_area = area(polygon)
    assert polygon_area > 1, (stem, "collapsed filled band", polygon_area)
    return {"file": stem + ".svg", "data_rows": len(data),
            "x_range": [min(x), max(x)], "filled_band_vertices": len(polygon),
            "filled_band_area_svg_units2": polygon_area,
            "mean_line_max_error_svg_units": residual,
            "bound_vertex_max_error_svg_units": maximum_error,
            "axes_and_titles": ["".join(t.itertext()) for t in tree.iter()
                                if t.tag.endswith("}text")]}


def check_nose(family):
    filename = f"postestimation-{family}-nose-profile.svg"
    tree, paths = svg_read(filename)
    bands = [p for p in paths if "Z" in p.attrib.get("d", "").upper()
             and area(points(p)) > 1 and filled(tree, p)]
    assert not bands, (filename, "nose unexpectedly has a filled confidence band")
    lines = [points(p) for p in paths if len(points(p)) == 4
             and "Z" not in p.attrib.get("d", "").upper()]
    assert len(lines) == 1 and len({p[0] for p in lines[0]}) == 4, (filename, "missing four-point estimate profile")
    text = ["".join(t.itertext()) for t in tree.iter() if t.tag.endswith("}text")]
    assert "x0_1" in text, (filename, "coordinate-collision axis label")
    return {"file": filename, "filled_band_count": 0, "estimate_vertices": 4,
            "axes_and_titles": text}


def main():
    global ROOT
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--evidence-root", type=Path, default=ROOT)
    ROOT = parser.parse_args().evidence_root.resolve()
    receipt = {"description": "Actual SVG path geometry compared with exported numeric data; no Stata run by this checker.",
               "positive_bands": [], "nose_profiles": []}
    for family in ("cate", "np", "clate"):
        receipt["positive_bands"].append(check_band(f"postestimation-{family}-interval-profile"))
    receipt["positive_bands"].append(check_band("postestimation-np-outcome-profile", outcome=True))
    for family in ("cate", "np", "clate"):
        receipt["nose_profiles"].append(check_nose(family))
    output = ROOT / "postestimation-svg-receipt.json"
    output.write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps({"positive_bands": len(receipt["positive_bands"]),
                      "nose_profiles": len(receipt["nose_profiles"]),
                      "receipt": str(output)}, indent=2))


if __name__ == "__main__":
    main()
