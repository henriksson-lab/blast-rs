#!/usr/bin/env python3
"""Audit Rust function provenance comments and naming.

Auditable Rust functions should carry one of these comments directly above the
function definition:

  /// NCBI: Blast_HSPFree (blast_hits.c:130).
  /// blast-rs: Native helper; not a direct NCBI C port.

The script also understands the older prose comments already present in the
tree, such as `Port of NCBI `Blast_HSPFree``.  That keeps the audit useful
while comments are migrated to the shorter standard form.
"""

from __future__ import annotations

import argparse
import json
import re
from collections import Counter, defaultdict
from pathlib import Path


FN_RE = re.compile(
    r"^(?:pub\s+(?:\([^)]*\)\s+)?)?"
    r"(?:unsafe\s+)?(?:extern\s+\"C\"\s+)?"
    r"fn\s+([A-Za-z_][A-Za-z0-9_]*)\b"
)

STANDARD_NCBI_RE = re.compile(r"NCBI:\s*`?([A-Za-z_][A-Za-z0-9_]*)`?")
NATIVE_RE = re.compile(r"blast-rs:\s*(.*)", re.IGNORECASE)
NAMING_NOTE_RE = re.compile(r"naming:\s*(.*)", re.IGNORECASE)
LEGACY_REF_RES = [
    re.compile(r"Port(?:-[a-z]+)? of (?:the [^`]* )?(?:NCBI )?`([^`]+)`"),
    re.compile(r"1-1 (?:port|translation) of (?:NCBI )?`([^`]+)`"),
    re.compile(r"(?:Rust ownership|Rust-owned|Rust) equivalent of (?:NCBI )?`([^`]+)`"),
    re.compile(r"(?:counterpart|wrapper) (?:of|for) (?:NCBI )?`([^`]+)`"),
    re.compile(r"NCBI `([A-Za-z_][A-Za-z0-9_:]*)`"),
]


def load_function_names(path: Path | None) -> set[str]:
    if path is None or not path.exists():
        return set()
    with path.open("r", encoding="utf-8") as handle:
        report = json.load(handle)
    return {fn["name"] for fn in report.get("functions", [])}


def attached_comment(lines: list[str], fn_index: int) -> str:
    comments: list[str] = []
    idx = fn_index - 1
    while idx >= 0:
        stripped = lines[idx].strip()
        if stripped == "" or stripped.startswith("#["):
            idx -= 1
            continue
        if stripped.startswith(("///", "//", "//!")):
            comments.append(stripped)
            idx -= 1
            continue
        break
    return "\n".join(reversed(comments))


def attached_attrs(lines: list[str], fn_index: int) -> list[str]:
    attrs: list[str] = []
    idx = fn_index - 1
    while idx >= 0:
        stripped = lines[idx].strip()
        if stripped == "":
            idx -= 1
            continue
        if stripped.startswith("#["):
            attrs.append(stripped)
            idx -= 1
            continue
        if stripped.startswith(("///", "//", "//!")):
            idx -= 1
            continue
        break
    return list(reversed(attrs))


def test_module_lines(lines: list[str]) -> set[int]:
    test_lines: set[int] = set()
    depth: int | None = None
    pending_cfg_test = False
    for idx, line in enumerate(lines):
        stripped = line.strip()
        if stripped == "#[cfg(test)]":
            pending_cfg_test = True
        elif pending_cfg_test and re.match(r"(?:pub\s+)?mod\s+\w+\s*\{", stripped):
            depth = 0
            pending_cfg_test = False
        elif stripped and not stripped.startswith("#["):
            pending_cfg_test = False

        if depth is not None:
            test_lines.add(idx)
            depth += line.count("{") - line.count("}")
            if depth <= 0:
                depth = None
    return test_lines


def clean_ref(ref: str) -> str | None:
    if ref.endswith((".c", ".h")) or ".cpp::" in ref or "::" in ref:
        return None
    if not re.match(r"^[A-Za-z_][A-Za-z0-9_]*$", ref):
        return None
    return ref


def extract_ncbi_refs(comment: str) -> tuple[list[str], bool]:
    refs: list[str] = []
    standard = False
    for match in STANDARD_NCBI_RE.finditer(comment):
        standard = True
        cleaned = clean_ref(match.group(1))
        if cleaned:
            refs.append(cleaned)
    for regex in LEGACY_REF_RES:
        for ref in regex.findall(comment):
            cleaned = clean_ref(ref)
            if cleaned:
                refs.append(cleaned)
    return sorted(set(refs)), standard


def is_native(comment: str) -> bool:
    return bool(NATIVE_RE.search(comment))


def has_naming_note(comment: str) -> bool:
    return bool(NAMING_NOTE_RE.search(comment))


def c_to_rust_name(name: str) -> str:
    name = name.strip("_")
    name = re.sub(r"([A-Z]+)([A-Z][a-z])", r"\1_\2", name)
    name = re.sub(r"([a-z0-9])([A-Z])", r"\1_\2", name)
    name = name.replace("__", "_")
    snake = re.sub(r"[^A-Za-z0-9]+", "_", name).strip("_").lower()
    # Common BLAST acronyms should stay as one Rust token.
    snake = snake.replace("hs_ps", "hsps")
    snake = snake.replace("hs_p", "hsp")
    snake = snake.replace("ps_sm", "pssm")
    snake = snake.replace("rp_s", "rps")
    snake = snake.replace("ph_i", "phi")
    return snake


def scan_source(root: Path) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for path in sorted(root.rglob("*.rs")):
        lines = path.read_text(encoding="utf-8", errors="ignore").splitlines()
        test_lines = test_module_lines(lines)
        for idx, line in enumerate(lines):
            match = FN_RE.match(line.strip())
            if not match:
                continue
            attrs = attached_attrs(lines, idx)
            if idx in test_lines or any(attr == "#[test]" for attr in attrs):
                continue
            comment = attached_comment(lines, idx)
            refs, standard_ncbi = extract_ncbi_refs(comment)
            native = is_native(comment)
            naming_note = has_naming_note(comment)
            rows.append(
                {
                    "file": str(path),
                    "line": idx + 1,
                    "rust": match.group(1),
                    "refs": refs,
                    "standard_ncbi": standard_ncbi,
                    "native": native,
                    "naming_note": naming_note,
                    "tagged": bool(refs or native),
                }
            )
    return rows


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--src", default="src", type=Path)
    parser.add_argument(
        "--file",
        action="append",
        type=Path,
        help="Limit output to one Rust source file. May be passed multiple times.",
    )
    parser.add_argument("--rust-report", default=".tmp/ccc-rust.json", type=Path)
    parser.add_argument("--other-report", default=".tmp/ccc-ncbi-core.json", type=Path)
    parser.add_argument("--json", action="store_true")
    parser.add_argument("--limit", default=40, type=int)
    args = parser.parse_args()

    rows = scan_source(args.src)
    if args.file:
        wanted = {str(path) for path in args.file}
        rows = [row for row in rows if row["file"] in wanted]
    rust_report_names = load_function_names(args.rust_report)
    other_report_names = load_function_names(args.other_report)

    tagged = [row for row in rows if row["tagged"]]
    ncbi_rows = [row for row in rows if row["refs"]]
    native_rows = [row for row in rows if row["native"]]
    untagged = [row for row in rows if not row["tagged"]]
    standard_ncbi = [row for row in ncbi_rows if row["standard_ncbi"]]

    ref_counts = Counter(ref for row in ncbi_rows for ref in row["refs"])
    missing_other = sorted(
        {
            ref
            for row in ncbi_rows
            for ref in row["refs"]
            if other_report_names and ref not in other_report_names
        }
    )
    missing_rust = sorted(
        str(row["rust"]) for row in tagged if rust_report_names and row["rust"] not in rust_report_names
    )
    duplicate_refs = sorted((ref, count) for ref, count in ref_counts.items() if count > 1)
    naming_mismatches = []
    for row in ncbi_rows:
        if len(row["refs"]) != 1:
            continue
        if row["naming_note"]:
            continue
        expected = c_to_rust_name(str(row["refs"][0]))
        if row["rust"] != expected:
            naming_mismatches.append((str(row["file"]), int(row["line"]), str(row["rust"]), row["refs"][0], expected))

    if args.json:
        print(json.dumps(rows, indent=2, sort_keys=True))
        return 0

    print(f"Rust functions scanned: {len(rows)}")
    print(f"Functions with provenance tag: {len(tagged)}")
    print(f"Functions with NCBI refs: {len(ncbi_rows)}")
    print(f"Functions using standard NCBI: tag: {len(standard_ncbi)}")
    print(f"Functions marked blast-rs native: {len(native_rows)}")
    print(f"Functions missing provenance tag: {len(untagged)}")
    for row in untagged[: args.limit]:
        print(f"  - {row['file']}:{row['line']} {row['rust']}")
    if len(untagged) > args.limit:
        print(f"  ... {len(untagged) - args.limit} more")

    print(f"Referenced NCBI names absent from original function report: {len(missing_other)}")
    for name in missing_other[: args.limit]:
        print(f"  - {name}")

    print(f"Tagged Rust names absent from Rust function report: {len(missing_rust)}")
    print(f"Duplicate referenced NCBI names: {len(duplicate_refs)}")
    for name, count in duplicate_refs[: args.limit]:
        print(f"  - {name} ({count})")

    print(f"NCBI-tagged naming mismatches: {len(naming_mismatches)}")
    for file, line, rust, other, expected in naming_mismatches[: args.limit]:
        print(f"  - {file}:{line} {rust} -> {other}; expected {expected}")

    by_file: dict[str, int] = defaultdict(int)
    for row in untagged:
        by_file[str(row["file"])] += 1
    print("Top untagged files:")
    for path, count in sorted(by_file.items(), key=lambda item: (-item[1], item[0]))[: args.limit]:
        print(f"  - {path}: {count}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
