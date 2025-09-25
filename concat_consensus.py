#!/usr/bin/env python3
"""
Concatenate multiple consensus FASTA files into a single multi-FASTA.

Typical usage
-------------
python concat_consensus_fasta.py \
  --in-dir /path/to/consensus \
  --pattern "*_consensus.fasta" \
  --out all_consensus.fasta

Options
-------
--recursive               Recurse into subdirectories (default: off)
--use-header              Keep record IDs as-is from each FASTA header (default: derive from filename)
--sample-regex REGEX      Regex with (?P<sample>...) to extract sample name from filename
                          Default: '^(?P<sample>[^_]+)_' (text before first underscore)
--fail-on-duplicate       Error if duplicate IDs are encountered (default: auto-deduplicate by appending _1, _2, ...)
--mapping TSV             Optional path to write a mapping table: filename → final_id → length

Dependencies: Biopython
  pip install biopython
"""
from __future__ import annotations
import argparse
import sys
from pathlib import Path
import re
from typing import Dict, List, Tuple

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def infer_sample_from_filename(name: str, regex: str) -> str:
    m = re.search(regex, name)
    if m and 'sample' in m.groupdict():
        return m.group('sample')
    # Fallbacks: before _consensus, then before first underscore, then stem
    stem = name.rsplit('.', 1)[0]
    if stem.endswith('_consensus'):
        stem = stem[:-10]
    return stem.split('_')[0] if '_' in stem else stem


def unique_id(base: str, seen: Dict[str, int], fail_on_duplicate: bool) -> str:
    if base not in seen:
        seen[base] = 0
        return base
    if fail_on_duplicate:
        raise SystemExit(f"Duplicate ID encountered: {base}. Use a different --sample-regex or omit --fail-on-duplicate to auto-deduplicate.")
    seen[base] += 1
    return f"{base}_{seen[base]}"


def collect_fastas(in_dir: Path, pattern: str, recursive: bool) -> List[Path]:
    if recursive:
        return sorted(p for p in in_dir.rglob(pattern))
    else:
        return sorted(p for p in in_dir.glob(pattern))


def main():
    ap = argparse.ArgumentParser(description="Concatenate consensus FASTA files into one multi-FASTA")
    ap.add_argument('--in_dir', type=Path, required=True, help='Directory containing consensus FASTA files')
    ap.add_argument('--pattern', default='*_consensus.fasta', help='Glob pattern for input FASTA files (default: *_consensus.fasta)')
    ap.add_argument('--out', type=Path, required=True, help='Output multi-FASTA path')
    ap.add_argument('--recursive', action='store_true', help='Recurse into subdirectories')

    ap.add_argument('--use-header', action='store_true', help='Keep IDs from FASTA headers; otherwise derive from filename')
    ap.add_argument('--sample-regex', default=r'^(?P<sample>[^_]+)_', help='Regex with (?P<sample>...) to extract sample name (default: before first underscore)')
    ap.add_argument('--fail-on-duplicate', action='store_true', help='Fail if duplicate IDs are found')
    ap.add_argument('--mapping', type=Path, help='Write TSV mapping: filename\tfinal_id\tlength')

    args = ap.parse_args()

    files = collect_fastas(args.in_dir, args.pattern, args.recursive)
    if not files:
        raise SystemExit('No FASTA files found. Check --in-dir and --pattern.')

    args.out.parent.mkdir(parents=True, exist_ok=True)

    seen: Dict[str, int] = {}
    mapping_rows: List[str] = []
    all_records: List[SeqRecord] = []

    for fp in files:
        # Expect one record per file; if more, include all
        records = list(SeqIO.parse(str(fp), 'fasta'))
        if not records:
            print(f"[WARN] Empty FASTA: {fp}", file=sys.stderr)
            continue
        for i, rec in enumerate(records):
            if args.use_header:
                base_id = rec.id
            else:
                base_id = infer_sample_from_filename(fp.name, args.sample_regex)
                # If file has multiple records, suffix with index for stability
                if len(records) > 1:
                    base_id = f"{base_id}_{i+1}"
            final_id = unique_id(base_id, seen, args.fail_on_duplicate)
            # Normalize: strip description to keep clean headers
            rec.id = final_id
            rec.name = final_id
            rec.description = ''
            all_records.append(rec)
            mapping_rows.append(f"{fp.name}\t{final_id}\t{len(rec.seq)}")

    with open(args.out, 'w') as out_fh:
        SeqIO.write(all_records, out_fh, 'fasta')

    if args.mapping:
        with open(args.mapping, 'w') as mfh:
            mfh.write("filename\tfinal_id\tlength\n")
            mfh.write("\n".join(mapping_rows) + "\n")

    print(f"Wrote {len(all_records)} records → {args.out}")
    if args.mapping:
        print(f"Mapping table → {args.mapping}")

if __name__ == '__main__':
    main()
