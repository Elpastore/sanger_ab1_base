#!/usr/bin/env python3
"""
AB1 Forward+Reverse Preprocessing → Consensus FASTA (paired across directories)

Features
- Read Sanger chromatograms (.ab1) for each sample from separate forward and reverse directories
- Infer sample name automatically (default: text before first underscore) or via a user regex with a named group (?P<sample>...)
- Quality trimming (sliding window, Q20 default) and masking of low-Q bases to 'N'
- Reverse-complement the reverse read
- Global alignment of FWD vs REV_RC using Bio.Align.PairwiseAligner (pairwise2 deprecated)
- Consensus by quality/IUPAC, with robust fallbacks when one read is too short/empty
- Output {sample}_consensus.fasta to --out-dir and a summary report.tsv

Examples
  python ab1_consensus_dirpair.py \
    --fwd-dir ./FWD --rev-dir ./REV --out-dir ./consensus

  # If your sample name pattern is custom, e.g., "^SMP(?P<sample>\\d+)_" → SMP12_...
  python ab1_consensus_dirpair.py \
    --fwd-dir ./FWD --rev-dir ./REV --out-dir ./consensus \
    --sample-regex '^(?P<sample>[^_]+)_'  # default; captures before first underscore

Dependencies: Biopython (pip install biopython)
"""
from __future__ import annotations
import argparse
import re
from pathlib import Path
from typing import Dict, List, Tuple, Optional

from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# --------------------- IUPAC helpers ---------------------
IUPAC_TO_SET: Dict[str, set] = {
    'A': {'A'}, 'C': {'C'}, 'G': {'G'}, 'T': {'T'}, 'U': {'T'},
    'R': {'A','G'}, 'Y': {'C','T'}, 'S': {'G','C'}, 'W': {'A','T'},
    'K': {'G','T'}, 'M': {'A','C'},
    'B': {'C','G','T'}, 'D': {'A','G','T'}, 'H': {'A','C','T'}, 'V': {'A','C','G'},
    'N': {'A','C','G','T'}, '-': set()
}
SET_TO_IUPAC = {frozenset(v): k for k, v in IUPAC_TO_SET.items()}

# --------------------- IO helpers ---------------------
def read_ab1(path: Path) -> Tuple[str, List[int]]:
    rec = SeqIO.read(str(path), 'abi')
    seq = str(rec.seq).upper()
    quals = rec.letter_annotations.get('phred_quality', [20]*len(seq))
    return seq, list(quals)

# --------------------- Trimming & masking ---------------------
def sliding_trim(seq: str, qual: List[int], qmin: int = 20, window: int = 10, min_len: int = 50) -> Tuple[str, List[int]]:
    if not seq:
        return seq, qual
    n = len(seq)
    # Left
    left = 0
    while left < n:
        w = qual[left:left+window]
        if sum(w) / max(1, len(w)) >= qmin:
            break
        left += 1
    # Right
    right = n
    while right > left:
        start = max(left, right-window)
        w = qual[start:right]
        if sum(w) / max(1, len(w)) >= qmin:
            break
        right -= 1
    seq2, qual2 = seq[left:right], qual[left:right]
    if len(seq2) < min_len:
        L = 0
        while L < n and qual[L] < qmin:
            L += 1
        R = n
        while R > L and qual[R-1] < qmin:
            R -= 1
        seq2, qual2 = seq[L:R], qual[L:R]
    return seq2, qual2

def quality_mask(seq: str, qual: List[int], qmin: int = 20) -> str:
    return ''.join((b if q >= qmin else 'N') for b, q in zip(seq, qual))

# --------------------- Alignment helpers ---------------------
def _reconstruct_gapped_strings_from_alignment(aln, s1: str, s2: str) -> Tuple[str, str]:
    """Reconstruct gapped strings from a PairwiseAligner alignment object."""
    a_gapped = []
    b_gapped = []
    i, j = 0, 0
    for (i0, i1), (j0, j1) in zip(aln.aligned[0], aln.aligned[1]):
        # gaps in s1 before next block
        while i < i0:
            a_gapped.append(s1[i]); b_gapped.append('-'); i += 1
        # gaps in s2 before next block
        while j < j0:
            a_gapped.append('-'); b_gapped.append(s2[j]); j += 1
        # matched block
        a_gapped.append(s1[i0:i1])
        b_gapped.append(s2[j0:j1])
        i, j = i1, j1
    # tail gaps
    while i < len(s1):
        a_gapped.append(s1[i]); b_gapped.append('-'); i += 1
    while j < len(s2):
        a_gapped.append('-'); b_gapped.append(s2[j]); j += 1
    return ''.join(a_gapped), ''.join(b_gapped)


def align_global(seq1: str, seq2: str) -> Tuple[str, str]:
    """Global alignment (Needleman–Wunsch-like) using Bio.Align.PairwiseAligner.
    Returns (aligned_seq1, aligned_seq2). Raises ValueError if no alignment.
    """
    aligner = PairwiseAligner()
    aligner.mode = 'global'
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    it = aligner.align(seq1, seq2)
    try:
        aln = next(iter(it))
    except StopIteration:
        raise ValueError('No alignment produced (sequences empty/too short?)')
    return _reconstruct_gapped_strings_from_alignment(aln, seq1, seq2)


def align_local(seq1: str, seq2: str) -> Tuple[str, str]:
    """Local alignment fallback using PairwiseAligner."""
    aligner = PairwiseAligner()
    aligner.mode = 'local'
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    it = aligner.align(seq1, seq2)
    try:
        aln = next(iter(it))
    except StopIteration:
        raise ValueError('No local alignment produced')
    return _reconstruct_gapped_strings_from_alignment(aln, seq1, seq2)


def expand_quals(aln_seq: str, raw_qual: List[int]) -> List[int]:
    out: List[int] = []
    i = 0
    for c in aln_seq:
        if c == '-':
            out.append(-1)
        else:
            out.append(raw_qual[i])
            i += 1
    return out

# --------------------- Consensus ---------------------
def consensus_from_alignment(a: str, b: str, qa: List[int], qb: List[int], qmin: int = 20) -> str:
    cons = []
    for x, y, qx, qy in zip(a, b, qa, qb):
        if x == '-' and y == '-':
            continue
        if x == '-':
            cons.append(y if (qy >= qmin and y != '-') else 'N'); continue
        if y == '-':
            cons.append(x if (qx >= qmin and x != '-') else 'N'); continue
        if x == y:
            cons.append(x if (qx >= qmin or qy >= qmin) else 'N'); continue
        if qx >= qy and qx >= qmin:
            cons.append(x); continue
        if qy > qx and qy >= qmin:
            cons.append(y); continue
        # both low: emit ambiguity if possible
        sx = IUPAC_TO_SET.get(x, {x}); sy = IUPAC_TO_SET.get(y, {y})
        cons.append(SET_TO_IUPAC.get(frozenset(sx | sy), 'N'))
    # trim leading/trailing Ns
    s = ''.join(cons)
    s = s.lstrip('N').rstrip('N')
    return s

# --------------------- Pairing logic ---------------------
def extract_sample_name(filename: str, sample_regex: Optional[str]) -> Optional[str]:
    base = filename
    if sample_regex:
        m = re.search(sample_regex, base)
        if m and 'sample' in m.groupdict():
            return m.group('sample')
        return None
    # default: take text before first underscore
    return base.split('_')[0] if '_' in base else base.rsplit('.', 1)[0]

# --------------------- Main pipeline per pair ---------------------
MIN_ALIGN_LEN = 30  # minimal non-N length to attempt alignment

def build_consensus(forward_path: Path, reverse_path: Path, qmin: int, window: int, min_len: int) -> Tuple[str, str]:
    fseq, fq = read_ab1(forward_path)
    rseq, rq = read_ab1(reverse_path)

    fseq_t, fq_t = sliding_trim(fseq, fq, qmin=qmin, window=window, min_len=min_len)
    rseq_t, rq_t = sliding_trim(rseq, rq, qmin=qmin, window=window, min_len=min_len)

    fseq_m = quality_mask(fseq_t, fq_t, qmin=qmin)
    rseq_m = quality_mask(rseq_t, rq_t, qmin=qmin)

    rseq_rc = str(Seq(rseq_m).reverse_complement())
    rq_rc = list(reversed(rq_t))

    # Pre-alignment guards
    nonN_f = len(fseq_m.replace('N',''))
    nonN_r = len(rseq_rc.replace('N',''))
    if nonN_f < MIN_ALIGN_LEN and nonN_r >= MIN_ALIGN_LEN:
        cons = rseq_rc
        stats = f"FWD too short; used REV only | len_f={len(fseq)} len_r={len(rseq)} len_cons={len(cons)}"
        return cons, stats
    if nonN_r < MIN_ALIGN_LEN and nonN_f >= MIN_ALIGN_LEN:
        cons = fseq_m
        stats = f"REV too short; used FWD only | len_f={len(fseq)} len_r={len(rseq)} len_cons={len(cons)}"
        return cons, stats
    if nonN_f < MIN_ALIGN_LEN and nonN_r < MIN_ALIGN_LEN:
        cons = fseq_m if len(fseq_m) >= len(rseq_rc) else rseq_rc
        stats = f"Both too short; returned longer masked read | len_f={len(fseq)} len_r={len(rseq)} len_cons={len(cons)}"
        return cons, stats

    # Align (global, then fallback to local)
    try:
        a, b = align_global(fseq_m, rseq_rc)
    except ValueError:
        try:
            a, b = align_local(fseq_m, rseq_rc)
        except ValueError:
            cons = fseq_m if len(fseq_m) >= len(rseq_rc) else rseq_rc
            stats = f"Global+local failed; fallback to longer strand | len_f={len(fseq)} len_r={len(rseq)} len_cons={len(cons)}"
            return cons, stats

    qa = expand_quals(a, fq_t)
    qb = expand_quals(b, rq_rc)

    cons = consensus_from_alignment(a, b, qa, qb, qmin=qmin)
    stats = f"len_f={len(fseq)} len_r={len(rseq)} len_cons={len(cons)}"
    return cons, stats

# --------------------- CLI ---------------------

def main():
    ap = argparse.ArgumentParser(description='Preprocess AB1 from separate forward/reverse directories into consensus FASTA per sample')
    ap.add_argument('--fwd-dir', type=Path, required=True, help='Directory containing forward .ab1 files')
    ap.add_argument('--rev-dir', type=Path, required=True, help='Directory containing reverse .ab1 files')
    ap.add_argument('--out-dir', type=Path, required=True, help='Output directory for consensus FASTA and report.tsv')
    ap.add_argument('--fwd-glob', default='*.ab1', help='Glob pattern for forward files (default: *.ab1)')
    ap.add_argument('--rev-glob', default='*.ab1', help='Glob pattern for reverse files (default: *.ab1)')
    ap.add_argument('--sample-regex', default=r'^(?P<sample>[^_]+)_', help='Regex with (?P<sample>...) to extract sample name (default: before first underscore)')

    ap.add_argument('--qmin', type=int, default=20, help='Minimum Phred; mask to N below this (default: 20)')
    ap.add_argument('--window', type=int, default=10, help='Sliding window size for trimming (default: 10)')
    ap.add_argument('--min-len', type=int, default=50, help='Minimum kept length after trimming (default: 50)')

    args = ap.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)

    # Map sample -> file path for forward and reverse
    fwd_map: Dict[str, Path] = {}
    rev_map: Dict[str, Path] = {}

    for p in args.fwd_dir.glob(args.fwd_glob):
        s = extract_sample_name(p.name, args.sample_regex)
        if s:
            prev = fwd_map.get(s)
            if not prev or ('Fwd' in p.name and 'Fwd' not in prev.name):
                fwd_map[s] = p
    for p in args.rev_dir.glob(args.rev_glob):
        s = extract_sample_name(p.name, args.sample_regex)
        if s:
            prev = rev_map.get(s)
            if not prev or ('Rev' in p.name and 'Rev' not in prev.name):
                rev_map[s] = p

    samples = sorted(set(fwd_map.keys()) & set(rev_map.keys()))
    if not samples:
        raise SystemExit('No overlapping samples between forward and reverse directories. Check --sample-regex and filenames.')

    report_lines = ['sample\tforward\treverse\tconsensus_fasta\tlen_cons']
    for s in samples:
        fwd = fwd_map[s]
        rev = rev_map[s]
        cons, stats = build_consensus(fwd, rev, qmin=args.qmin, window=args.window, min_len=args.min_len)
        out_fa = args.out_dir / f'{s}_consensus.fasta'
        rec = SeqRecord(Seq(cons), id=s, description='')
        with open(out_fa, 'w') as fh:
            SeqIO.write(rec, fh, 'fasta')
        print(f'{s}: {stats} -> {out_fa}')
        report_lines.append(f'{s}\t{fwd.name}\t{rev.name}\t{out_fa.name}\t{len(cons)}')

    with open(args.out_dir / 'report.tsv', 'w') as fh:
        fh.write('\n'.join(report_lines) + '\n')
    print(f'Summary report -> {args.out_dir / "report.tsv"}')

if __name__ == '__main__':
    main()
