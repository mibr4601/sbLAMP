import argparse
import bisect
import heapq
import os
import sys
import time
from dataclasses import dataclass, replace
import numpy as np
import primer3 as _p3

from single import (
    read_fasta,
    validate_alignment,
    reverse_complement,
    load_params,
    _watson_crick,
    _col_to_real,
    _real_to_col,
    _pair_ok,
    _P3_HOM1,
)
import repairSet

_COV_SCALE = 0.4
_FAIL_REPAIR = 4.0


def _target_credit(total_mismatches):
    # Assumed more than 8 differences is total mismatch
    return max(
        0.0,
        1.0 - min(total_mismatches, 8) / 8,
    )


_fail_pairs = set()
_pass_pairs = set()
_PASS_CAP = 50_000  # This is for the cache with secondary structure checking, not sets.


def _cached_check(key, compute):
    if key in _fail_pairs:
        return False
    if key in _pass_pairs:
        return True
    if not compute():
        _fail_pairs.add(key)
        return False
    if len(_pass_pairs) >= _PASS_CAP:
        _pass_pairs.clear()
    _pass_pairs.add(key)
    return True


_CORE_ROLES = ["F3", "F2", "F1c", "B1c", "B2", "B3"]


@dataclass
class Candidate:
    pos: int
    length: int
    rev: int
    seq: str
    regGC: float
    GC: float
    Tm: float
    dG5: float
    dG3: float
    score: float = 0.0
    rc_seq: str = ""
    src: list = None
    cov: float = 0.0


def _parse_candidate(line):
    f = {}
    for part in line.strip().split("\t"):
        k, v = part.split(":", 1)
        f[k] = v
    raw_src = f.get("src")
    return Candidate(
        pos=int(f["pos"]),
        length=int(f["len"]),
        rev=int(f["rev"]),
        seq=f["seq"],
        regGC=float(f["regGC"]),
        GC=float(f["GC"]),
        Tm=float(f["Tm"]),
        dG5=float(f["dG5"]),
        dG3=float(f["dG3"]),
        score=float(f.get("score", 0.0)),
        src=raw_src.split(",") if raw_src else None,
        cov=float(f.get("cov", 0.0)),
    )


def load_candidates(path):
    cands = []
    if not os.path.exists(path):
        return cands
    with open(path) as fh:
        for line in fh:
            if line.strip():
                c = _parse_candidate(line)
                c.rc_seq = reverse_complement(c.seq)
                cands.append(c)
    cands.sort(key=lambda c: c.pos)
    return cands


_LINKER = "TTTT"


# Default gBlock synthesis-constraint criteria.
# Criteria can be changed to accomodate especially for GC-rich sequences
DEFAULT_GBLOCK = {
    "wide_window": 100,  # no window this wide can exceed wide_gc%
    "wide_gc": 76.0,
    "poly_at": 10,
    "poly_gc": 7,
    "len_min": 240,
    "len_max": 300,
    "gc_max": 70.0,  # overall gBlock GC ceiling
    "tail_win": 60,
    "tail_gc": 72.0,
    "pad_max": 60,
    "pad_min": 10,
}


def _gc_percent_windows(sub, window):
    n = len(sub)
    w = min(window, n) if n else 0
    if w == 0:
        return []
    gc_count = sub[:w].count("G") + sub[:w].count("C")
    out = [100.0 * gc_count / w]
    for i in range(w, n):
        if sub[i - w] in "GC":
            gc_count -= 1
        if sub[i] in "GC":
            gc_count += 1
        out.append(100.0 * gc_count / w)
    return out


def _max_homopolymer_runs(sub):
    max_at, max_gc = 0, 0
    n = len(sub)
    i = 0
    while i < n:
        j = i
        while j < n and sub[j] == sub[i]:
            j += 1
        run_len = j - i
        if sub[i] in "AT":
            max_at = max(max_at, run_len)
        elif sub[i] in "GC":
            max_gc = max(max_gc, run_len)
        i = j
    return max_at, max_gc


def _gblock_fail_reason(seq, start, end):
    sub = seq[max(0, start) : min(len(seq), end)].upper()
    n = len(sub)
    if n == 0:
        return None

    if (
        max(_gc_percent_windows(sub, DEFAULT_GBLOCK["wide_window"]), default=0.0)
        > DEFAULT_GBLOCK["wide_gc"]
    ):
        return "wide_window"
    if max(_gc_percent_windows(sub, 20), default=0.0) > 90.0:
        return "tight_window"

    window_gc = _gc_percent_windows(sub, 50)
    if window_gc and max(window_gc) - min(window_gc) > 50.0:
        return "spread"

    max_at, max_gc_run = _max_homopolymer_runs(sub)
    if max_at > DEFAULT_GBLOCK["poly_at"] or max_gc_run > DEFAULT_GBLOCK["poly_gc"]:
        return "homopolymer"

    return None


def _region_gc(seq, start, end):
    n = end - start
    if n <= 0:
        return 0.0
    sub = seq[start:end]
    return 100.0 * (sub.count("G") + sub.count("C")) / n


def optimize_gblock(
    seq,
    amp_start,
    amp_end,
    lo_limit=0,
    hi_limit=None,
):
    """Checks room available on each side, constraints and picks survivor with GC% closest to 50%"""
    pad_max = DEFAULT_GBLOCK["pad_max"]
    pad_min = DEFAULT_GBLOCK["pad_min"]
    len_min = DEFAULT_GBLOCK["len_min"]
    len_max = DEFAULT_GBLOCK["len_max"]
    tail_win = DEFAULT_GBLOCK["tail_win"]
    tail_gc = DEFAULT_GBLOCK["tail_gc"]

    if hi_limit is None:
        hi_limit = len(seq)
    max_head = min(pad_max, amp_start - lo_limit)
    max_tail = min(pad_max, hi_limit - amp_end)
    if max_head < pad_min or max_tail < pad_min:
        return None

    lo_bound = amp_start - max_head
    hi_bound = amp_end + max_tail
    window = np.frombuffer(seq[lo_bound:hi_bound].upper().encode(), dtype=np.uint8)
    is_gc = (window == ord("G")) | (window == ord("C"))
    prefix = np.zeros(is_gc.size + 1, dtype=np.int32)
    np.cumsum(is_gc, out=prefix[1:])

    heads = np.arange(pad_min, max_head + 1)
    tails = np.arange(pad_min, max_tail + 1)
    s_vals = amp_start - heads
    e_vals = amp_end + tails
    s_idx = s_vals - lo_bound
    e_idx = e_vals - lo_bound

    length = e_idx[None, :] - s_idx[:, None]
    valid = (length >= len_min) & (length <= len_max)

    tail_start_idx = np.maximum(e_idx - tail_win, 0)
    tail_gc_pct = 100.0 * (prefix[e_idx] - prefix[tail_start_idx]) / tail_win
    valid &= (tail_gc_pct <= tail_gc)[None, :]

    if not valid.any():
        return None

    gc_pct = 100.0 * (prefix[e_idx][None, :] - prefix[s_idx][:, None]) / length
    score = np.where(valid, np.abs(gc_pct - 50.0), np.inf)
    i, j = np.unravel_index(np.argmin(score), score.shape)
    return (int(s_vals[i]), int(e_vals[j]))


def check_add(f1c_pos, accepted_sorted, seq_len):
    """Checks that distance between accepted f1c is large enough"""
    min_dist = 150 if seq_len > 50000 else 20
    i = bisect.bisect_left(accepted_sorted, f1c_pos)
    if i < len(accepted_sorted) and accepted_sorted[i] - f1c_pos < min_dist:
        return False
    if i > 0 and f1c_pos - accepted_sorted[i - 1] < min_dist:
        return False
    return True


def check_structure_loop(
    cands,
    lf,
    lb,
    gc,
):
    threshold = 49 if gc in (1, 4) else 45
    extras = [c for c in (lf, lb) if c is not None]
    all_cands = cands + extras

    for i in range(len(all_cands)):
        for j in range(i + 1, len(all_cands)):
            ci, cj = all_cands[i], all_cands[j]
            if not _cached_check(
                (ci.seq, cj.seq),
                lambda ci=ci, cj=cj: _pair_ok(
                    ci.seq, cj.seq, threshold, ci.rc_seq, cj.rc_seq
                ),
            ):
                return False

    fip = cands[2].seq + _LINKER + cands[1].seq
    bip = cands[3].seq + _LINKER + cands[4].seq

    def _hairpin_ok(seq):
        r = _p3.calc_hairpin(seq, **_P3_HOM1)
        return not (r.structure_found and r.tm > threshold)

    if not _cached_check((fip,), lambda: _hairpin_ok(fip)):
        return False
    if not _cached_check((bip,), lambda: _hairpin_ok(bip)):
        return False
    if not _cached_check((fip, bip), lambda: _pair_ok(fip, bip, threshold)):
        return False

    return True


def _src_masks(cands, registry):
    """Bitmask ecoding of which real input sequences is the primer compatible"""
    masks = np.empty(len(cands), dtype=object)
    for i, c in enumerate(cands):
        if not c.src:
            masks[i] = -1
            continue
        m = 0
        for name in c.src:
            if name not in registry:
                registry[name] = len(registry)
            m |= 1 << registry[name]
        masks[i] = m
    return masks


def _combo_chunk_gen(
    inner,
    outer,
    seq,
    use_loop,
    loop_cands=[],
):
    inner_plus = [c for c in inner if c.rev == 1]
    inner_minus = [c for c in inner if c.rev == -1]
    outer_plus = [c for c in outer if c.rev == 1]
    outer_minus = [c for c in outer if c.rev == -1]

    if not (inner_plus and inner_minus and outer_plus and outer_minus):
        return

    def _attr_array(lst, attr, dt=np.int32):
        return np.array([getattr(c, attr) for c in lst], dtype=dt)

    _src_registry = {}

    # Inner minus details
    im_pos = _attr_array(inner_minus, "pos")
    im_len = _attr_array(inner_minus, "length")
    im_end = im_pos + im_len
    im_tm = _attr_array(inner_minus, "Tm", np.float32)
    im_gc = _attr_array(inner_minus, "GC", np.float32)
    im_score = _attr_array(inner_minus, "score", np.float32)
    im_src = _src_masks(inner_minus, _src_registry)

    # Inner plus details
    ip_pos = _attr_array(inner_plus, "pos")
    ip_len = _attr_array(inner_plus, "length")
    ip_end = ip_pos + ip_len
    ip_tm = _attr_array(inner_plus, "Tm", np.float32)
    ip_gc = _attr_array(inner_plus, "GC", np.float32)
    ip_score = _attr_array(inner_plus, "score", np.float32)
    ip_src = _src_masks(inner_plus, _src_registry)

    # Outer plus details
    op_pos = _attr_array(outer_plus, "pos")
    op_len = _attr_array(outer_plus, "length")
    op_end = op_pos + op_len
    op_tm = _attr_array(outer_plus, "Tm", np.float32)
    op_gc = _attr_array(outer_plus, "GC", np.float32)
    op_score = _attr_array(outer_plus, "score", np.float32)
    op_src = _src_masks(outer_plus, _src_registry)
    op_end_ord = np.argsort(op_end, kind="stable")
    op_end_sorted = op_end[op_end_ord]
    op_tm_sorted = op_tm[op_end_ord]

    # Outer minus
    om_pos = _attr_array(outer_minus, "pos")
    om_end = om_pos + _attr_array(outer_minus, "length")
    om_tm = _attr_array(outer_minus, "Tm", np.float32)
    om_gc = _attr_array(outer_minus, "GC", np.float32)
    om_score = _attr_array(outer_minus, "score", np.float32)
    om_src = _src_masks(outer_minus, _src_registry)
    om_end_ord = np.argsort(om_end, kind="stable")
    om_end_sorted = om_end[om_end_ord]
    om_pos_ord = np.argsort(om_pos, kind="stable")
    om_pos_sorted = om_pos[om_pos_ord]

    # Loop forward and backwards
    if use_loop and loop_cands:
        lf_cands = [c for c in loop_cands if c.rev == -1]
        lb_cands = [c for c in loop_cands if c.rev == 1]
        lf_pos_arr = (
            np.array([c.pos for c in lf_cands], dtype=np.int32)
            if lf_cands
            else np.empty(0, np.int32)
        )
        lb_pos_arr = (
            np.array([c.pos for c in lb_cands], dtype=np.int32)
            if lb_cands
            else np.empty(0, np.int32)
        )
    else:
        lf_pos_arr = lb_pos_arr = np.empty(0, np.int32)

    seq_b = seq.upper().encode()
    raw = np.frombuffer(seq_b, dtype=np.uint8)
    is_gc = (raw == ord("G")) | (raw == ord("C"))
    gc_pfx = np.zeros(len(seq) + 1, dtype=np.int32)
    np.cumsum(is_gc, out=gc_pfx[1:])

    def _gcb_arr(gc_arr):
        d = np.maximum(0.0, np.maximum(50.0 - gc_arr, gc_arr - 60.0))
        return (1.0 / (1.0 + (d / 10.0) ** 2)).astype(np.float32)

    im_gcb = _gcb_arr(im_gc)
    ip_gcb = _gcb_arr(ip_gc)
    op_gcb = _gcb_arr(op_gc)
    om_gcb = _gcb_arr(om_gc)

    max_expand = 30_000_000

    def _expand(lo_arr, hi_arr):
        counts = (hi_arr.astype(np.int64) - lo_arr.astype(np.int64)).clip(min=0)
        if counts.sum() > max_expand:
            hi_arr = np.minimum(hi_arr, lo_arr + max(1, max_expand // len(lo_arr)))
            counts = (hi_arr.astype(np.int64) - lo_arr.astype(np.int64)).clip(min=0)
        total = int(counts.sum())
        if total == 0:
            empty = np.empty(0, dtype=np.int32)
            return empty, empty
        par = np.repeat(np.arange(len(lo_arr), dtype=np.int32), counts)
        cum = np.empty(len(lo_arr) + 1, dtype=np.int64)
        cum[0] = 0
        np.cumsum(counts, out=cum[1:])
        k = np.arange(total, dtype=np.int64)
        idx = (lo_arr[par].astype(np.int64) + k - cum[par]).astype(np.int32)
        return par, idx

    _f3_lo = np.searchsorted(op_end_sorted, op_pos - 19)
    _f3_hi = np.searchsorted(op_end_sorted, op_pos + 1)
    op_has_f3 = _f3_hi > _f3_lo

    _b3_lo = np.searchsorted(om_pos_sorted, om_end)
    _b3_hi = np.searchsorted(om_pos_sorted, om_end + 20)
    om_has_b3 = _b3_hi > _b3_lo

    b1c_lo_all = np.searchsorted(ip_pos, im_end + 1)
    b1c_hi_all = np.searchsorted(ip_pos, im_pos + 60)
    im_has_b1c = b1c_hi_all > b1c_lo_all

    f2_lo_all = np.searchsorted(op_pos, im_pos - 60)
    f2_hi_all = np.searchsorted(op_pos, im_pos - 39)
    op_f3_pfx = np.zeros(len(op_pos) + 1, dtype=np.int32)
    np.cumsum(op_has_f3.astype(np.int32), out=op_f3_pfx[1:])
    im_f2_ok = (op_f3_pfx[f2_hi_all] - op_f3_pfx[f2_lo_all]) > 0

    f1c_global = np.where(im_has_b1c & im_f2_ok)[0].astype(np.int32)
    if f1c_global.size == 0:
        return

    gp_sc = np.empty(0, dtype=np.float32)
    gp_f1c = np.empty(0, dtype=np.int32)
    gp_b1c = np.empty(0, dtype=np.int32)

    stride = 400
    max_pairs = 10_000_000
    max_quad = 10_000_000
    max_quint = 10_000_000
    max_sext = 10_000_000

    n_windows = 0
    n_capped = 0
    n_pairs = 0
    im_pos_valid = im_pos[f1c_global]
    pos_start = int(im_pos_valid[0])
    pos_end = int(im_pos_valid[-1])

    for _wpos in range(pos_start, pos_end + 1, stride):
        w_lo = int(np.searchsorted(im_pos_valid, _wpos))
        w_hi = int(np.searchsorted(im_pos_valid, _wpos + stride))
        if w_lo >= w_hi:
            continue
        win_f1c = f1c_global[w_lo:w_hi]

        f1cv_loc, b1cv = _expand(b1c_lo_all[win_f1c], b1c_hi_all[win_f1c])
        if f1cv_loc.size == 0:
            continue
        f1cv = win_f1c[f1cv_loc]
        keep = (ip_pos[b1cv] - im_pos[f1cv] + ip_len[b1cv]) <= 100
        keep &= np.abs(ip_tm[b1cv] - im_tm[f1cv]) <= 1.5
        f1cv = f1cv[keep]
        b1cv = b1cv[keep]
        if f1cv.size == 0:
            continue

        pair_score = (
            -2.0 * np.abs(ip_tm[b1cv] - im_tm[f1cv])
            + (im_gcb[f1cv] + ip_gcb[b1cv]) * np.float32(0.8)
            + (im_score[f1cv] + ip_score[b1cv]) * np.float32(0.75)
        )

        n_windows += 1
        n_pairs += int(pair_score.size)
        gp_sc = np.concatenate([gp_sc, pair_score])
        gp_f1c = np.concatenate([gp_f1c, f1cv])
        gp_b1c = np.concatenate([gp_b1c, b1cv])

        if gp_sc.size > max_pairs:
            n_capped += 1
            top_k = np.argpartition(gp_sc, -max_pairs)[-max_pairs:]
            gp_sc = gp_sc[top_k]
            gp_f1c = gp_f1c[top_k]
            gp_b1c = gp_b1c[top_k]

    if gp_sc.size == 0:
        return

    print(
        f"  Phase 1: pair cap triggered {n_capped:,} times,"
        f"{gp_sc.size:,} pairs retained"
    )

    f1cv = gp_f1c
    b1cv = gp_b1c
    ia_m3_p = (im_tm[f1cv] + ip_tm[b1cv]) * 0.5 - 5.0

    # Add f2 and b2 first
    f2_lo_p = np.searchsorted(op_pos, im_pos[f1cv] - 60)
    f2_hi_p = np.searchsorted(op_pos, im_pos[f1cv] - 39)
    pairv, f2v = _expand(f2_lo_p, f2_hi_p)
    if pairv.size == 0:
        return
    f2_end_v = op_pos[f2v] + op_len[f2v]
    keep = np.abs(op_tm[f2v] - ia_m3_p[pairv]) <= 2.0
    keep &= op_has_f3[f2v]
    if use_loop:
        keep &= (im_pos[f1cv[pairv]] - f2_end_v) >= 20
    pairv = pairv[keep]
    f2v = f2v[keep]
    if pairv.size == 0:
        return
    f1ct = f1cv[pairv]
    b1ct = b1cv[pairv]
    ia_m3_t = ia_m3_p[pairv]

    b1c_end_t = ip_end[b1ct]
    f2_pos_t = op_pos[f2v]
    b2_lo_end = np.maximum(b1c_end_t + 40, f2_pos_t + 120)
    b2_hi_end = np.minimum(b1c_end_t + 61, f2_pos_t + 181)
    valid_t = b2_lo_end <= b2_hi_end
    b2_lo_s = np.where(valid_t, np.searchsorted(om_end_sorted, b2_lo_end), 0)
    b2_hi_s = np.where(valid_t, np.searchsorted(om_end_sorted, b2_hi_end + 1), 0)
    triplev, b2sv = _expand(b2_lo_s, b2_hi_s)
    if triplev.size == 0:
        return
    b2v = om_end_ord[b2sv]
    ia_m3_q = ia_m3_t[triplev]
    keep = np.abs(om_tm[b2v] - ia_m3_q) <= 2.0
    keep &= om_has_b3[b2v]
    triplev = triplev[keep]
    b2v = b2v[keep]
    ia_m3_q = ia_m3_q[keep]
    if triplev.size == 0:
        return
    f1cq = f1ct[triplev]
    b1cq = b1ct[triplev]
    f2q = f2v[triplev]

    quad_score = (
        -2.0 * np.abs(ip_tm[b1cq] - im_tm[f1cq])
        + (im_gcb[f1cq] + ip_gcb[b1cq]) * np.float32(0.8)
        + (im_score[f1cq] + ip_score[b1cq]) * np.float32(0.75)
        - 2.0 * (np.abs(op_tm[f2q] - ia_m3_q) + np.abs(om_tm[b2v] - ia_m3_q))
        + (op_gcb[f2q] + om_gcb[b2v]) * np.float32(0.15)
        + (op_score[f2q] + om_score[b2v]) * np.float32(0.3)
    )
    if f1cq.size > max_quad:  # if too many sets
        top_k = np.argpartition(quad_score, -max_quad)[-max_quad:]
        f1cq = f1cq[top_k]
        b1cq = b1cq[top_k]
        f2q = f2q[top_k]
        b2v = b2v[top_k]
        ia_m3_q = ia_m3_q[top_k]

    f3_lo_q = np.searchsorted(op_end_sorted, op_pos[f2q] - 19)
    f3_hi_q = np.searchsorted(op_end_sorted, op_pos[f2q] + 1)
    qv, f3sv = _expand(f3_lo_q, f3_hi_q)
    if qv.size == 0:
        return
    ia_m3_5 = ia_m3_q[qv]
    keep = np.abs(op_tm_sorted[f3sv] - ia_m3_5) <= 2.0
    qv = qv[keep]
    f3sv = f3sv[keep]
    ia_m3_5 = ia_m3_5[keep]
    if qv.size == 0:
        return
    if qv.size > max_quint:
        top_k = np.argpartition(np.abs(op_tm_sorted[f3sv] - ia_m3_5), max_quint)[
            :max_quint
        ]
        qv = qv[top_k]
        f3sv = f3sv[top_k]
        ia_m3_5 = ia_m3_5[top_k]
    f3v = op_end_ord[f3sv]
    b2v_5 = b2v[qv]

    b3_lo_5 = np.searchsorted(om_pos_sorted, om_end[b2v_5])
    b3_hi_5 = np.searchsorted(om_pos_sorted, om_end[b2v_5] + 20)
    p5v, b3sv = _expand(b3_lo_5, b3_hi_5)
    if p5v.size == 0:
        return
    b3v_raw = om_pos_ord[b3sv]
    keep = np.abs(om_tm[b3v_raw] - ia_m3_5[p5v]) <= 2.0
    p5v = p5v[keep]
    b3v = b3v_raw[keep]
    if p5v.size == 0:
        return

    quad_idx = qv[p5v]
    f1c6 = f1cq[quad_idx]
    b1c6 = b1cq[quad_idx]
    f2_6 = f2q[quad_idx]
    b2_6 = b2v_5[p5v]
    f3_6 = f3v[p5v]
    ia_m3_f = ia_m3_5[p5v]

    amp_s6 = op_pos[f3_6]
    amp_e6 = om_end[b3v]
    gc_num6 = gc_pfx[amp_e6] - gc_pfx[amp_s6]
    amp_ln6 = amp_e6 - amp_s6
    safe6 = np.where(amp_ln6 > 0, amp_ln6, 1)
    gc_pct6 = np.where(amp_ln6 > 0, 100.0 * gc_num6 / safe6, 50.0)
    gc_code6 = np.where(
        gc_pct6 >= 60, np.int8(1), np.where(gc_pct6 <= 45, np.int8(2), np.int8(4))
    )
    d6 = np.maximum(0.0, np.maximum(50.0 - gc_pct6, gc_pct6 - 60.0))
    gc_bonus6 = (1.0 / (1.0 + (d6 / 10.0) ** 2)).astype(np.float32)

    scores6 = (
        -2.0 * np.abs(ip_tm[b1c6] - im_tm[f1c6])
        + (im_gcb[f1c6] + ip_gcb[b1c6]) * np.float32(0.8)
        + (im_score[f1c6] + ip_score[b1c6]) * np.float32(0.75)
        - 2.0 * (np.abs(op_tm[f2_6] - ia_m3_f) + np.abs(om_tm[b2_6] - ia_m3_f))
        + (op_gcb[f2_6] + om_gcb[b2_6]) * np.float32(0.15)
        + (op_score[f2_6] + om_score[b2_6]) * np.float32(0.3)
        - (np.abs(op_tm[f3_6] - ia_m3_f) + np.abs(om_tm[b3v] - ia_m3_f))
        + (op_gcb[f3_6] + om_gcb[b3v]) * np.float32(0.15)
        + (op_score[f3_6] + om_score[b3v]) * np.float32(0.12)
        + gc_bonus6 * np.float32(0.1)
    )

    if _src_registry:  # Ensures it comes from one viable sequence
        combined_mask = (
            im_src[f1c6]
            & ip_src[b1c6]
            & op_src[f2_6]
            & om_src[b2_6]
            & op_src[f3_6]
            & om_src[b3v]
        )
        same_src = np.fromiter(
            (m != 0 for m in combined_mask), dtype=bool, count=len(combined_mask)
        )
        if not same_src.all():
            scores6 = scores6[same_src]
            f1c6 = f1c6[same_src]
            b1c6 = b1c6[same_src]
            f2_6 = f2_6[same_src]
            b2_6 = b2_6[same_src]
            f3_6 = f3_6[same_src]
            b3v = b3v[same_src]
            gc_code6 = gc_code6[same_src]
        if scores6.size == 0:
            return

    if scores6.size > max_sext:
        top_k = np.argpartition(scores6, -max_sext)[-max_sext:]
        scores6 = scores6[top_k]
        f1c6 = f1c6[top_k]
        b1c6 = b1c6[top_k]
        f2_6 = f2_6[top_k]
        b2_6 = b2_6[top_k]
        f3_6 = f3_6[top_k]
        b3v = b3v[top_k]
        gc_code6 = gc_code6[top_k]

    if lf_pos_arr.size > 0:
        lf_count6 = np.searchsorted(lf_pos_arr, im_pos[f1c6]) - np.searchsorted(
            lf_pos_arr, op_pos[f2_6] + op_len[f2_6]
        )
    else:
        lf_count6 = np.zeros(f1c6.size, dtype=np.int32)
    if lb_pos_arr.size > 0:
        lb_count6 = np.searchsorted(lb_pos_arr, om_pos[b2_6]) - np.searchsorted(
            lb_pos_arr, ip_end[b1c6]
        )
    else:
        lb_count6 = np.zeros(f1c6.size, dtype=np.int32)
    has_loop6 = ((lf_count6 > 0) | (lb_count6 > 0)).astype(np.int8)

    print(f"  Phase 2: {scores6.size:,} sextets retained")

    if use_loop:
        order = np.lexsort((-scores6, -has_loop6))
    else:
        order = np.argsort(scores6)[::-1]

    for k in order:
        yield (
            outer_plus[int(f3_6[k])],
            outer_plus[int(f2_6[k])],
            inner_minus[int(f1c6[k])],
            inner_plus[int(b1c6[k])],
            outer_minus[int(b2_6[k])],
            outer_minus[int(b3v[k])],
            int(gc_code6[k]),
            float(scores6[k]),
            bool(has_loop6[k]),
        )


def _has_loop_candidates(f2, f1c, b1c, b2, lf_pool, lf_pos, lb_pool, lb_pos, inner_avg):
    f2_end = f2.pos + f2.length
    lo = bisect.bisect_left(lf_pos, f2_end)
    hi = bisect.bisect_left(lf_pos, f1c.pos)
    for i in range(lo, hi):
        c = lf_pool[i]
        if c.pos + c.length <= f1c.pos and abs(c.Tm - inner_avg) <= 1.0:
            return True
    b1c_end = b1c.pos + b1c.length
    lo = bisect.bisect_left(lb_pos, b1c_end)
    hi = bisect.bisect_left(lb_pos, b2.pos)
    for i in range(lo, hi):
        c = lb_pool[i]
        if c.pos + c.length <= b2.pos and abs(c.Tm - inner_avg) <= 1.0:
            return True
    return False


def _find_loops(
    f2, f1c, b1c, b2, lf_pool, lf_pos, lb_pool, lb_pos, inner_avg, gc, check, core_cands
):
    f2_end = f2.pos + f2.length
    b1c_end = b1c.pos + b1c.length
    lfs = [
        c
        for c in lf_pool[
            bisect.bisect_left(lf_pos, f2_end) : bisect.bisect_left(lf_pos, f1c.pos)
        ]
        if c.pos + c.length <= f1c.pos and abs(c.Tm - inner_avg) <= 1.0
    ]
    lbs = [
        c
        for c in lb_pool[
            bisect.bisect_left(lb_pos, b1c_end) : bisect.bisect_left(lb_pos, b2.pos)
        ]
        if c.pos + c.length <= b2.pos and abs(c.Tm - inner_avg) <= 1.0
    ]
    # Rank to pick best loop candidates
    rank = lambda c: 2.0 * abs(c.Tm - inner_avg) - c.score - 4.0 * c.cov
    lfs.sort(key=rank)
    lbs.sort(key=rank)
    for lf in lfs:
        for lb in lbs:
            if not check or check_structure_loop(core_cands, lf, lb, gc):
                return lf, lb, 1.0
    for lf in lfs:
        if not check or check_structure_loop(core_cands, lf, None, gc):
            return lf, None, 0.0
    for lb in lbs:
        if not check or check_structure_loop(core_cands, None, lb, gc):
            return None, lb, 0.0
    return None, None, -50.0


def _positions_ok(f3, f2, f1c, b1c, b2, b3, lf=None, lb=None):
    f2_end = f2.pos + f2.length
    f1c_end = f1c.pos + f1c.length
    b1c_end = b1c.pos + b1c.length
    b2_end = b2.pos + b2.length

    if not (f3.pos + f3.length) <= f2.pos:
        return False
    if lf is not None:
        if not f2_end <= lf.pos:
            return False
        if not (lf.pos + lf.length) <= f1c.pos:
            return False
    if not f1c_end < b1c.pos:
        return False
    if lb is not None:
        if not b1c_end <= lb.pos:
            return False
        if not (lb.pos + lb.length) <= b2.pos:
            return False
    if not b2_end <= b3.pos:
        return False

    if not (120 <= b2_end - f2.pos <= 180):
        return False
    if not (40 <= f1c.pos - f2.pos <= 60):
        return False
    if not (40 <= b2_end - b1c_end <= 60):
        return False

    return True


def _search_seq(cand):
    return cand.seq if cand.rev == 1 else reverse_complement(cand.seq)


def _orient_mismatch(cand, positions, target_bases):
    if cand.rev == 1:
        return positions, target_bases
    L = cand.length
    return [L - 1 - p for p in positions], [
        _watson_crick.get(b, b) for b in target_bases
    ]


def _set_seq_key(r):
    roles = ("F3", "F2", "F1c", "B1c", "B2", "B3", "LF", "LB")
    return tuple(r[role].seq if r.get(role) is not None else None for role in roles)


def _sets_multiplex_compatible(a, b):
    """False if any primer interactions"""
    roles = ("F3", "F2", "F1c", "B1c", "B2", "B3", "LF", "LB")
    cands_a = [a[r] for r in roles if a.get(r) is not None]
    cands_b = [b[r] for r in roles if b.get(r) is not None]
    threshold = 45 if (a["gc"] == 2 or b["gc"] == 2) else 49
    for ca in cands_a:
        for cb in cands_b:
            if not _cached_check(
                (ca.seq, cb.seq),
                lambda ca=ca, cb=cb: _pair_ok(
                    ca.seq, cb.seq, threshold, ca.rc_seq, cb.rc_seq
                ),
            ):
                return False
    return True


def _build_one_multiplex(pool, all_targets, used, used_seqs, set_cap):
    remaining = [
        r
        for r in pool
        if id(r) not in used
        and _set_seq_key(r) not in used_seqs
        and r["rank_score"] > 0
    ]
    if not remaining:
        return []
    remaining.sort(key=lambda r: r["rank_score"], reverse=True)
    best = remaining.pop(0)
    selected = [best]
    group_seqs = {_set_seq_key(best)}
    covered = all_targets - set(best.get("failing_targets", []))
    while len(selected) < set_cap and remaining and covered != all_targets:

        def _gain(r):
            return len(all_targets - set(r.get("failing_targets", [])) - covered)

        remaining.sort(key=lambda r: (_gain(r), r["rank_score"]), reverse=True)
        picked = None
        for cand in remaining:
            if _gain(cand) <= 0:
                break
            if _set_seq_key(cand) in group_seqs:
                continue
            if all(_sets_multiplex_compatible(cand, s) for s in selected):
                picked = cand
                break
        if picked is None:
            break
        remaining.remove(picked)
        group_seqs.add(_set_seq_key(picked))
        covered |= all_targets - set(picked.get("failing_targets", []))
        selected.append(picked)
    return selected


def _multiplex_select(pool, all_targets, num_bundles, set_cap):
    if not pool or not all_targets:
        return []
    bundles = []
    used = set()
    used_seqs = set()
    for i in range(num_bundles):
        group = _build_one_multiplex(pool, all_targets, used, used_seqs, set_cap)
        if not group:
            break
        used.update(id(r) for r in group)
        used_seqs.update(_set_seq_key(r) for r in group)
        failing_sets = [set(r.get("failing_targets", [])) for r in group]
        combined_failing = (
            set.intersection(*failing_sets) if failing_sets else all_targets
        )
        cov = 100.0 * (len(all_targets) - len(combined_failing)) / len(all_targets)
        for r in group:
            r["multiplex_coverage"] = cov
        print(f"  multiplex {i + 1}: {len(group)} set(s), combined coverage={cov:.1f}%")
        bundles.append(group)
    return bundles


def assemble(
    inner,
    outer,
    loop_cands,
    seq,
    expect,
    check,
    use_loop,
    segments=None,
    targets=None,
    max_repairs=0,
    min_coverage=0.0,
    use_gblock=True,
    multiplex=0,
    seq_names=None,
    tm_params=None,
    col_maps=None,
):
    segments = segments or [(0, len(seq))]
    seg_starts = [s for s, _ in segments]

    def _segment_bounds(pos):
        i = bisect.bisect_right(seg_starts, pos) - 1
        return segments[max(0, i)]

    def _local(pos):
        i = max(0, bisect.bisect_right(seg_starts, pos) - 1)
        name = seq_names[i] if seq_names else None
        return pos - seg_starts[i], name

    t0 = time.time()
    n_checked = 0
    checks_last_set = 0
    gb_reject = {
        "padding": 0,
        "overall_gc": 0,
        "wide_window": 0,
        "tight_window": 0,
        "spread": 0,
        "homopolymer": 0,
    }
    reject_counts = {
        "core_duplicate": 0,
        "f1c_spacing": 0,
        "b1c_spacing": 0,
        "structure_fail": 0,
        "positions_invalid": 0,
    }
    gblock_cache = {}
    gb_unseen = object()
    t_last_set = time.time()
    accepted_f1c = []
    accepted_b1c = []
    accepted_core = set()
    accepted_seqs = set()
    results = []

    lf_pool = [c for c in loop_cands if c.rev == -1]
    lb_pool = [c for c in loop_cands if c.rev == 1]
    lf_pos = [c.pos for c in lf_pool]
    lb_pos = [c.pos for c in lb_pool]

    targets_idx = {name: True for name in targets} if targets else {}
    real_to_col_maps = (
        {name: _real_to_col(cm) for name, cm in col_maps.items()} if col_maps else {}
    )
    targets_arr = (
        {
            name: np.frombuffer(s.upper().encode(), dtype=np.uint8)
            for name, s in targets.items()
        }
        if targets
        else {}
    )
    primer_cache = {}

    def _mismatches(cand, name):
        key = (id(cand), name)
        cached = primer_cache.get(key)
        if cached is not None:
            return cached
        primer = _search_seq(cand)
        m = len(primer)
        local_pos, src_name = _local(cand.pos)
        start_col = real_to_col_maps[src_name][local_pos]
        end_col = real_to_col_maps[src_name][local_pos + m - 1] + 1
        target_col_map = col_maps[name]
        real_positions = [
            target_col_map[c]
            for c in range(start_col, end_col)
            if target_col_map[c] is not None
        ]
        if len(real_positions) != m or real_positions[-1] - real_positions[0] != m - 1:
            result = (m, None)
        else:
            best_off = real_positions[0]
            primer_arr = np.frombuffer(primer.encode(), dtype=np.uint8)
            n_mismatch = int(
                (targets_arr[name][best_off : best_off + m] != primer_arr).sum()
            )
            result = (n_mismatch, best_off)
        primer_cache[key] = result
        return result

    def _set_coverage(r):
        roles = [r["F3"], r["F2"], r["F1c"], r["B1c"], r["B2"], r["B3"]]
        if use_loop:
            if r["LF"] is not None:
                roles.append(r["LF"])
            if r["LB"] is not None:
                roles.append(r["LB"])

        names = list(targets_idx.keys())
        n = len(names)
        totals = {name: sum(_mismatches(c, name)[0] for c in roles) for name in names}
        all_ok = {name: totals[name] == 0 for name in names}
        strict_cov = 100.0 * sum(all_ok.values()) / n
        graded_cov = 100.0 * sum(_target_credit(totals[name]) for name in names) / n
        failing_targets = [name for name, ok in all_ok.items() if not ok]
        return strict_cov, graded_cov, totals, failing_targets

    def _set_mismatches(r, failing_targets):
        role_labels = list(_CORE_ROLES)
        if use_loop:
            if r["LF"] is not None:
                role_labels.append("LF")
            if r["LB"] is not None:
                role_labels.append("LB")

        mismatch_data = {}
        for label in role_labels:
            cand = r[label]
            primer = _search_seq(cand)
            m = len(primer)
            role_fail = {}
            for name in failing_targets:
                _, best_off = _mismatches(cand, name)
                if best_off is None:
                    role_fail[name] = repairSet.MismatchInfo(None, [])
                    continue
                target_seq = targets[name].upper()
                positions = [
                    i for i in range(m) if primer[i] != target_seq[best_off + i]
                ]
                target_bases = [target_seq[best_off + i] for i in positions]
                role_fail[name] = repairSet.MismatchInfo(positions, target_bases)
            if role_fail:
                mismatch_data[label] = role_fail
        return mismatch_data

    def _score_coverage(r):
        if not targets_idx:
            return
        strict_cov, graded_cov, totals, failing = _set_coverage(r)
        r["coverage"] = strict_cov
        r["mismatch_totals"] = totals
        if failing:
            r["failing_targets"] = failing
            r["mismatch_data"] = _set_mismatches(r, failing)
        else:
            r.pop("failing_targets", None)
            r.pop("mismatch_data", None)
        r["rank_score"] = r["score"] - _COV_SCALE * (100.0 - graded_cov)
        if failing:
            r["rank_score"] -= _FAIL_REPAIR

    def _commit(r):
        nonlocal t_last_set, checks_last_set
        seq_key = _set_seq_key(r)
        if seq_key in accepted_seqs:
            return
        accepted_seqs.add(seq_key)
        r["turn"] = len(results) + 1
        n_loops = (r["LF"] is not None) + (r["LB"] is not None)
        loop_tag = " [loop]" if n_loops > 0 else ""
        since_last = time.time() - t_last_set
        new_checks = n_checked - checks_last_set
        cov_tag = f" coverage={r['coverage']:.1f}%" if targets_idx else ""
        failing = r.get("failing_targets")
        uncov_tag = f"  ({len(failing)} target(s) uncovered)" if failing else ""
        ms_tag = f" [multiplex {r['multiplex']}]" if r.get("multiplex") else ""
        degen_tag = f" [degen x{r['degenerate']}]" if r.get("degenerate") else ""
        local_f3, src_name = _local(r["F3"].pos)
        src_tag = f" [{src_name}]" if src_name else ""
        print(
            f"  [Set {r['turn']}]{ms_tag}{degen_tag}{src_tag} score={r['rank_score']:.2f} F3@{local_f3}"
            f" GC={r['amp_gc']:.1f}%{loop_tag}{cov_tag}{uncov_tag}  (+{since_last:.2f}s, {new_checks:,} checks)"
        )
        t_last_set = time.time()
        checks_last_set = n_checked
        results.append(r)

    pending_cap = max(20_000, expect)
    pending_heap = []
    push_order = 0
    score_floor = float("-inf")
    stale_limit = 30_000_000
    stale = 0
    scan_time_limit = 600
    scan_start = time.time()
    n_scanned = 0

    def _offer(r, core, f1c_pos, b1c_pos):
        nonlocal score_floor, stale, push_order
        push_order += 1
        if len(pending_heap) < pending_cap:
            accepted_core.add(core)
            bisect.insort(accepted_f1c, f1c_pos)
            bisect.insort(accepted_b1c, b1c_pos)
            heapq.heappush(
                pending_heap, (r["score"], push_order, r, core, f1c_pos, b1c_pos)
            )
            stale = 0
            if len(pending_heap) == pending_cap:
                score_floor = pending_heap[0][0]
            return
        if r["score"] <= pending_heap[0][0]:
            return
        _, _, _evicted, ev_core, ev_f1c, ev_b1c = heapq.heapreplace(
            pending_heap, (r["score"], push_order, r, core, f1c_pos, b1c_pos)
        )
        accepted_core.discard(ev_core)
        accepted_f1c.pop(bisect.bisect_left(accepted_f1c, ev_f1c))
        accepted_b1c.pop(bisect.bisect_left(accepted_b1c, ev_b1c))
        accepted_core.add(core)
        bisect.insort(accepted_f1c, f1c_pos)
        bisect.insort(accepted_b1c, b1c_pos)
        score_floor = pending_heap[0][0]
        stale = 0

    group_has_loop = None
    group_done = False
    for f3, f2, f1c, b1c, b2, b3, gc, combo_score, has_loop in _combo_chunk_gen(
        inner, outer, seq, use_loop, loop_cands
    ):
        stale += 1
        n_scanned += 1
        if n_scanned % 200_000 == 0 and time.time() - scan_start > scan_time_limit:
            print(
                f"  Scan time budget ({scan_time_limit}s) reached after "
                f"stopping early with {len(pending_heap)} candidate set(s)"
            )
            break
        if use_loop:
            if has_loop != group_has_loop:
                group_has_loop = has_loop
                group_done = False
                stale = 0
            if group_done:
                continue
            if len(pending_heap) >= pending_cap:
                bound = combo_score + (5.5 if has_loop else 0.0)
                if bound <= score_floor:
                    group_done = True
                    continue
            elif stale > stale_limit:
                group_done = True
                continue
        elif len(pending_heap) >= pending_cap and combo_score <= score_floor:
            break
        elif stale > stale_limit:
            print(
                f"  Stale limit ({stale_limit:,}) reached after "
                f"{n_scanned:,} combos scanned -- stopping scan"
            )
            break
        core = (f1c.pos, b1c.pos)
        if core in accepted_core:
            reject_counts["core_duplicate"] += 1
            continue
        if not check_add(f1c.pos, accepted_f1c, len(seq)):
            reject_counts["f1c_spacing"] += 1
            continue
        if not check_add(b1c.pos, accepted_b1c, len(seq)):
            reject_counts["b1c_spacing"] += 1
            continue
        inner_avg = (f1c.Tm + b1c.Tm) * 0.5 if use_loop else 0.0
        loop_possible = (
            use_loop
            and has_loop
            and _has_loop_candidates(
                f2, f1c, b1c, b2, lf_pool, lf_pos, lb_pool, lb_pos, inner_avg
            )
        )
        core_cands = [f3, f2, f1c, b1c, b2, b3]
        gb_start = gb_end = None
        if use_gblock:
            gb_key = (f3.pos, b3.pos, b3.length)
            cached = gblock_cache.get(gb_key, gb_unseen)
            if cached is gb_unseen:
                seg_lo, seg_hi = _segment_bounds(f3.pos)
                gb_opt = optimize_gblock(
                    seq,
                    f3.pos,
                    b3.pos + b3.length,
                    lo_limit=seg_lo,
                    hi_limit=seg_hi,
                )
                if gb_opt is None:
                    cached = ("padding", None, None)
                else:
                    gs, ge = gb_opt
                    if _region_gc(seq, gs, ge) > DEFAULT_GBLOCK["gc_max"]:
                        cached = ("overall_gc", None, None)
                    else:
                        gb_reason = _gblock_fail_reason(seq, gs, ge)
                        cached = (
                            (gb_reason, gs, ge)
                            if gb_reason is None
                            else (gb_reason, None, None)
                        )
                gblock_cache[gb_key] = cached
            gb_reason, gb_start, gb_end = cached
            if gb_reason is not None:
                gb_reject[gb_reason] += 1
                continue
        if check:
            n_checked += 1
            if not check_structure_loop(core_cands, None, None, gc):
                reject_counts["structure_fail"] += 1
                continue
        lf_cand = lb_cand = None
        loop_bonus = 0.0
        if loop_possible:
            lf_cand, lb_cand, loop_bonus = _find_loops(
                f2,
                f1c,
                b1c,
                b2,
                lf_pool,
                lf_pos,
                lb_pool,
                lb_pos,
                inner_avg,
                gc,
                check,
                core_cands,
            )
        if not _positions_ok(f3, f2, f1c, b1c, b2, b3, lf_cand, lb_cand):
            reject_counts["positions_invalid"] += 1
            continue
        amp_start = f3.pos
        amp_end = b3.pos + b3.length
        amp_sub = seq[amp_start:amp_end].upper()
        amp_gc = 100.0 * (amp_sub.count("G") + amp_sub.count("C")) / len(amp_sub)
        loop_cand_score = (lf_cand.score if lf_cand else 0.0) + (
            lb_cand.score if lb_cand else 0.0
        )
        loop_weight = 5
        final_score = combo_score + loop_weight * loop_bonus + 0.1 * loop_cand_score
        if use_loop and lf_cand is None and lb_cand is None:
            final_score -= 10.0
        r = {
            "score": final_score,
            "rank_score": final_score,
            "gc": gc,
            "amp_gc": amp_gc,
            "F3": f3,
            "F2": f2,
            "F1c": f1c,
            "B1c": b1c,
            "B2": b2,
            "B3": b3,
            "LF": lf_cand,
            "LB": lb_cand,
        }
        if gb_start is not None:
            r["gb_start"] = gb_start
            r["gb_end"] = gb_end
        _offer(r, core, f1c.pos, b1c.pos)

    pending = [item[2] for item in pending_heap]
    # print(f"  Pool after scan: {len(pending):,} distinct cores retained")
    # print(f"  Rejects: {reject_counts}")
    # print(f"  gBlock rejects: {gb_reject}")

    def _degen_penalty(n):
        return 0.3 * n

    def _to_seq_orientation(cand, m):
        if m.pos is None:
            return repairSet.MismatchInfo(None, m.target_bases)
        positions, target_bases = _orient_mismatch(cand, m.pos, m.target_bases)
        return repairSet.MismatchInfo(positions, target_bases)

    def _run_degenerate_repair(r):
        mismatch_data = r["mismatch_data"]
        primers = {}
        hits = {}
        role_mismatches_seq = {}
        for role, per_target in mismatch_data.items():
            cand = r[role]
            primers[role] = repairSet.PrimerInfo(
                role=role,
                pos=cand.pos,
                length=cand.length,
                rev=cand.rev,
                seq=cand.seq,
                Tm=cand.Tm,
                GC=cand.GC,
            )
            remapped = {
                name: _to_seq_orientation(cand, m) for name, m in per_target.items()
            }
            role_mismatches_seq[role] = remapped
            hits[role] = {
                name: (m.pos is not None and not m.pos) for name, m in remapped.items()
            }

        threshold = 49 if r["gc"] in (1, 4) else 45
        config = repairSet.RepairConfig(
            max_mismatches=0,
            coverage_goal=100.0,
            repair_budget=max_repairs,
            tm_penalty_scale=0.75,
            compat_threshold=threshold,
        )
        return repairSet.repair_set(
            primers, hits, role_mismatches_seq, tm_params, config
        )

    def _attempt_repair(r):
        result = _run_degenerate_repair(r)
        n_degen = sum(
            sum(1 for c in result.primers[role].seq if c not in "ACGT")
            for role in result.modified_roles
        )
        repaired = bool(result.modified_roles) and n_degen <= max_repairs
        new_failing = sorted(
            result.unresolved_targets if repaired else r["failing_targets"]
        )
        n_total = len(targets_idx)
        new_coverage = 100.0 * (n_total - len(new_failing)) / n_total
        new_score = r["score"]
        if repaired:
            new_score -= _degen_penalty(n_degen) + result.thermo_penalty
        new_rank_score = new_score - _COV_SCALE * (100.0 - new_coverage)
        if new_failing:
            new_rank_score -= _FAIL_REPAIR
        if repaired:
            for role in result.modified_roles:
                r[role] = replace(r[role], seq=result.primers[role].seq)
            r["degenerate"] = n_degen
        r["coverage"] = new_coverage
        if new_failing:
            r["failing_targets"] = new_failing
            new_failing_set = set(new_failing)
            filtered = {}
            for role, per_target in r["mismatch_data"].items():
                role_hits = result.hits.get(role, {}) if repaired else {}
                kept = {
                    name: m
                    for name, m in per_target.items()
                    if name in new_failing_set and not role_hits.get(name, False)
                }
                if kept:
                    filtered[role] = kept
            r["mismatch_data"] = filtered
        else:
            r.pop("failing_targets", None)
            r.pop("mismatch_data", None)
        r["score"] = new_score
        r["rank_score"] = new_rank_score

    repair_limit = 200
    repair_count = 0

    def _finalize_one(r):
        """Score coverage for pending set. Attempt repair if on as well."""
        nonlocal repair_count
        _score_coverage(r)
        if (
            targets_idx
            and tm_params is not None
            and max_repairs > 0
            and r.get("failing_targets")
        ):
            if repair_count < repair_limit:
                repair_count += 1
                _attempt_repair(r)

    if multiplex and targets_idx:
        for r in pending:
            _finalize_one(r)
        if min_coverage > 0:
            pending = [r for r in pending if r.get("coverage", 0.0) >= min_coverage]
        bundles = _multiplex_select(pending, set(targets_idx.keys()), expect, multiplex)
        for group_i, group in enumerate(bundles, 1):
            for r in group:
                r["multiplex"] = group_i
                _commit(r)
    else:
        pending.sort(key=lambda p: p["score"], reverse=True)
        n_pending = len(pending)
        idx = 0
        rank_heap = []
        final_core = set()
        final_f1c = []
        final_b1c = []
        while len(results) < expect:
            while idx < n_pending and (
                not rank_heap or -rank_heap[0][0] <= pending[idx]["score"]
            ):
                r = pending[idx]
                idx += 1
                _finalize_one(r)
                if (
                    targets_idx
                    and min_coverage > 0
                    and r.get("coverage", 0.0) < min_coverage
                ):
                    continue
                heapq.heappush(rank_heap, (-r["rank_score"], id(r), r))
            if not rank_heap:
                break
            _, _, r = heapq.heappop(rank_heap)
            core = (r["F1c"].pos, r["B1c"].pos)
            if core in final_core:
                continue
            if not check_add(r["F1c"].pos, final_f1c, len(seq)):
                continue
            if not check_add(r["B1c"].pos, final_b1c, len(seq)):
                continue
            final_core.add(core)
            bisect.insort(final_f1c, r["F1c"].pos)
            bisect.insort(final_b1c, r["B1c"].pos)
            _commit(r)

    if multiplex and targets_idx:
        for i, r in enumerate(results, 1):
            r["turn"] = i
    else:
        results.sort(key=lambda r: r["rank_score"], reverse=True)
        for i, r in enumerate(results, 1):
            r["turn"] = i

    print(f"Sets found in {time.time()-t0:.2f}s")

    return results


def write_results(
    results,
    fp,
    use_loop,
    seq=None,
    segments=None,
    seq_names=None,
):
    seg_starts = [s for s, _ in segments] if segments else None

    def _local(pos):
        """Change to be position within individual sequence instead of overall"""
        if not seg_starts:
            return pos, None
        i = max(0, bisect.bisect_right(seg_starts, pos) - 1)
        name = seq_names[i] if seq_names else None
        return pos - seg_starts[i], name

    coverages = [r["coverage"] for r in results if r.get("coverage") is not None]
    if coverages:
        fp.write(
            f"{len(results)} set(s)  avg coverage={sum(coverages) / len(coverages):.1f}%"
            f"  best coverage={max(coverages):.1f}%\n\n"
        )
    last_bundle = None
    bundle_turn = 0
    for r in results:
        ms = r.get("multiplex")
        if ms is not None:
            if ms != last_bundle:
                mcov = r.get("multiplex_coverage")
                cov_str = f" (coverage = {mcov:.1f}%)" if mcov is not None else ""
                fp.write(f"--- Multiplex {ms}{cov_str} ---\n\n")
                bundle_turn = 0
            bundle_turn += 1
            display_turn = bundle_turn
        else:
            display_turn = r["turn"]
        last_bundle = ms
        f3, b3 = r["F3"], r["B3"]
        amp_len = b3.pos + b3.length - f3.pos
        local_f3, src_name = _local(f3.pos)
        local_b3_end, _ = _local(b3.pos + b3.length)
        src_tag = f"  source:{src_name}" if src_name else ""
        cov = r.get("coverage")
        cov_tag = f"  coverage={cov:.1f}%" if cov is not None else ""
        degen_tag = (
            f"  degenerate positions:{r['degenerate']}" if r.get("degenerate") else ""
        )
        fp.write(
            f"LAMP set {display_turn}  score={r['rank_score']:.2f}"
            f"  amplicon {local_f3}-{local_b3_end} ({amp_len} bp)"
            f"  regional GC={r['amp_gc']:.1f}%{cov_tag}{degen_tag}{src_tag}\n"
        )
        for key, label in (
            ("F3", "F3 "),
            ("F2", "F2 "),
            ("F1c", "F1c"),
            ("B1c", "B1c"),
            ("B2", "B2 "),
            ("B3", "B3 "),
        ):
            c = r[key]
            local_pos, _ = _local(c.pos)
            fp.write(
                f"  {label}: pos:{local_pos}  len:{c.length} bp"
                f"  Tm:{c.Tm:.1f}C  GC:{c.GC:.1f}%"
                f"  seq(5'-3'):{c.seq}\n"
            )
        if use_loop:
            for key, label in (("LF", "LF "), ("LB", "LB ")):
                c = r[key]
                if c is None:
                    fp.write(f"  {label}: none found\n")
                else:
                    local_pos, _ = _local(c.pos)
                    fp.write(
                        f"  {label}: pos:{local_pos}  len:{c.length} bp"
                        f"  Tm:{c.Tm:.1f}C  GC:{c.GC:.1f}%"
                        f"  seq(5'-3'):{c.seq}\n"
                    )
        mismatch_data = r.get("mismatch_data")
        if mismatch_data:
            failing = r.get("failing_targets", [])
            fp.write(f"  Uncovered targets ({len(failing)}): {', '.join(failing)}\n")
            for role, per_target in mismatch_data.items():
                for name, m in per_target.items():
                    if m.pos is None:
                        fp.write(f"    {role} vs {name}: no seedable match")
                    elif m.pos:
                        positions, target_bases = _orient_mismatch(
                            r[role], m.pos, m.target_bases
                        )
                        detail = ", ".join(
                            f"{p}:{b}" for p, b in zip(positions, target_bases)
                        )
                        fp.write(
                            f"    {role} vs {name}: {len(m.pos)} mismatch(es)"
                            f" at pos [{detail}]\n"
                        )

        if seq is not None:
            gb_start = r.get("gb_start", max(0, f3.pos - 40))
            gb_end = r.get("gb_end", min(len(seq), b3.pos + b3.length + 40))
            gblock_seq = seq[gb_start:gb_end]
            local_gb_start, _ = _local(gb_start)
            local_gb_end, _ = _local(gb_end)
            amb_pos = [i for i, ch in enumerate(gblock_seq) if ch.upper() not in "ACGT"]
            amb_tag = f"  ambiguous base positions:{amb_pos}" if amb_pos else ""
            fp.write(
                f"  gblock_{display_turn} region: {local_gb_start}-{local_gb_end}({gb_end - gb_start})"
                f"  ambiguous bases:{len(amb_pos)}{amb_tag}\n"
            )
            fp.write(f"  gblock_{display_turn}: {gblock_seq}\n")
        fp.write("\n")


def parse_args():
    p = argparse.ArgumentParser(
        prog="LAMP",
        usage="LAMP -in <ref_genome> [options]*",
        add_help=False,
    )
    p.add_argument("-in", dest="input", metavar="<ref_genome>", default=None)
    p.add_argument("-num", dest="num", type=int, default=10, metavar="<int>")
    p.add_argument(
        "-check",
        dest="check",
        type=int,
        default=1,
        metavar="<int>",
        help="Secondary structure check (0=off, other=on)",
    )
    p.add_argument(
        "-loop",
        dest="loop",
        type=lambda x: x.lower() not in ("false", "0"),
        default=True,
        metavar="<True|False>",
    )
    p.add_argument(
        "-gblock",
        dest="gblock",
        type=lambda x: x.lower() not in ("false", "0"),
        default=True,
        metavar="<True|False>",
        help="Design/require a synthesizable gBlock for each set (default: True)",
    )
    p.add_argument(
        "-multiplex",
        dest="multiplex",
        type=int,
        default=0,
        metavar="<int>",
        help="Instead of single sets, allow multiplex bundles of up to -multiplex sets each",
    )
    p.add_argument("-par", dest="par", metavar="<par_directory>", default=None)
    p.add_argument("-out", dest="out", metavar="<output_file>", default=None)
    p.add_argument(
        "-max_repairs",
        dest="max_repairs",
        type=int,
        default=0,
        metavar="<int>",
        help="How many bases in a primer set that can be repaired. Defaults "
        "to no repairs (0)",
    )
    p.add_argument(
        "-min_coverage",
        dest="min_coverage",
        type=float,
        default=0.0,
        metavar="<0-100>",
        help="Sets must be above this coverage %% (default: 0)",
    )
    p.add_argument("-h", "-help", dest="help", action="store_true", default=False)
    args = p.parse_args()

    if args.help:
        print("USAGE:")
        print("  LAMP -in <ref_genome> [options]*\n")
        print("ARGUMENTS:")
        print("  -in <ref_genome>   Reference genome (FASTA format)")
        print("  -num <int>         Number of primer sets to find (default: 10)")
        print("  -check <int>       Secondary structure check; 0=off (default: 1)")
        print("  -loop <bool>       Design loop primers (default: True)")
        print(
            "  -gblock <bool>     Design/require a synthesizable gBlock (default: True)"
        )
        print("  -par <directory>   Parameter file directory (default: Par/)")
        print("  -out <file>        Output file (default: stdout)")
        print(
            "  -max_repairs <int>  How many bases in a primer set that can be repaired. Defaults to no repairs (0)"
        )
        print(
            "  -min_coverage <0-100>  Sets must be above this coverage % (default: 0)"
        )
        print(
            "  -multiplex <int>   Instead of single sets, allow multiplex bundles of up to -multiplex sets each. Default: 0"
        )
        sys.exit(0)

    if args.input is None:
        print("Needs sequence file (-in .fasta)")
        p.print_usage()
        sys.exit(1)

    par_path = args.par
    if par_path is None:
        par_path = os.path.join(os.getcwd(), "Par") + os.sep
    elif not par_path.endswith(("/", os.sep)):
        par_path += os.sep

    return {
        "input": args.input,
        "num": args.num,
        "check": bool(args.check),
        "loop": args.loop,
        "gblock": args.gblock,
        "par_path": par_path,
        "out": args.out,
        "max_repairs": args.max_repairs,
        "min_coverage": args.min_coverage,
        "multiplex": args.multiplex,
    }


def main():
    args = parse_args()
    t_total = time.time()

    t0 = time.time()
    gapped_sequences = read_fasta(args["input"], strip_gaps=False)
    if not gapped_sequences:
        print("Error: no sequences found in input FASTA", file=sys.stderr)
        sys.exit(1)
    if len(gapped_sequences) > 1:
        validate_alignment(gapped_sequences)
        sequences = {n: s.replace("-", "") for n, s in gapped_sequences.items()}
        seq = "".join(sequences.values())
        segments = []
        offset = 0
        for s in sequences.values():
            segments.append((offset, offset + len(s)))
            offset += len(s)
        targets = sequences
        col_maps = {name: _col_to_real(s) for name, s in gapped_sequences.items()}
    else:
        sequences = {n: s.replace("-", "") for n, s in gapped_sequences.items()}
        seq = next(iter(sequences.values()))
        segments = [(0, len(seq))]
        targets = None
        col_maps = None
    seq_names = list(sequences.keys()) if len(sequences) > 1 else None
    print(f"Loaded reference sequence: {len(seq):,} bp  ({time.time()-t0:.2f}s)")

    results_dir = os.path.join(os.getcwd(), "Results")

    t0 = time.time()
    inner = load_candidates(os.path.join(results_dir, "InnerCandidates"))
    print(f"Loaded InnerCandidates: {len(inner):,} candidates  ({time.time()-t0:.2f}s)")

    t0 = time.time()
    outer = load_candidates(os.path.join(results_dir, "OuterCandidates"))
    print(f"Loaded OuterCandidates: {len(outer):,} candidates  ({time.time()-t0:.2f}s)")

    if args["loop"]:
        t0 = time.time()
        loop_cands = load_candidates(os.path.join(results_dir, "LoopCandidates"))
        print(
            f"Loaded LoopCandidates: {len(loop_cands):,} candidates  ({time.time()-t0:.2f}s)"
        )
    else:
        loop_cands = []

    if not inner or not outer:
        print("No candidates found in Results/. Run single.py first.", file=sys.stderr)
        sys.exit(1)

    print(
        f"\nAssembling LAMP primer sets (check={args['check']}, loop={args['loop']})..."
    )
    t0 = time.time()
    results = assemble(
        inner,
        outer,
        loop_cands,
        seq,
        expect=args["num"],
        check=args["check"],
        use_loop=args["loop"],
        segments=segments,
        targets=targets,
        max_repairs=args["max_repairs"],
        min_coverage=args["min_coverage"],
        use_gblock=args["gblock"],
        multiplex=args["multiplex"],
        seq_names=seq_names,
        tm_params=load_params(args["par_path"]) if targets else None,
        col_maps=col_maps,
    )
    print(f"Assembly complete: {len(results)} sets found in {time.time()-t0:.2f}s")

    t0 = time.time()
    out_path = args["out"] or os.path.join(results_dir, "results.txt")
    out_fp = open(out_path, "w")
    try:
        write_results(
            results, out_fp, args["loop"], seq, segments=segments, seq_names=seq_names
        )
    finally:
        out_fp.close()
    print(f"Results written to {out_path}  ({time.time()-t0:.2f}s)")

    print(f"\nTotal runtime: {time.time()-t_total:.2f}s")

    if not results:
        print("No LAMP primer sets found.", file=sys.stderr)


if __name__ == "__main__":
    t0 = time.time()
    main()
