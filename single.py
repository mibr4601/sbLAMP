import argparse
import hashlib
import io
import json
import math
import os
import re
import sys

import numpy as np
import primer3 as _p3
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform

_P3_HOM1 = {
    "mv_conc": 30,
    "dv_conc": 4.0,
    "dntp_conc": 1.6,
    "dna_conc": 800,
    "temp_c": 65.0,
}
_P3_HOM2 = {
    "mv_conc": 30,
    "dv_conc": 4.0,
    "dntp_conc": 1.6,
    "dna_conc": 200,
    "temp_c": 65.0,
}


def read_fasta(filepath, strip_gaps=True):
    sequences = {}
    name = None
    sequence = []

    def _finalize(name, sequence):
        seq = "".join(sequence)
        if strip_gaps:
            n_gaps = seq.count("-")
            if n_gaps:
                print(
                    f"Warning: {filepath}: '{name}' contained {n_gaps} gap"
                    f" character(s) -- stripped."
                )
                seq = seq.replace("-", "")
        return seq

    with open(filepath, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    sequences[name] = _finalize(name, sequence)
                name = line[1:].split()[0]
                sequence = []
            else:
                sequence.append(line.upper())
    if name is not None:
        sequences[name] = _finalize(name, sequence)
    return sequences


def load_tm_params(par_path):
    deltah = [0.0] * 16
    deltas = [0.0] * 16
    with open(os.path.join(par_path, "santalucia_tm_nn_parameters.txt")) as f:
        for line in f:
            if line.lstrip().startswith("#"):
                continue
            parts = line.split()
            if len(parts) != 3:
                continue
            nn = parts[0].upper()
            h = float(parts[1])
            s = float(parts[2])
            rc = _watson_crick[nn[1]] + _watson_crick[nn[0]]
            i1 = _base[nn[0]] * 4 + _base[nn[1]]
            i2 = _base[rc[0]] * 4 + _base[rc[1]]
            deltah[i1] = deltah[i2] = h
            deltas[i1] = deltas[i2] = s
    return deltah, deltas


def load_stab(deltah, deltas):
    inv = "ATCG"
    t = (273.15 + 37.0) / 1000.0
    at_init = 2.3 - t * 4.1
    gc_init = 0.1 + t * 2.8
    stab = [0.0] * 4096
    for i in range(4096):
        h = [inv[(i >> (2 * k)) & 3] for k in range(5, -1, -1)]
        g = (at_init if h[0] in "AT" else gc_init) + (
            at_init if h[5] in "AT" else gc_init
        )
        for k in range(5):
            j = _base[h[k]] * 4 + _base[h[k + 1]]
            g += t * deltas[j] - deltah[j]
        stab[i] = abs(g)
    return stab


def load_params(par_path):
    deltah, deltas = load_tm_params(par_path)
    return {
        "deltah": deltah,
        "deltas": deltas,
        "stab": load_stab(deltah, deltas),
    }


_complement = str.maketrans("ACGT", "TGCA")
_base = {"A": 0, "T": 1, "C": 2, "G": 3}
_watson_crick = {"A": "T", "T": "A", "C": "G", "G": "C"}
_NON_ACGT = re.compile(r"[^ACGT]")


def reverse_complement(seq):
    return seq.translate(_complement)[::-1]


def gc_content(seq):
    return 100.0 * (seq.count("G") + seq.count("C")) / len(seq)


def calc_tm(seq, deltah, deltas):
    n = len(seq)
    total_h = -sum(deltah[_base[seq[i]] * 4 + _base[seq[i + 1]]] for i in range(n - 1))
    total_s = -sum(deltas[_base[seq[i]] * 4 + _base[seq[i + 1]]] for i in range(n - 1))
    for end in (seq[0], seq[-1]):
        if end in "AT":
            total_h += 2.3
            total_s += 4.1
        else:
            total_h += 0.1
            total_s -= 2.8
    return 1000.0 * total_h / (total_s - 0.51986 * (n - 1) - 36.70381) - 273.15


def calc_stability(seq, stab):
    pos = 0
    for b in seq[:6]:
        pos = pos * 4 + _base[b]
    s5 = stab[pos]
    pos = 0
    for b in seq[-6:]:
        pos = pos * 4 + _base[b]
    s3 = stab[pos]
    return s5, s3


def dimer_ok(primer):
    if primer[-1] == primer[-2] == primer[-3] == primer[-4]:
        return False
    if (
        _watson_crick.get(primer[-1]) == primer[-6]
        and _watson_crick.get(primer[-2]) == primer[-5]
        and _watson_crick.get(primer[-3]) == primer[-4]
    ):
        return False
    return True


def _needs_bimolecular(seq):
    n = len(seq)
    alpha_idx = {"A": 0, "C": 1, "G": 2, "T": 3}
    num1 = [alpha_idx[b] for b in seq]
    num2 = [alpha_idx[b] for b in reversed(seq)]
    for dj in range(-(n - 1), n):
        run = 0
        for i in range(n):
            j = i + dj
            if 0 <= j < n and num1[i] + num2[j] == 3:
                run += 1
                if run > 6:
                    return True
            else:
                run = 0
    return False


def _gc_clamp_bonus(seq):
    tail = seq[-6:].upper()
    run = 0
    for b in reversed(tail):
        if b in "GCS":
            run += 1
        else:
            break
    if run <= 1:
        return 0.0
    if run <= 3:
        return 3.0
    return -2.0


def _primer_score(Tm, gc, s5, s3, target_tm, seq=""):
    if 50.0 <= gc <= 60.0:
        gc_score = 3.0
    else:
        gc_score = max(0.0, 3.0 - min(abs(gc - 50.0), abs(gc - 60.0)) / 5.0)
    clamp_score = _gc_clamp_bonus(seq) if seq else 0.0
    tm_dev = abs(Tm - target_tm)
    tm_score = 1.0 if tm_dev <= 5.0 else max(0.0, 1.0 - (tm_dev - 5.0) * 0.4)
    stab_score = 0.2 * (s5 + s3)
    return gc_score + clamp_score + tm_score + stab_score


def _secondary_ok(primer, Tm):
    thresh = Tm - 10
    worst_dg = 0.0
    r = _p3.calc_hairpin(primer, **_P3_HOM1)
    if r.structure_found:
        if r.dg < -3000:
            return False, 0.0
        if r.tm > thresh:
            return False, 0.0
        worst_dg = min(worst_dg, r.dg)
    if _needs_bimolecular(primer):
        r1 = _p3.calc_homodimer(primer, **_P3_HOM1)
        if r1.structure_found:
            if r1.dg < -8000:
                return False, 0.0
            if r1.tm > thresh:
                return False, 0.0
            worst_dg = min(worst_dg, r1.dg)
        r2 = _p3.calc_end_stability(primer, primer, **_P3_HOM2)
        if r2.structure_found:
            if r2.dg < -2000:
                return False, 0.0
            if r2.tm > thresh:
                return False, 0.0
            worst_dg = min(worst_dg, r2.dg)
    return True, max(0.0, (-6000 - worst_dg) / 1000)


def _thal_p3(seq1, seq2, thal_type):
    if thal_type == 3:
        return _thal_p3(seq2, seq1, 2)
    if min(len(seq1), len(seq2)) < 4:
        return 0.0
    if thal_type == 1:
        r = _p3.calc_heterodimer(seq1, seq2, **_P3_HOM1)
    else:
        r = _p3.calc_end_stability(seq1, seq2, **_P3_HOM2)
    return r.tm if r.structure_found else 0.0


def _pair_ok(a, b, threshold, rc_a="", rc_b=""):
    if not rc_a:
        rc_a = reverse_complement(a)
    if not rc_b:
        rc_b = reverse_complement(b)
    if _thal_p3(a, b, 1) > threshold or _thal_p3(a, b, 2) > threshold:
        return False
    if _thal_p3(b, a, 2) > threshold:
        return False
    if _thal_p3(rc_b, rc_a, 2) > threshold:
        return False
    if _thal_p3(rc_a, rc_b, 2) > threshold:
        return False
    return True


DEFAULT_DESIGN = {
    "region_thresholds": {
        "at_rich_upper": 47.5,
        "gc_rich_lower": 57.5,
        "neutral_lower": 42.5,
        "neutral_upper": 62.5,
    },
    "regions": {
        "AT-rich": {
            "outer": {
                "gc_min": 30,
                "gc_max": 65,
                "tm_min": 55.0,
                "tm_max": 58.0,
                "len_min": 18,
                "len_max": 25,
            },
            "inner": {
                "gc_min": 30,
                "gc_max": 65,
                "tm_min": 60.0,
                "tm_max": 63.0,
                "len_min": 20,
                "len_max": 25,
            },
            "loop": {
                "gc_min": 30,
                "gc_max": 65,
                "tm_min": 60.0,
                "tm_max": 63.0,
                "len_min": 20,
                "len_max": 25,
            },
        },
        "GC-rich": {
            "outer": {
                "gc_min": 40,
                "gc_max": 70,
                "tm_min": 59.0,
                "tm_max": 63.0,
                "len_min": 15,
                "len_max": 20,
            },
            "inner": {
                "gc_min": 40,
                "gc_max": 70,
                "tm_min": 64.0,
                "tm_max": 68.0,
                "len_min": 15,
                "len_max": 22,
            },
            "loop": {
                "gc_min": 40,
                "gc_max": 70,
                "tm_min": 64.0,
                "tm_max": 68.0,
                "len_min": 15,
                "len_max": 22,
            },
        },
        "Neutral": {
            "outer": {
                "gc_min": 40,
                "gc_max": 65,
                "tm_min": 59.0,
                "tm_max": 61.0,
                "len_min": 18,
                "len_max": 20,
            },
            "inner": {
                "gc_min": 40,
                "gc_max": 65,
                "tm_min": 64.0,
                "tm_max": 66.0,
                "len_min": 20,
                "len_max": 22,
            },
            "loop": {
                "gc_min": 40,
                "gc_max": 65,
                "tm_min": 64.0,
                "tm_max": 66.0,
                "len_min": 20,
                "len_max": 22,
            },
        },
    },
    "global": {
        "min_len": 15,
        "max_len": 25,
        "tm_min": 55.0,
        "tm_max": 68.0,
        "gc_min": None,  # scales with sequence length by default
        "gc_max": None,
    },
}


def load_design(path):
    if not path:
        return DEFAULT_DESIGN
    with open(path) as fh:
        return {**DEFAULT_DESIGN, **json.load(fh)}


def _region_ok(gc, Tm, length, rgc, regions, ptype, rt):
    b = regions["AT-rich"][ptype]
    if (
        rgc < rt["at_rich_upper"]
        and b["gc_min"] <= gc <= b["gc_max"]
        and b["tm_min"] <= Tm <= b["tm_max"]
        and b["len_min"] <= length <= b["len_max"]
    ):
        return True
    b = regions["GC-rich"][ptype]
    if (
        rgc > rt["gc_rich_lower"]
        and b["gc_min"] <= gc <= b["gc_max"]
        and b["tm_min"] <= Tm <= b["tm_max"]
        and b["len_min"] <= length <= b["len_max"]
    ):
        return True
    b = regions["Neutral"][ptype]
    if (
        rt["neutral_lower"] <= rgc <= rt["neutral_upper"]
        and b["gc_min"] <= gc <= b["gc_max"]
        and b["tm_min"] <= Tm <= b["tm_max"]
        and b["len_min"] <= length <= b["len_max"]
    ):
        return True
    return False


def _dedup(cands, target_tm, pos_tol=2, tm_tol=0.25):
    if len(cands) <= 1:
        return cands

    def _quality(c):
        _, _, _, seq, _, gc, Tm, s5, s3, dg_pen = c
        return _primer_score(Tm, gc, s5, s3, target_tm, seq=seq) - dg_pen

    result = []
    for rev_dir in (1, -1):
        group = sorted([c for c in cands if c[2] == rev_dir], key=lambda c: c[0])
        if not group:
            continue
        kept = [group[0]]
        for curr in group[1:]:
            last = kept[-1]
            if abs(curr[0] - last[0]) <= pos_tol and abs(curr[6] - last[6]) <= tm_tol:
                if _quality(curr) > _quality(last):
                    kept[-1] = curr
            else:
                kept.append(curr)
        result.extend(kept)
    return result


def candidate_primer(
    seq,
    params,
    out_inner,
    out_outer,
    out_loop=None,
    offset=0,
    src_name=None,
    design=None,
    sec_cache=None,
):
    deltah = params["deltah"]
    deltas = params["deltas"]
    stab = params["stab"]
    counts = [0, 0, 0]

    design = design or DEFAULT_DESIGN
    g = design["global"]
    regions = design["regions"]
    rt = design["region_thresholds"]

    MIN_LEN, MAX_LEN = g["min_len"], g["max_len"]
    TM_MIN, TM_MAX = g["tm_min"], g["tm_max"]

    inner_list = []
    outer_list = []
    loop_list = []

    if sec_cache is None:
        sec_cache = {}
    n_windows = 0
    n = len(seq)
    _t = max(
        0.0,
        min(
            1.0,
            (math.log(n) - math.log(50_000)) / (math.log(1_000_000) - math.log(50_000)),
        ),
    )
    GC_MIN = g["gc_min"] if g["gc_min"] is not None else int(30 + 20 * _t)
    GC_MAX = g["gc_max"] if g["gc_max"] is not None else int(70 - 10 * _t)

    for pos in range(n - MIN_LEN + 1):
        out_rgc_f3_p = gc_content(seq[pos : min(n, pos + 220)])
        out_rgc_f2_p = gc_content(seq[max(0, pos - 30) : min(n, pos + 190)])
        inn_rgc_f1c_p = gc_content(seq[max(0, pos - 80) : min(n, pos + 140)])
        lp_rgc_lf_p = gc_content(seq[max(0, pos - 55) : min(n, pos + 165)])

        for length in range(MIN_LEN, MAX_LEN + 1):
            if pos + length > n:
                continue

            primer = seq[pos : pos + length]
            n_windows += 1
            if n_windows % 250000 == 0:
                print(f"Processed {n_windows} windows")
            if _NON_ACGT.search(primer):
                break

            gc = gc_content(primer)
            if gc < GC_MIN or gc > GC_MAX:
                continue

            Tm = calc_tm(primer, deltah, deltas)
            if Tm < TM_MIN or Tm > TM_MAX:
                continue

            end = pos + length
            out_rgc_b3_m = gc_content(seq[max(0, end - 220) : end])
            out_rgc_b2_m = gc_content(seq[max(0, end - 190) : min(n, end + 30)])
            inn_rgc_b1c_m = gc_content(seq[max(0, end - 140) : min(n, end + 80)])
            lp_rgc_lb_m = gc_content(seq[max(0, end - 165) : min(n, end + 55)])

            fouter_ok = _region_ok(gc, Tm, length, out_rgc_f3_p, regions, "outer", rt)
            fouter_rgc = out_rgc_f3_p if fouter_ok else out_rgc_f2_p
            if not fouter_ok:
                fouter_ok = _region_ok(
                    gc, Tm, length, out_rgc_f2_p, regions, "outer", rt
                )
            f1c_ok = _region_ok(gc, Tm, length, inn_rgc_f1c_p, regions, "inner", rt)
            lf_ok = _region_ok(gc, Tm, length, lp_rgc_lf_p, regions, "loop", rt)
            bouter_ok = _region_ok(gc, Tm, length, out_rgc_b3_m, regions, "outer", rt)
            bouter_rgc = out_rgc_b3_m if bouter_ok else out_rgc_b2_m
            if not bouter_ok:
                bouter_ok = _region_ok(
                    gc, Tm, length, out_rgc_b2_m, regions, "outer", rt
                )
            b1c_ok = _region_ok(gc, Tm, length, inn_rgc_b1c_m, regions, "inner", rt)
            lb_ok = _region_ok(gc, Tm, length, lp_rgc_lb_m, regions, "loop", rt)

            if not (
                fouter_ok
                or b1c_ok
                or f1c_ok
                or bouter_ok
                or (bool(out_loop) and (lf_ok or lb_ok))
            ):
                continue

            rev = reverse_complement(primer)
            plus_ok = dimer_ok(primer)
            minus_ok = dimer_ok(rev)
            if not plus_ok and not minus_ok:
                continue

            s5p = s3p = s5m = s3m = 0
            if plus_ok:
                s5p, s3p = calc_stability(primer, stab)
            if minus_ok:
                s5m, s3m = calc_stability(rev, stab)

            fouter = plus_ok and fouter_ok and s5p >= 3 and 4 <= s3p <= 8
            bouter = (
                pos >= MIN_LEN and minus_ok and bouter_ok and s5m >= 3 and 4 <= s3m <= 8
            )
            ip = plus_ok and b1c_ok and s5p >= 4 and 3 <= s3p <= 8
            im = pos >= MIN_LEN and minus_ok and f1c_ok and s5m >= 4 and 3 <= s3m <= 8
            lp = plus_ok and lb_ok and s5p >= 3 and 4 <= s3p <= 8
            lm = pos >= MIN_LEN and minus_ok and lf_ok and s5m >= 3 and 4 <= s3m <= 8

            is_inner = ip or im
            is_outer = fouter or bouter
            is_loop = bool(out_loop) and (lp or lm)
            if not is_inner and not is_outer and not is_loop:
                continue

            if ip or fouter or lp:
                if primer not in sec_cache:
                    sec_cache[primer] = _secondary_ok(primer, Tm)
                if not sec_cache[primer][0]:
                    ip = fouter = lp = False
            if im or bouter or lm:
                if rev not in sec_cache:
                    sec_cache[rev] = _secondary_ok(rev, Tm)
                if not sec_cache[rev][0]:
                    im = bouter = lm = False

            if ip:
                inner_list.append(
                    (
                        pos,
                        length,
                        1,
                        primer,
                        inn_rgc_b1c_m,
                        gc,
                        Tm,
                        s5p,
                        s3p,
                        sec_cache[primer][1],
                    )
                )
            if im:
                inner_list.append(
                    (
                        pos,
                        length,
                        -1,
                        rev,
                        inn_rgc_f1c_p,
                        gc,
                        Tm,
                        s5m,
                        s3m,
                        sec_cache[rev][1],
                    )
                )
            if fouter:
                outer_list.append(
                    (
                        pos,
                        length,
                        1,
                        primer,
                        fouter_rgc,
                        gc,
                        Tm,
                        s5p,
                        s3p,
                        sec_cache[primer][1],
                    )
                )
            if bouter:
                outer_list.append(
                    (
                        pos,
                        length,
                        -1,
                        rev,
                        bouter_rgc,
                        gc,
                        Tm,
                        s5m,
                        s3m,
                        sec_cache[rev][1],
                    )
                )
            if out_loop:
                if lp:
                    loop_list.append(
                        (
                            pos,
                            length,
                            1,
                            primer,
                            lp_rgc_lb_m,
                            gc,
                            Tm,
                            s5p,
                            s3p,
                            sec_cache[primer][1],
                        )
                    )
                if lm:
                    loop_list.append(
                        (
                            pos,
                            length,
                            -1,
                            rev,
                            lp_rgc_lf_p,
                            gc,
                            Tm,
                            s5m,
                            s3m,
                            sec_cache[rev][1],
                        )
                    )

    if not inner_list or not outer_list:
        return counts

    _DEDUP_THRESHOLD = 6000  # if more than 6000 primers, remove close duplicates
    if len(inner_list) > _DEDUP_THRESHOLD:
        inner_list = _dedup(inner_list, target_tm=65.0)
    if len(outer_list) > _DEDUP_THRESHOLD:
        outer_list = _dedup(outer_list, target_tm=60.0)
    if loop_list and len(loop_list) > _DEDUP_THRESHOLD:
        loop_list = _dedup(loop_list, target_tm=65.0)

    inner_nbr = [0] * len(inner_list)
    if inner_list:
        outer_pos = [e[0] for e in outer_list]
        loop_pos = [e[0] for e in loop_list]
        out_lo = out_hi = lp_lo = lp_hi = 0
        for i, (pos, *_) in enumerate(inner_list):
            wlo, whi = pos - 150, pos + 150
            while out_lo < len(outer_pos) and outer_pos[out_lo] < wlo:
                out_lo += 1
            while out_hi < len(outer_pos) and outer_pos[out_hi] <= whi:
                out_hi += 1
            inner_nbr[i] += out_hi - out_lo
            while lp_lo < len(loop_pos) and loop_pos[lp_lo] < wlo:
                lp_lo += 1
            while lp_hi < len(loop_pos) and loop_pos[lp_hi] <= whi:
                lp_hi += 1
            inner_nbr[i] += lp_hi - lp_lo

    def _src_field(pos):
        return f"\tsrc:{src_name}\tlocalpos:{pos}" if src_name is not None else ""

    max_nbr = max(inner_nbr) or 1
    for i, (pos, length, rev, pseq, rgc, gc, Tm, s5, s3, dg_pen) in enumerate(
        inner_list
    ):
        base = _primer_score(Tm, gc, s5, s3, 65.0, seq=pseq) - dg_pen
        score = base + 2.0 * (inner_nbr[i] / max_nbr)
        out_inner.write(
            f"pos:{pos + offset}\tlen:{length}\trev:{rev}\tseq:{pseq}\tregGC:{rgc:.1f}\tGC:{gc:.1f}\tTm:{Tm:.2f}\tdG5:{-s5:.2f}\tdG3:{-s3:.2f}\tscore:{score:.2f}{_src_field(pos)}\n"
        )
    for pos, length, rev, pseq, rgc, gc, Tm, s5, s3, dg_pen in outer_list:
        score = _primer_score(Tm, gc, s5, s3, 60.0, seq=pseq) - dg_pen
        out_outer.write(
            f"pos:{pos + offset}\tlen:{length}\trev:{rev}\tseq:{pseq}\tregGC:{rgc:.1f}\tGC:{gc:.1f}\tTm:{Tm:.2f}\tdG5:{-s5:.2f}\tdG3:{-s3:.2f}\tscore:{score:.2f}{_src_field(pos)}\n"
        )
    if out_loop:
        for pos, length, rev, pseq, rgc, gc, Tm, s5, s3, dg_pen in loop_list:
            score = _primer_score(Tm, gc, s5, s3, 65.0, seq=pseq) - dg_pen
            out_loop.write(
                f"pos:{pos + offset}\tlen:{length}\trev:{rev}\tseq:{pseq}\tregGC:{rgc:.1f}\tGC:{gc:.1f}\tTm:{Tm:.2f}\tdG5:{-s5:.2f}\tdG3:{-s3:.2f}\tscore:{score:.2f}{_src_field(pos)}\n"
            )

    counts[0], counts[1], counts[2] = len(inner_list), len(outer_list), len(loop_list)
    return counts


_MIN_SEGMENT_LEN = 200
_WINDOW = 300
_STRIDE = 50


def parse_args():
    parser = argparse.ArgumentParser(
        prog="Single",
        usage="Single -in <ref_genome> [options]*",
        add_help=False,
    )

    parser.add_argument(
        "-in",
        dest="input",
        metavar="<ref_genome>",
        default=None,
        help="Can be MSA or single sequence fasta",
    )
    parser.add_argument(
        "-loop",
        dest="loop",
        type=lambda x: x.lower() not in ("false", "0"),
        default=True,
        metavar="<True|False>",
    )
    parser.add_argument("-par", dest="par", metavar="<par_directory>", default=None)
    parser.add_argument(
        "-design",
        dest="design",
        metavar="<design.json>",
        default=None,
        help="JSON file of primer-design criteria (region GC/Tm/length bounds",
    )
    parser.add_argument("-h", "-help", dest="help", action="store_true", default=False)

    args = parser.parse_args()

    if args.help:
        print("USAGE:")
        print("  Single -in <ref_genome> [options]*\n")
        print("ARGUMENTS:")
        print("  -in <ref_genome>")
        print("    Can be MSA or single sequence fasta")
        print("  -loop")
        print(
            "    identify candidate single primer regions for loop primers (default: on)"
        )
        print("  -par <par_directory>")
        print("    parameter file directory (default: Par/)")
        print("  -design <design.json>")
        print("    primer-design criteria file")
        sys.exit(0)

    if args.input is None:
        print("Needs sequence file")
        parser.print_usage()
        sys.exit(1)

    if args.par is not None:
        par_path = args.par
        if not par_path.endswith("/") and not par_path.endswith(os.sep):
            par_path += "/"
    else:
        par_path = os.getcwd() + "/Par/"

    return {
        "input": args.input,
        "loop": args.loop,
        "par_path": par_path,
        "design_path": args.design,
    }


def validate_alignment(sequences):
    if len(sequences) <= 1:
        return
    lengths = {len(s) for s in sequences.values()}
    if len(lengths) != 1:
        print(
            "Error: Input must be a pre-aligned MSA",
            file=sys.stderr,
        )
        sys.exit(1)


_MIN_CLUSTER_SIZE = 2
_MIN_SILHOUETTE = 0.3
_MAX_K = 4


def _silhouette(dist, labels):
    n = len(labels)
    uniq = np.unique(labels)
    if len(uniq) < 2:
        return -1.0
    s = np.zeros(n)
    for i in range(n):
        same = labels == labels[i]
        same[i] = False
        a = dist[i, same].mean() if same.any() else 0.0
        b = min(dist[i, labels == lbl].mean() for lbl in uniq if lbl != labels[i])
        s[i] = (b - a) / max(a, b) if max(a, b) > 0 else 0.0
    return float(s.mean())


def cluster_sequences(names, dist):
    n = len(names)
    if n <= 2:
        return {1: list(names)}

    Z = linkage(squareform(dist, checks=False), method="complete")
    best_k, best_score = 1, -1.0
    for k in range(2, min(_MAX_K, n - 1) + 1):
        labels = fcluster(Z, t=k, criterion="maxclust")
        sizes = np.unique(labels, return_counts=True)[1]
        if sizes.min() < _MIN_CLUSTER_SIZE:
            continue
        score = _silhouette(dist, labels)
        if score > best_score:
            best_k, best_score = k, score

    if best_score < _MIN_SILHOUETTE:
        return {1: list(names)}

    labels = fcluster(Z, t=best_k, criterion="maxclust")
    clusters = {}
    for name, label in zip(names, labels):
        clusters.setdefault(int(label), []).append(name)
    return clusters


def cluster_tightness(names, dist, members):
    """Mean pairwise similarity (1 - distance) within a cluster."""
    if len(members) < 2:
        return 1.0
    idx = [names.index(m) for m in members]
    pairs = [dist[i, j] for a, i in enumerate(idx) for j in idx[a + 1 :]]
    return 1.0 - (sum(pairs) / len(pairs))


def _window_distance(seq_arrs, lo, hi):
    """Pairwise Hamming distance"""
    n = len(seq_arrs)
    width = hi - lo
    dist = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            d = np.count_nonzero(seq_arrs[i][lo:hi] != seq_arrs[j][lo:hi]) / width
            dist[i, j] = dist[j, i] = d
    return dist


def _qualifying_clusters(names, dist):
    clusters = cluster_sequences(names, dist)
    scored = [
        (frozenset(members), cluster_tightness(names, dist, members))
        for members in clusters.values()
    ]
    scored.sort(key=lambda t: (len(t[0]), t[1]), reverse=True)
    return scored


def _col_to_real(aligned_seq):
    mapping = [None] * len(aligned_seq)
    real_i = 0
    for col, ch in enumerate(aligned_seq):
        if ch != "-":
            mapping[col] = real_i
            real_i += 1
    return mapping


def _real_subseq(raw_seq, col_map, start, end):
    real_positions = [col_map[c] for c in range(start, end) if col_map[c] is not None]
    if not real_positions:
        return ""
    return raw_seq[real_positions[0] : real_positions[-1] + 1]


def _real_to_col(col_map):
    out = {}
    for col, r in enumerate(col_map):
        if r is not None:
            out[r] = col
    return out


def _apply_compatibility(
    text, rep_name, real_to_col, names, col_maps, ungapped, win_start
):
    """identify all genomes that contain primer's exact sequence"""
    if not text:
        return text
    out_lines = []
    for line in text.splitlines():
        if not line:
            continue
        parts = line.split("\t")
        fields = {}
        for i, part in enumerate(parts):
            k, v = part.split(":", 1)
            fields[k] = (i, v)

        real_pos = win_start + int(fields["localpos"][1])
        length = int(fields["len"][1])
        primer_seq = fields["seq"][1]

        start_col = real_to_col[real_pos]
        end_col = real_to_col[real_pos + length - 1] + 1
        compatible = [rep_name]
        for n in names:
            if n == rep_name:
                continue
            if _real_subseq(ungapped[n], col_maps[n], start_col, end_col) == primer_seq:
                compatible.append(n)

        src_idx = fields["src"][0]
        parts[src_idx] = "src:" + ",".join(sorted(compatible))
        out_lines.append("\t".join(parts))
    return "\n".join(out_lines) + "\n" if out_lines else ""


def _primer_cov(fields, real_to_col, names, col_maps, ungapped, win_start):
    real_pos = win_start + int(fields["localpos"][1])
    length = int(fields["len"][1])
    primer_seq = fields["seq"][1]
    start_col = real_to_col[real_pos]
    end_col = real_to_col[real_pos + length - 1] + 1

    covered = sum(
        _real_subseq(ungapped[n], col_maps[n], start_col, end_col) == primer_seq
        for n in names
    )
    return covered / len(names)


def _apply_cov_tag(text, real_to_col, names, col_maps, ungapped, win_start):
    if not text:
        return text
    out_lines = []
    for line in text.splitlines():
        if not line:
            continue
        parts = line.split("\t")
        fields = {}
        for i, part in enumerate(parts):
            k, v = part.split(":", 1)
            fields[k] = (i, v)
        cov = _primer_cov(fields, real_to_col, names, col_maps, ungapped, win_start)
        parts.append(f"cov:{cov:.4f}")
        out_lines.append("\t".join(parts))
    return "\n".join(out_lines) + "\n" if out_lines else ""


def _apply_flat_cov_tag(text, cov):
    if not text:
        return text
    out_lines = [f"{line}\tcov:{cov:.4f}" for line in text.splitlines() if line]
    return "\n".join(out_lines) + "\n" if out_lines else ""


def _apply_cluster_info(text, cluster_idx, n_clusters):
    if not text:
        return text
    tag = f"\tcluster:{cluster_idx}\tnclusters:{n_clusters}"
    return "".join(line + tag + "\n" for line in text.splitlines() if line)


def _dedup_lines(text, seen):
    if not text:
        return text
    out_lines = []
    for line in text.splitlines():
        fields = dict(part.split(":", 1) for part in line.split("\t"))
        key = (fields["rev"], fields["seq"], fields["pos"])
        if key in seen:
            continue
        seen.add(key)
        out_lines.append(line)
    return "\n".join(out_lines) + "\n" if out_lines else ""


_SKETCH_SIZE = 2000
_SKETCH_THRESHOLD = 100_000


def _choose_k(length, margin=2, k_min=8, k_max=24):
    """Don't want too large of a k for a small sequence."""
    k = round(math.log2(max(length, 2))) + margin
    return max(k_min, min(k_max, k))


def _hash_kmer(kmer):
    return int.from_bytes(hashlib.blake2b(kmer.encode(), digest_size=8).digest(), "big")


def _kmer_profile(seq, k, use_sketch):
    hashes = {_hash_kmer(seq[i : i + k]) for i in range(len(seq) - k + 1)}
    if use_sketch and len(hashes) > _SKETCH_SIZE:
        return set(sorted(hashes)[:_SKETCH_SIZE])
    return hashes


def _jaccard(a, b, use_sketch):
    if not a or not b:
        return 0.0
    if not use_sketch:
        return len(a & b) / len(a | b)
    union_bottom = sorted(a | b)[:_SKETCH_SIZE]
    if not union_bottom:
        return 0.0
    in_both = sum(1 for x in union_bottom if x in a and x in b)
    return in_both / len(union_bottom)


def pick_most_similar(sequences):
    """Ranks based on kmer similarity"""
    names = list(sequences.keys())
    if len(names) == 1:
        name = names[0]
        return name, sequences[name], None

    max_len = max(len(s) for s in sequences.values())
    k = _choose_k(max_len)
    use_sketch = max_len > _SKETCH_THRESHOLD
    profiles = {n: _kmer_profile(s, k, use_sketch) for n, s in sequences.items()}

    avg_sim = {}
    for ni in names:
        total = sum(
            _jaccard(profiles[ni], profiles[nj], use_sketch) for nj in names if nj != ni
        )
        avg_sim[ni] = total / (len(names) - 1)

    best_name = max(avg_sim, key=avg_sim.get)
    return best_name, sequences[best_name], avg_sim[best_name]


if __name__ == "__main__":
    import time

    t0 = time.time()

    args = parse_args()
    sequences = read_fasta(args["input"], strip_gaps=False)
    params = load_params(args["par_path"])
    design = load_design(args["design_path"])
    validate_alignment(sequences)

    os.makedirs("Results", exist_ok=True)

    f_inner = open(os.path.join("Results", "InnerCandidates"), "w")
    f_outer = open(os.path.join("Results", "OuterCandidates"), "w")
    f_loop = (
        open(os.path.join("Results", "LoopCandidates"), "w") if args["loop"] else None
    )

    total = [0, 0, 0]

    if len(sequences) == 1:
        for _, seq in sequences.items():
            buf_loop = io.StringIO() if f_loop else None
            counts = candidate_primer(
                seq, params, f_inner, f_outer, buf_loop, design=design
            )
            if f_loop:
                f_loop.write(_apply_flat_cov_tag(buf_loop.getvalue(), 1.0))
            total[0] += counts[0]
            total[1] += counts[1]
            total[2] += counts[2]
        f_inner.close()
        f_outer.close()
        if f_loop:
            f_loop.close()
    else:
        names = list(sequences.keys())
        aln_len = len(next(iter(sequences.values())))
        seq_arrs = [np.frombuffer(sequences[n].encode(), dtype=np.uint8) for n in names]
        col_maps = {n: _col_to_real(sequences[n]) for n in names}
        ungapped = {n: sequences[n].replace("-", "") for n in names}

        window, stride = _WINDOW, _STRIDE

        seq_offset = {}
        cum = 0
        for n in names:
            seq_offset[n] = cum
            cum += len(ungapped[n])

        seen_inner, seen_outer, seen_loop = set(), set(), set()
        sec_cache = {}
        w = 0
        while w < aln_len:
            hi = min(w + window, aln_len)
            if hi - w < max(window // 2, 1):
                break

            dist = _window_distance(seq_arrs, w, hi)
            clusters = _qualifying_clusters(names, dist)
            n_clusters = len(clusters)

            for cluster_idx, (members, _) in enumerate(clusters):
                subset = {
                    n: _real_subseq(ungapped[n], col_maps[n], w, hi)
                    for n in sorted(members)
                }
                rep_name, rep_seq, _sim = pick_most_similar(subset)
                if len(rep_seq) < _MIN_SEGMENT_LEN:
                    continue
                src_name = f"{rep_name}_w{w}-{hi}_c{cluster_idx}"
                real_start = next(
                    col_maps[rep_name][c]
                    for c in range(w, hi)
                    if col_maps[rep_name][c] is not None
                )
                offset = seq_offset[rep_name] + real_start

                buf_inner = io.StringIO()
                buf_outer = io.StringIO()
                buf_loop = io.StringIO() if f_loop else None
                candidate_primer(
                    rep_seq,
                    params,
                    buf_inner,
                    buf_outer,
                    buf_loop,
                    offset=offset,
                    src_name=src_name,
                    sec_cache=sec_cache,
                    design=design,
                )
                real_to_col = _real_to_col(col_maps[rep_name])

                def _finalize(buf, seen, is_loop=False):
                    text = _apply_compatibility(
                        buf.getvalue(),
                        rep_name,
                        real_to_col,
                        names,
                        col_maps,
                        ungapped,
                        real_start,
                    )
                    if is_loop:
                        text = _apply_cov_tag(
                            text, real_to_col, names, col_maps, ungapped, real_start
                        )
                    text = _apply_cluster_info(text, cluster_idx, n_clusters)
                    return _dedup_lines(text, seen)

                inner_text = _finalize(buf_inner, seen_inner)
                outer_text = _finalize(buf_outer, seen_outer)
                loop_text = (
                    _finalize(buf_loop, seen_loop, is_loop=True) if f_loop else ""
                )
                f_inner.write(inner_text)
                f_outer.write(outer_text)
                if f_loop:
                    f_loop.write(loop_text)
                total[0] += inner_text.count("\n")
                total[1] += outer_text.count("\n")
                total[2] += loop_text.count("\n")

            w += stride

        f_inner.close()
        f_outer.close()
        if f_loop:
            f_loop.close()

    print(f"Total candidates: {total[0]} inner, {total[1]} outer, {total[2]} loop")
    print(f"Done in {time.time() - t0:.1f}s")
