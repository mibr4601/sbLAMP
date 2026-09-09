from dataclasses import dataclass

from single import calc_tm, _secondary_ok, _pair_ok

_MAX_ROUNDS = 40  # Safety cap
_MAX_TARGET_MM = 5
_FOLD_BUDGET = 8
_TM_TOL = 2.0
_TM_CEIL = 5.0

IUPAC_BY_BASES = {
    frozenset("A"): "A",
    frozenset("C"): "C",
    frozenset("G"): "G",
    frozenset("T"): "T",
    frozenset("AG"): "R",
    frozenset("CT"): "Y",
    frozenset("GC"): "S",
    frozenset("AT"): "W",
    frozenset("GT"): "K",
    frozenset("AC"): "M",
    frozenset("CGT"): "B",
    frozenset("AGT"): "D",
    frozenset("ACT"): "H",
    frozenset("ACG"): "V",
    frozenset("ACGT"): "N",
}
IUPAC_TO_BASES = {v: set(k) for k, v in IUPAC_BY_BASES.items()}


@dataclass
class PrimerInfo:
    role: str
    pos: int
    length: int
    rev: int
    seq: str
    Tm: float
    GC: float


@dataclass
class MismatchInfo:
    pos: list
    target_bases: list


@dataclass
class RepairConfig:
    max_mismatches: int = 2
    coverage_goal: float = 90.0
    tm_penalty_scale: float = 1.5
    repair_budget: int = 0
    compat_threshold: float = 45.0


@dataclass
class RepairResult:
    primers: dict
    modified_roles: list
    unresolved_targets: set
    hits: dict
    thermo_penalty: float = 0.0


def _set_coverage_from_hits(hits, target_names):
    all_ok = {name: True for name in target_names}
    for per_target in hits.values():
        for name in target_names:
            all_ok[name] = all_ok[name] and per_target.get(name, False)
    return all_ok


def _eligible_and_excluded(mismatches, max_target_mm):
    eligible = {}
    excluded = set()
    for name, m in mismatches.items():
        if m.pos is None or len(m.pos) > max_target_mm:
            excluded.add(name)
        else:
            eligible[name] = m
    return eligible, excluded


def _fold_for_positions(seq, pos, mismatches):
    fold = 1
    for p in pos:
        bases = set(IUPAC_TO_BASES.get(seq[p], {seq[p]}))
        for m in mismatches.values():
            if p in m.pos:
                c = m.target_bases[m.pos.index(p)]
                bases.update(IUPAC_TO_BASES.get(c, {c}))
        fold *= len(bases)
    return fold


def _best_degen_positions(
    seq, mismatches, max_mismatches, fold_budget, protect_3prime
):
    primer_len = len(seq)
    candidate_pos = sorted(
        {
            p
            for m in mismatches.values()
            for p in m.pos
            if p < primer_len - protect_3prime
        }
    )

    def rescued_set(chosen):
        return {
            name
            for name, m in mismatches.items()
            if sum(1 for p in m.pos if p not in chosen) <= max_mismatches
        }

    best_pos = ()
    best_rescued = rescued_set(set())
    best_key = (len(best_rescued), 0)

    max_size = min(len(candidate_pos), max(0, fold_budget.bit_length() - 1))

    chosen = set()
    remaining = list(candidate_pos)
    while len(chosen) < max_size and remaining:
        step_p = None
        step_key = None
        step_rescued = None
        for p in remaining:
            trial = chosen | {p}
            fold = _fold_for_positions(seq, trial, mismatches)
            if fold > fold_budget:
                continue
            rescued = rescued_set(trial)
            key = (len(rescued), -fold)
            if step_key is None or key > step_key:
                step_key = key
                step_p = p
                step_rescued = rescued
        if step_p is None or step_key <= best_key:
            break
        chosen.add(step_p)
        remaining.remove(step_p)
        best_key = step_key
        best_pos = tuple(sorted(chosen))
        best_rescued = step_rescued

    return best_pos, best_rescued


def _apply_degen(seq, pos, mismatches):
    chars = list(seq)
    for p in pos:
        bases = set(IUPAC_TO_BASES.get(seq[p], {seq[p]}))
        for m in mismatches.values():
            if p in m.pos:
                c = m.target_bases[m.pos.index(p)]
                bases.update(IUPAC_TO_BASES.get(c, {c}))
        chars[p] = IUPAC_BY_BASES[frozenset(bases)]
    return "".join(chars)


_MAX_EXP = 24


def _expand_iupac(seq):
    expansions = [""]
    for c in seq:
        bases = IUPAC_TO_BASES.get(c, {c})
        expansions = [e + b for e in expansions for b in bases]
        if len(expansions) > _MAX_EXP:
            return []
    return expansions


def _verify_thermo(seq, params, ref_tm, tm_tolerance, tm_ceiling, tm_penalty_scale):
    degen_pos = [i for i, c in enumerate(seq) if c not in "ACGT"]
    if not degen_pos:
        return True, 0.0

    expansions = _expand_iupac(seq)
    if not expansions:
        return False, 0.0

    max_shift = 0.0
    for exp in expansions:
        tm = calc_tm(exp, params["deltah"], params["deltas"])
        shift = abs(tm - ref_tm)
        if shift > tm_ceiling:
            return False, 0.0
        max_shift = max(max_shift, shift)
        ok, _ = _secondary_ok(exp, tm)
        if not ok:
            return False, 0.0
    penalty = tm_penalty_scale * max(0.0, max_shift - tm_tolerance)
    return True, penalty


def _compatible_with_set(primers, role, seq, threshold):
    for other_role, p in primers.items():
        if other_role == role:
            continue
        if not _pair_ok(seq, p.seq, threshold):
            return False
    return True


def _round_candidates(
    mismatch_data,
    primers,
    hits,
    role_excluded,
    failing,
    config,
    remaining_budget,
):
    candidates = []
    free_progress = False

    for role, role_mismatches in mismatch_data.items():
        role_failing = {
            n: role_mismatches[n]
            for n in failing
            if n in role_mismatches
            and not hits[role].get(n, False)
            and n not in role_excluded[role]
        }
        if not role_failing:
            continue

        eligible, excluded = _eligible_and_excluded(role_failing, _MAX_TARGET_MM)
        role_excluded[role].update(excluded)
        if not eligible:
            continue

        pos, rescued = _best_degen_positions(
            primers[role].seq,
            eligible,
            config.max_mismatches,
            _FOLD_BUDGET,
            3,
        )
        if not rescued:
            continue

        if not pos:
            for name in rescued:
                hits[role][name] = True
            free_progress = True
            continue

        cost = len(pos)
        if cost > remaining_budget:
            continue
        gain = len(rescued)
        if gain <= 0:
            continue
        candidates.append((gain / cost, gain, role, pos, rescued, cost, eligible))

    return candidates, free_progress


def repair_set(
    primers,
    hits,
    mismatch_data,
    params,
    config,
):
    target_names = sorted({n for per_target in hits.values() for n in per_target})

    primers = dict(primers)
    hits = {role: dict(per_target) for role, per_target in hits.items()}
    modified = []
    unresolved = set()
    role_excluded = {role: set() for role in mismatch_data}
    remaining_budget = config.repair_budget
    total_thermo_penalty = 0.0

    for _ in range(_MAX_ROUNDS):
        all_ok = _set_coverage_from_hits(hits, target_names)
        ok_set = {n for n, ok in all_ok.items() if ok}
        coverage = 100.0 * len(ok_set) / len(target_names) if target_names else 100.0
        if coverage >= config.coverage_goal:
            break

        failing = [n for n in target_names if n not in ok_set]
        if not failing:
            break

        candidates, free_progress = _round_candidates(
            mismatch_data,
            primers,
            hits,
            role_excluded,
            failing,
            config,
            remaining_budget,
        )
        if free_progress:
            continue
        if not candidates:
            break

        candidates.sort(key=lambda c: (c[0], c[1]), reverse=True)

        applied = False
        for _, _, role, pos, rescued, cost, eligible in candidates:
            new_seq = _apply_degen(primers[role].seq, pos, eligible)
            thermo_ok, thermo_penalty = _verify_thermo(
                new_seq,
                params,
                primers[role].Tm,
                _TM_TOL,
                _TM_CEIL,
                config.tm_penalty_scale,
            )
            if not thermo_ok:
                continue

            variants = _expand_iupac(new_seq)
            if not variants or not all(
                _compatible_with_set(primers, role, v, config.compat_threshold)
                for v in variants
            ):
                continue

            primers[role] = PrimerInfo(
                role=role,
                pos=primers[role].pos,
                length=primers[role].length,
                rev=primers[role].rev,
                seq=new_seq,
                Tm=primers[role].Tm,
                GC=primers[role].GC,
            )
            total_thermo_penalty += thermo_penalty
            for name in rescued:
                hits[role][name] = True
            modified.append(role)
            remaining_budget -= cost
            applied = True
            break

        if not applied:
            break

    for name, ok in _set_coverage_from_hits(hits, target_names).items():
        if not ok:
            unresolved.add(name)

    return RepairResult(
        primers,
        modified,
        unresolved,
        hits,
        total_thermo_penalty,
    )
