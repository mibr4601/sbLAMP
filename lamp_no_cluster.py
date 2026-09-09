import argparse
import os
import sys
import time

from single_no_cluster import (
    read_fasta,
    validate_alignment,
    load_params,
    _col_to_real,
    select_reference,
)
from lamp import load_candidates, assemble, write_results


def parse_args():
    p = argparse.ArgumentParser(
        prog="LAMPNoCluster",
        usage="LAMPNoCluster -in <ref_genome> [options]*",
        add_help=False,
    )
    p.add_argument("-in", dest="input", metavar="<ref_genome>", default=None)
    p.add_argument(
        "-ref",
        dest="ref",
        type=int,
        default=0,
        metavar="<int>",
        help="Reference sequence: 0 picks the most central one, N picks Nth. "
        " Nth Must match single_no_cluster.",
    )
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
        print("  LAMPNoCluster -in <ref_genome> [options]*\n")
        print("ARGUMENTS:")
        print("  -in <ref_genome>   Reference genome (FASTA format)")
        print("  -ref <int>         0 picks the sequence closest to all the others;")
        print("                     N picks the Nth. Must match single_no_cluster")
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
        "ref": args.ref,
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

    ungapped = {n: s.replace("-", "") for n, s in gapped_sequences.items()}

    if len(gapped_sequences) > 1:
        validate_alignment(gapped_sequences)
        ref_name, ref_index = select_reference(gapped_sequences, args["ref"])
        seq = ungapped[ref_name]
        segments = [(0, len(seq))]
        targets = ungapped
        col_maps = {name: _col_to_real(s) for name, s in gapped_sequences.items()}
        seq_names = [ref_name]
        print(f"Reference: {ref_name}  (#{ref_index + 1} of {len(gapped_sequences)})")
    else:
        seq = next(iter(ungapped.values()))
        segments = [(0, len(seq))]
        targets = None
        col_maps = None
        seq_names = None

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
        print(
            "No candidates found in Results/. Run single_no_cluster.py first.",
            file=sys.stderr,
        )
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
    main()
