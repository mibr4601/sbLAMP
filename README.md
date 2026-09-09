# sbLAMP

sbLAMP is a score-based LAMP primer set designer. Given a single target sequence
or an MSA, it produces sets of 4–6 primers. For an MSA it can also allow
degenerate primers and multiplex-compatible sets.

`interface.py` is a simple way to run `single.py` and `lamp.py`. `single.py`
runs first and creates candidate regions; `lamp.py` then builds sets from those
candidate regions.

## Files

- `interface.py` — GUI wrapper for `single.py` and `lamp.py`


- `single.py` — creates candidate regions from the target sequence(s)


- `lamp.py` — combines candidate regions into primer sets


- `repairSet.py` — degenerate primer repairs


- `Par/santalucia_tm_nn_parameters.txt` — nearest-neighbour parameters for Tm calculation


- `Results/`
  - `InnerCandidates`, `OuterCandidates`, `LoopCandidates` — candidate regions written by `single.py`
  - `results.txt` — default location for the final sets


- `example/` — example sequences from NCBI


- `single_no_cluster.py`, `lamp_no_cluster.py` — generate sets from a single
  reference sequence instead of clustering (used to test `single.py`)

## How to run

Install dependencies:

```
pip install -r requirements.txt
```

Run all commands from the repo root.

### interface.py (GUI)

```
python interface.py
```

### single.py (stage 1)

Identifies candidate primer regions.

```
python single.py -in <input.fasta> [-loop True|False] [-par Par/]
```

Arguments:

- `-in` — a single sequence or MSA FASTA file
- `-loop` — default `True`; when true, loop primers are designed with the core set
- `-par` — default `Par/`; the directory must contain `santalucia_tm_nn_parameters.txt`

### lamp.py (stage 2)

Combines candidate regions into full sets. Reads candidates from
`Results/InnerCandidates`, `Results/OuterCandidates`, `Results/LoopCandidates`.

```
python lamp.py -in <input.fasta> [-num 10] [-loop True|False] [-gblock True|False]
               [-multiplex 0] [-max_repairs 0] [-min_coverage 0] [-check 1]
               [-out Results/results.txt] [-par Par/]
```

Arguments:

- `-num` — number of sets to return
- `-loop` — include LF/LB loop primers
- `-gblock` — require a gBlock template per set (default `True`)
- `-multiplex N` — return bundles of up to N mutually compatible sets chosen for combined coverage
- `-max_repairs N` — allow up to N degenerate bases per set to recover uncovered targets
- `-min_coverage P` — drop sets below P% target coverage
- `-check 0` — skip secondary-structure checks (on by default)

## Comparison: single reference

`single_no_cluster.py` and `lamp_no_cluster.py` skip clustering and generate
candidates from a single reference sequence. Everything else is the same.

- `-ref 0` — reference to design from; `0` picks the most similar reference,
  any other number picks the Nth reference in the FASTA file

```
python single_no_cluster.py -in <MSA.fasta> -ref 0
python lamp_no_cluster.py   -in <MSA.fasta> -ref 0
```

## Output format (`Results/results.txt`)

```
10 set(s)  avg coverage=86.7%  best coverage=100.0%

LAMP set 1  score=17.75  amplicon 274-508 (234 bp)  regional GC=53.4%  coverage=100.0%  source:AM922330.1
  F3 : pos:274  len:20 bp  Tm:59.3C  GC:55.0%  seq(5'-3'):CGCTTAACTGGTCTGAGAGG
  F2 : pos:305  len:19 bp  Tm:59.3C  GC:57.9%  seq(5'-3'):CACTGGAACTGAGACACGG
  F1c: pos:361  len:20 bp  Tm:64.4C  GC:60.0%  seq(5'-3'):TCAGGGTTTCCCCCATTGCG
  B1c: pos:391  len:20 bp  Tm:64.1C  GC:60.0%  seq(5'-3'):CCGCGTGGAGGATGACACTT
  B2 : pos:450  len:21 bp  Tm:57.5C  GC:47.6%  seq(5'-3'):GTGCTTATTCCTTAGGTACCG
  B3 : pos:490  len:18 bp  Tm:60.8C  GC:61.1%  seq(5'-3'):CTCCGTATTACCGCGGCT
  LF : pos:334  len:20 bp  Tm:64.8C  GC:65.0%  seq(5'-3'):CCCTACTGCTGCCTCCCGTA
  LB : none found
  gblock_1 region: 224-523(299)  ambiguous bases:2  ambiguous base positions:[202, 203]
  gblock_1: ACTATATAGTATCAGCTAGTTGGTGAGGTAATGGCTCACCAAGGCTATGACGC...
```
