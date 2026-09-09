# sbLAMP

sbLAMP is a score-based LAMP primer set designer. Given a single target sequence or MSA, it creates between sets of between 4-6 primers. For MSA, there is the option of allowing degenerate primers and multiplex compatible sets.

interface.py is provided as a simple way to run single.py and lamp.py. Single.py is run first, creating candidate regions. Lamp.py follows by creating sets from the candidate regions.

## Files

interface.py
    Wrapper for single.py, lamp.py
single.py
    Creates candiddate regions from target sequence(s)
lamp.py
    Combines candidate regions into primer sets
repairSet.py
    Degenerate primer repairs
Par/santalucia_tm_nn_parameters.txt
    Numbers for Tm calculation
Results/
    InnerCandidates, OuterCandidates, LoopCandidates, 
        Inner candidates generated from single.py
    results.txt
        Default location to save results
example/
    Example sequences from NCBI
single_no_cluster.py, lamp_no_cluster.py
    Generate sets using single reference sequence instead of clustering. (Used to test single.py)
    


## How to run

Install dependencies with 
    pip install -r requirements.txt
and run all commands from the repo root 

### interface.py (GUI)

python interface.py

### single.py (stage 1)
Identifies candidate primer regions

python single.py -in <input.fasta> [-loop True|False] [-par Par/]

Arguments:
-in
    Input is a single sequence or MSA fasta file
-loop
    Defaults to true. If true, loop primers are designed with core set
-par
    Defaults to Par/. Location must contain santalucia_tm_nn_parameters.txt


### lamp.py (stage 2)
Combines candidate regions into full sets. Reads in candidates from Results/InnerCandidates, OuterCandidates, LoopCandidates

python lamp.py -in <input.fasta> [-num 10] [-loop True|False] [-gblock True|False]
               [-multiplex 0] [-max_repairs 0] [-min_coverage 0] [-check 1]
               [-out Results/results.txt] [-par Par/]

Arguments

-num
    number of sets to return
-loop
    include LF/LB loop primers
-gblock
    require a gBlock template per set. Defaults true
-multiplex N
    return bundles of up to N mutually compatible sets chosen for combined coverage
-max_repairs N
    allow up to N degenerate bases per set to recover uncovered targets
-min_coverage P
    drop sets below P% target coverage
-check 0
    skip secondary-structure checks. Defaults false

## Comparison: Single Reference

single_no_cluster.py and lamp_no_cluster.py 

No clustering, generates candidate based on single reference sequence. Rest stays the same

-ref 0
    Single reference to design from. 0 picks most similar reference, any other number picks the Nth reference in the fasta file.


python single_no_cluster.py -in <MSA.fasta> -ref 0
python lamp_no_cluster.py   -in <MSA.fasta> -ref 0

## Output format (`Results/results.txt`)

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
