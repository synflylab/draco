# DRACO
### Direct Repeat Aware Codon Optimizer

DRACO aims to ease the design and synthesis of gene fragments encoding for
highly repetitive protein sequence, such as TALE DNA binding domains by
utilizing repeat-aware codon optimization.

```
usage: draco [-h] [--upstream UPSTREAM] [--downstream DOWNSTREAM]
             [--check-repeats CHECK_REPEATS] [--repeat-len REPEAT_LEN]
             [--check-inv-repeats CHECK_INV_REPEATS]
             [--inv-repeat-len INV_REPEAT_LEN] [--check-stretch CHECK_STRETCH]
             [--stretch-len STRETCH_LEN] [--check-gc CHECK_GC]
             [--max-gc MAX_GC] [--gc-window GC_WINDOW]
             [--repeat-attempts REPEAT_ATTEMPTS] [--seq-attempts SEQ_ATTEMPTS]
             [--handicap HANDICAP] [--progress] [--log LOG] [--debug DEBUG]
             sequence

Direct Repeat Aware Codon Optimizer

positional arguments:
  sequence              TALE binding sequence (DNA)

optional arguments:
  -h, --help            show this help message and exit
  --upstream UPSTREAM   Protein sequence to include upstream of the repeats
  --downstream DOWNSTREAM
                        Protein sequence to include downstream of the repeats
  --check-repeats CHECK_REPEATS
                        Check the optimized sequence for repeats
  --repeat-len REPEAT_LEN
                        Max allowed length of a repeat
  --check-inv-repeats CHECK_INV_REPEATS
                        Check the optimized sequence for inverse repeats
  --inv-repeat-len INV_REPEAT_LEN
                        Max allowed length of an inverted repeat
  --check-stretch CHECK_STRETCH
                        Check for base stretches
  --stretch-len STRETCH_LEN
                        Max allowed length of a base stretch
  --check-gc CHECK_GC   Check GC content
  --max-gc MAX_GC       Max allowed GC content
  --gc-window GC_WINDOW
                        GC calculation window
  --repeat-attempts REPEAT_ATTEMPTS
                        Max attempts to optimize a repeat
  --seq-attempts SEQ_ATTEMPTS
                        Max attempts to optimize the sequence
  --handicap HANDICAP   Max codon handicap for the optimized sequence
  --progress            Show progress indicator
  --log LOG             Logging level
  --debug DEBUG         Enable debug logging
```