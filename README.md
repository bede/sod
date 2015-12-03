# Simulator of Diversity.

 Generates subsequences of a specified length from a supplied reference sequence, applying base insertions, deletions and substitutions at specified rates

Simulates:
- Sequencing of a randomly mutated reference genome
- Single end reads of a defined length
- Substitutions, insertions and deletions at specified rates
    - Insertions and deletions of geometrically distributed length (p=0.5)
        - Insertions of random composition 
- Forward and reverse complement reads

Does not simulate:
- Read quality
- Realistic sequencing error profiles
- CIGAR values

## Requirements
- Python3
- Argh (`pip install argh`)
- NumPy (`pip install numpy`)
- Biopython (`pip install biopython`)

## Example
$ ./sod.py 
usage: sod.py [-h] [-n N_READS] [-l LEN_READS] [-s SUB_RATE] [-i INS_RATE]
              [-d DEL_RATE]
              path-to-ref
$ ./sod.py -n 1 -l 250 -s 0.1 -i 0.05 -d 0.05 tests/hxb2.fasta
>fwd_s2328_e2568
CAGGAGCATGTACAGTTAGAGAAgATGAGTTTGCtgacaAGGAAGATGGAACCAAAAATG
ATAGGGGgcGAATTGAGGTTTATCAAAAAGAGTTCAGATACTCAAGAAACGTGGACATAA
gaAATAGGTACAGATTAcTAGGACCAACCTGTAAATAATTGGAAAAACTGTTGACTCAGA
TTGTTGCACTgcTTAATTTCCATTAggcattGCCTAtgctTTGAGA