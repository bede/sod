# Simulator of Diversity

 Simply generates subsequences of a specified length from a supplied reference sequence, applying base insertions, deletions and substitutions at specified rates

Simulates:
- Sequencing of a randomly mutated reference genome
- Single reads of a defined length
- Substitutions, insertions and deletions at specified rates
    - Insertions and deletions of geometrically distributed length (default p=0.5)
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
```
$ ./sod.py -h
usage: sod.py [-h] [-n N_READS] [-l LEN_READS] [-s SUB_RATE]
              [--ins-rate INS_RATE] [-d DEL_RATE]
              [--indel-ext-rate INDEL_EXT_RATE] [-f] [-u]
              path-to-ref

positional arguments:
  path-to-ref           -

optional arguments:
  -h, --help            show this help message and exit
  -n N_READS, --n-reads N_READS
                        1
  -l LEN_READS, --len-reads LEN_READS
                        250
  -s SUB_RATE, --sub-rate SUB_RATE
                        0.0
  --ins-rate INS_RATE   0.0
  -d DEL_RATE, --del-rate DEL_RATE
                        0.0
  --indel-ext-rate INDEL_EXT_RATE
                        0.5
  -f, --fastq           False
  -u, --uppercase       False
$ ./sod.py -n 1 -l 250 -s 0.1 -i 0.05 -d 0.05 tests/hxb2.fasta
>read1_r7937_s0_i27_d0
CCAGACTGTGAGcTTGCaAACgAGATGCTGTTGCGCCTCAAtTAaGCCCTCAGCAAATcT
GtTTCTGCTGCTGCACTAaTACCAGaACAAccTAATTGaTCTGGCCTGTAgCCGTCAGCa
GTCATTGAGGCgTGCGCCCATAGTGCTTCCgTGcCTGaCTCCCAAtGAACCCAAGGAACA
AAGCTCCaTATTCCagCAaCTGCtTCTTTTTTCTCTCTGCACCaACTCTTCcTCTTTGCC
TTtGgGaTGG
```