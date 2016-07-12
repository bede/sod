#!/usr/bin/env python3

from sod import sod

output = sod.simulate_reads('tests/hxb2.fasta', n_reads=1, len_reads=250,
         sub_rate=0.0, ins_rate=0.0, del_rate=0.0,
         fastq=False, uppercase=False)

def test_simulate():
    assert output