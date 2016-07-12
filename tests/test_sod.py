#!/usr/bin/env python3

from sod import sod

len_reads = 250

records = sod.simulate_reads('tests/hxb2.fasta', n_reads=1000, len_reads=len_reads,
         				    sub_rate=0.1, ins_rate=0.1, del_rate=0.1,
                            fastq=False, uppercase=False)

def test_simulate():
	assert records

def test_read_lens():
	for record in records:
		assert len(record.seq) == len_reads


