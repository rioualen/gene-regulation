"""Snakefile to test basic functions of snakemake.
"""

rule hello:
   """Write HELLO in a text file named hello.txt.
"""
    output: "hello.txt"
    benchmark: "hello_benchmark.json"
    shell: "echo HELLO > {output}"

rule bye:
   """Write BYE in a text file named bye.txt.
"""
    input: "hello.txt"
    output: "bye.txt"
    benchmark: "bye_benchmark.json"
    shell: "echo BYE > {output}"

TO_CLEAN="hello.txt bye.txt".split()
rule clean:
   """Delete the files hello.txt and bye.txt
"""
    shell: "rm -f {TO_CLEAN}"
