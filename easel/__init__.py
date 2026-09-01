"""easel - selection screens on thousands of coding sequences.

Layout
------
  cli.py        argument parsing, DEF.yaml construction, snakemake invocation
  validate.py   flag-combination and path checks, run before anything happens
  formats.py    validation of the biological input files themselves
  preprocess.py per-transcript prefiltering before the workflow starts
  runner.py     lives in cli.py for now; split out when it grows

The workflow files (Snakefile_*, rules/, envs/, scripts/, prof/) ship as
package data. scripts/ is deliberately NOT importable from here: those run
inside per-rule conda environments that do not have easel installed, so they
share code through scripts/_seqio.py instead.
"""

__version__ = "0.1.0"
