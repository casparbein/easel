# Changelog

All notable changes to this project are documented here.
Format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [0.1.0] — 2026-09-09

First tagged pre-release.

### Added
- Two-phase run: alignment/validation, followed by trees and selection screens
- Free mode (`-free`) and TOGA2 mode (`-toga2`)
- Aligners: PRANK (codon/nucleotide), MACSE2, MUSCLE5
- Selection screens: aBSREL, BUSTED (+error-sink), MEME, RELAX, BayesCode
- Per-rule conda environments and container-based HmmCleaner and TOGA2 instances

### Known limitations
- `-toga` (TOGA v1) is not implemented; easel exits if given
- `-st` (ASTRAL species tree) is not implemented; easel exits if given
- `pip install .` (non-editable) does not yet locate the workflow files —
  use `pip install -e .` as documented
