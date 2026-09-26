## Submission

First CRAN submission of sbwadjust (0.2.0).

## Test environments

* local macOS (R release): R CMD check --as-cran
* win-builder (devel): R Under development (unstable) (2026-09-21 r90579 ucrt)
* R-hub (R-devel): linux, windows, macos, macos-arm64 -- all OK
* GitHub Actions: macOS, Windows (R release); Ubuntu (R devel, release,
  oldrel-1) -- all OK

## R CMD check results

0 errors | 0 warnings | 1 note

The note is the CRAN incoming feasibility check, covering two items:

* "New submission" -- this is the package's first release.

* "Possibly misspelled words in DESCRIPTION: Luedtke, SBW, Zubizarreta,
  estimands". Luedtke and Zubizarreta are co-author surnames, SBW is the
  standard abbreviation for stable balancing weights (spelled out in the
  same sentence), and "estimands" is standard statistical terminology.

The local run additionally reports "unable to verify current time", a
sandbox artifact of the machine the check was run on.
