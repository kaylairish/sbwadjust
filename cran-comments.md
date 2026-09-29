## Submission

First CRAN submission of sbwadjust (0.3.0).

## Test environments

* local macOS, R 4.3.2: R CMD check --as-cran (with CRAN incoming checks)
* R-hub (R-devel): linux, windows, macos-arm64 -- all OK
* GitHub Actions: macOS, Windows (R release); Ubuntu (R devel, release,
  oldrel-1) -- all OK

## R CMD check results

0 errors | 0 warnings | 1 note

The note is the CRAN incoming feasibility check:

* "New submission" -- this is the package's first release.

* It may also list "Possibly misspelled words in DESCRIPTION: Luedtke, SBW,
  Zubizarreta, estimands". Luedtke and Zubizarreta are co-author surnames,
  SBW is the standard abbreviation for stable balancing weights (spelled out
  in the same sentence), and "estimands" is standard statistical terminology.

The local run additionally reports "unable to verify current time" and HTML
Tidy warnings on the HTML manual; both come from the local machine (no
time-server access; macOS's bundled 2006 HTML Tidy) and do not appear on
R-hub.
