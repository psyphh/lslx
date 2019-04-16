## lslx ver.0.6.6

## Test environments
* x86_64-apple-darwin15.6.0 (64-bit), R 3.5.0
* x86_64-w64-mingw32 (64-bit) via devtools::build_win()

## R CMD check results
no ERRORs, WARNINGs or NOTEs

## devtools::build_win() results
no ERRORs or WARNINGs, but 2 NOTEs
* checking installed package size ... NOTE
  installed size is  5.3Mb
  sub-directories of 1Mb or more:
    doc    2.5Mb
    libs   2.4Mb
* checking sizes of PDF files under 'inst/doc' ... NOTE
  'qpdf' made some significant size reductions:
     compacted 'vignette-lslx.pdf' from 744Kb to 639Kb
  consider running tools::compactPDF() on these files
  'gs+qpdf' made some significant size reductions:
     compacted 'vignette-lslx.pdf' from 639Kb to 333Kb
