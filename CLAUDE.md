# mismeasured

R package for bias correction in GLMs with measurement error and misclassification.

## References policy

**All academic references MUST be verified before committing.** Run `/verify-refs` on any file that contains citations (DESCRIPTION, README.Rmd, R/*.R roxygen blocks, vignettes). Every reference must have:

- Correct author names, year, title, journal
- A valid DOI (or arXiv ID for preprints)
- Consistency between the DOI metadata and what's written in the file

## Build and test

```bash
R CMD INSTALL .
Rscript -e 'devtools::test()'
```
