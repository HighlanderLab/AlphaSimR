How I set up an R package using the tskit C API:

1. Create a basic R package structure (in an existing R project folder)

I removed other files and created a package directory with:
```
AlphaSimRTmp/
  DESCRIPTION
  R/
  src/
```

(R will not recognise the directory as a package unless DESCRIPTION, R/, and src/ all exist)

2. Write a minimal DESCRIPTION

I added a minimal DESCRIPTION file.

```
Package: AlphaSimRTmp
Type: Package
Version: 0.0.1
Imports: Rcpp
LinkingTo: Rcpp
```

3. Vendor tskit and kastore (only C API files)

Tskit and karstore are from: https://github.com/tskit-dev/tskit/archive/refs/tags/1.0.0.zip
Inside src/, I created a deps/ directory and copied in only the C API parts:

```
src/deps/
  tskit/
  kastore/
  tskit.h
```
tskit folder from: tskit-1.0.0/c/tskit
karstore folder from: tskit-1.0.0/c/subprojects/kastore
tskit.h from: tskit-1.0.0/c/tskit.h

(I did not copy meson.build, examples, Python code and documentation etc.)

4. Create a minimal C++ test file

I added src/minimal.cpp with:

#include <tskit.h>

two exported Rcpp functions:

one to report the tskit version

one smoke test that loads a .trees file using
tsk_table_collection_init → load → tsk_treeseq_init

5. Create an initial NAMESPACE so Rcpp::compileAttributes() could run

Before anything would compile, I created a minimal NAMESPACE in R:

```
writeLines(c(
  "useDynLib(AlphaSimRTmp, .registration=TRUE)",
  "importFrom(Rcpp, evalCpp)"
), "NAMESPACE")
```

6. Write Makevars to compile vendored C code

In src/Makevars, I:

added include paths for deps, deps/tskit, and deps/kastore;

explicitly listed all tskit and kastore .c files;

added custom rules to compile .c files in subdirectories;

linked the resulting .o files into the package shared library.

7. Generate Rcpp and roxygen outputs

From the package root, in R:
```
Rcpp::compileAttributes()
roxygen2::roxygenise()
```
After this, R/RcppExports.R, src/RcppExports.cpp were created, but the minimal handwritten NAMESPACE was not replaced.

8.  Add zzz.R so roxygen2 can generate a correct NAMESPACE

Before generating documentation, I created R/zzz.R with the following contents:
```
#' @useDynLib AlphaSimRTmp, .registration = TRUE
#' @importFrom Rcpp evalCpp
NULL
```

This ensures that roxygen2 writes the required useDynLib() and importFrom(Rcpp, evalCpp) entries into NAMESPACE.

9. I remove the handwritten NAMESPACE and generated it with running `roxygen2::roxygenise()` again.

10. Clean install the package (not necessary; just because of bugs in lazy-load database caused by repeated installs)

```
rm -rf ~/Library/R/arm64/4.5/library/AlphaSimRTmp
rm -rf ~/Library/R/arm64/4.5/library/00LOCK-AlphaSimRTmp
R CMD INSTALL --preclean AlphaSimRTmp
```

11. Test

In a fresh R session:
```
> library(AlphaSimRTmp)
> tskit_version_test()
major minor patch 
    1     3     0 
> 
> tskit_smoke_load_free("...Projects/test_TSK2ASR/data/simulations/normal/msprime_chr1.trees")
[1] 1
```
