# Argument list handler

Infer whether the argument should be applied to all loci or matched by
index.

## Usage

``` r
arg_list_handler(arg, i, loci, use_names = FALSE, error = TRUE)
```

## Arguments

- arg:

  Name of the argument to evaluate.

- i:

  Iterator index.

- loci:

  A vector of loci that are being iterated over.

- use_names:

  Whether to identify locus-specific value by name instead of index.

- error:

  Throw an error (`TRUE`) instead of a warning (`FALSE`).

## Value

Selected argument value.
