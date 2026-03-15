# Argument types

Categorize arguments of a given function by type.

## Usage

``` r
arg_types(
  fun,
  keep_types = c("value", "function", "NULL", "no_default"),
  verbose = FALSE
)
```

## Arguments

- fun:

  Function to evaluate.

- keep_types:

  Argument types to keep in output.

- verbose:

  Print messages.

## Value

Named list of argument types.
