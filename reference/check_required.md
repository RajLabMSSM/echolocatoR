# Check required arguments

Check whether any required arguments to a function are missing. An
argument is considered required if there is no default.

## Usage

``` r
check_required(
  fun = echolocatoR::finemap_loci,
  args = match.call(call = sys.call(sys.parent(2)), expand.dots = FALSE),
  return_no_default = FALSE
)
```

## Arguments

- fun:

  Function to check.

- args:

  Argument calls to assess.

## Value

Null
