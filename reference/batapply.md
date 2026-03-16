# echolocatoR-themed progress bar

Iterator function with echolocatoR-themed progress bar.

## Usage

``` r
batapply(
  X,
  FUN = function(l) Sys.sleep(5/100),
  apply_func = lapply,
  total = length(X),
  name = NULL,
  .envir = parent.frame(),
  cli.progress_show_after = 0,
  clear = FALSE,
  color1 = cli::col_br_cyan,
  color2 = cli::col_br_magenta,
  ...
)
```

## Arguments

- X:

  a vector (atomic or list) or an
  [`expression`](https://rdrr.io/r/base/expression.html) object. Other
  objects (including classed objects) will be coerced by
  `base::`[`as.list`](https://rdrr.io/r/base/list.html).

- FUN:

  the function to be applied to each element of `X`: see ‘Details’. In
  the case of functions like `+`, `%*%`, the function name must be
  backquoted or quoted.

- apply_func:

  Iterator function to use (default:
  [lapply](https://rdrr.io/r/base/lapply.html)).

- total:

  Passed to
  [`cli_progress_bar()`](https://cli.r-lib.org/reference/cli_progress_bar.html).

- name:

  Name of the progress bar, a label, passed to
  [`cli_progress_bar()`](https://cli.r-lib.org/reference/cli_progress_bar.html).

- .envir:

  Passed to
  [`cli_progress_bar()`](https://cli.r-lib.org/reference/cli_progress_bar.html).

- cli.progress_show_after:

  How long to wait before showing the progress bar.

- clear:

  Whether to remove the progress bar from the screen after it has
  terminated. Defaults to the `cli.progress_clear` option, or `TRUE` if
  unset.

- color1:

  First color to use in the palette.

- color2:

  Second color to use in the palette.

- ...:

  Additional arguments passed to `apply_func`.

## Value

A (named) list.

## Examples

``` r
out <- batapply(X = seq_len(30))
#> ⠙   [  0%] 🧬 ))>))>))>))>))>))>))>))>))>))>))>))>))>))>))>))>))>))>))>))>))>))…
#> ⠹   [100%] 🧬  🦇
#> 
```
