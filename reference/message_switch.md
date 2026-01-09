# switch message output based in input

switch message output based in input

## Usage

``` r
message_switch(type, verbose = TRUE, ...)
```

## Arguments

- type:

  name of message output to produce

- verbose:

  logical if message should be outputed

- ...:

  any other parameter

## Details

Essentially, all message types are listed in this giant switch. The type
defines what message you want to output. We select named arguments from
\`...\` and parse them when required. Adding message is straight forward
by simply using the template in other switch types. This also means that
you can add and remove messages without causing errors. While this might
lead to dead code, better dead code than buggy code with undefined
functions...
