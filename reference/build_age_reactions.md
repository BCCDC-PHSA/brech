# get reaction list for an age-structured SEIR model

Construct reactions from number of age groups m and an m x m
contact_matrix for an SEIR model with case detection, hospitalizations,
and pre-exposure prophylaxis.

## Usage

``` r
build_age_reactions(contact_matrix)
```

## Arguments

- contact_matrix:

  matrix m x m

## Value

list of reactions
