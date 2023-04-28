
<!-- README.md is generated from README.Rmd. Please edit that file -->

# WalkingDead2

<!-- badges: start -->

[![R-CMD-check](https://github.com/UofUEpiBio/WalkingDead2/workflows/R-CMD-check/badge.svg)](https://github.com/UofUEpiBio/WalkingDead2/actions)
<!-- badges: end -->

The goal of WalkingDead2 is to simulate the events of infection (of
susceptible subjects) and immunization (of infected subjects) in a
population of susceptible, infected, and immunized individuals. It is an
extension of the
[**WalkingDead**](https://github.com/UofUEpiBio/PHS7045-advanced-programming/tree/main/projects/10-walking-dead)
project by [George Vegayon](https://github.com/gvegayon).

## Installation

**WalkingDead2** is available through the author’s Github ([click
here](https://github.com/Daniel-K-Addo)). In R, you can install and load
the package using the following code:

``` r
# Installing WalkingDead2 from github 
devtools::install_github("Daniel-K-Addo/WalkingDead2")
```

## Brief Details

Currently, **WalkingDead2** accommodates the transition of subjects from
the susceptible to infected state upon contact with an infected subject,
as well as the transition of subjects from the infected to immunized
state upon contact with an immunized subject. Immunized subjects do not
undergo any transitions. This is practical for a communicable disease
for which there is a cure that immunizes the infected but no vaccine for
the susceptible. This scenario can very well be extended to different
fields.

## Functions and Arguments

There are eight arguments in, `wd2`, the primary function of
**WalkingDead2**. Five of these have default values. Briefly, these
arguments perform the following roles:

- `n, m, k` `n` represents the total population size; `m` represents the
  infected population size; `k` represents the immunized population
  size. Note: $n \geq 3$ and ${m,k} > 0$.
- `fixedAllocation` a logical which determines whether the initial sub
  population sizes are randomly determined using multinomial
  distribution OR fixed given `{n,m,k}`.
- `infect` numeric value between $0$ and $1$ representing the infection
  rate.
- `infect_radius` numeric value between $0$ and $1$ representing the
  radius of a transition-causing subject within which an eligible
  subject will undergo a status change when exposed.
- `nIterations` represents a number of runs during which each subject’s
  sub population membership and location will potentially be updated.
  This may be viewed as some unit of time.
- `flight` a logical which determines if the infected population are
  programmed to move away from the immunized population (flight) OR move
  towards the susceptible (fighting to keep population numbers up). By
  default, the susceptible population will always move away from the
  infected population whereas the immunized population will always move
  towards the infected population.

The output from `wd2` has class $wd2$ and may be passed through `print`,
`summary`, `plot`, and `snapshot` functions for varied results as
demonstrated in the example below.

## Example

This is a basic example which shows you how to solve a common problem:

``` r
# Load the package
library(WalkingDead2)
set.seed(2608)
# Primary function
temp <- wd2(n = 100, # Total Population size 
            m = 16, # Initial Infected Population size
            k = 4, # Initial Immunized Population size
            fixedAllocation = FALSE, # Allocation type: Randomize if FALSE
            infect = 0.75, 
            infect_radius = 0.075, # Transition-causing radius
            nIterations = 100, # Total number of runs
            flight = FALSE # Infected will flee from immunized if TRUE
            )

# Demonstrate print.wd2
print(temp) # or print.wd2(temp)
#> [1] "The pathogen was eliminated successfully after run 77"
#> $snapshot
#>        [,1] [,2] [,3]
#>  [96,]   33    0   67
#>  [97,]   33    0   67
#>  [98,]   33    0   67
#>  [99,]   33    0   67
#> [100,]   33    0   67
#> [101,]   33    0   67

# Demonstrate summary.wd2 and view first and last few rows
summary(temp) |> head(10) # or summary.wd2(temp)
#>     Run Susceptible Infected Immune ChangeSI ChangeII Transitions
#>  1:   0        0.83     0.12   0.05       NA       NA          NA
#>  2:   1        0.73     0.20   0.07       10        2          12
#>  3:   2        0.68     0.25   0.07        5        0           5
#>  4:   3        0.65     0.27   0.08        3        1           4
#>  5:   4        0.62     0.29   0.09        3        1           4
#>  6:   5        0.60     0.31   0.09        2        0           2
#>  7:   6        0.60     0.31   0.09        0        0           0
#>  8:   7        0.60     0.31   0.09        0        0           0
#>  9:   8        0.60     0.31   0.09        0        0           0
#> 10:   9        0.59     0.32   0.09        1        0           1
summary(temp) |> tail(5)
#>    Run Susceptible Infected Immune ChangeSI ChangeII Transitions
#> 1:  96        0.33        0   0.67        0        0           0
#> 2:  97        0.33        0   0.67        0        0           0
#> 3:  98        0.33        0   0.67        0        0           0
#> 4:  99        0.33        0   0.67        0        0           0
#> 5: 100        0.33        0   0.67        0        0           0

# Demonstrate plot.wd2
plot(temp) # or plot.wd2(temp)
```

<img src="man/figures/README-example-1.gif" width="100%" />

``` r

# Demonstrate snapshot
snapshot(temp)
#>   SampleSize InfectionRate Run PropSusceptible PropInfected PropRecovered
#> 1        100          0.75   1               1            0             0
#>   Equilibrium
#> 1        TRUE
```

    #> [1] "The pathogen was eliminated successfully after run 77"

This second plot shows the case when flight is set to TRUE.

<img src="man/figures/README-unnamed-chunk-3-1.gif" width="100%" />

    #> [1] "The pathogen was NOT eliminated. Infections are still possible."
