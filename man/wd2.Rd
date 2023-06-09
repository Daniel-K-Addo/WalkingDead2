% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wd2Function.R
\name{wd2}
\alias{wd2}
\alias{summary.wd2}
\alias{plot.wd2}
\alias{print.wd2}
\alias{snapshot}
\title{wd2 - Walking Dead 2}
\usage{
wd2(
  n,
  m,
  k,
  fixedAllocation = FALSE,
  infect = 0.75,
  infect_radius = 0.075,
  nIterations = 10,
  flight = FALSE
)

\method{summary}{wd2}(object, ...)

\method{plot}{wd2}(x, ...)

\method{print}{wd2}(x, ...)

snapshot(x)
}
\arguments{
\item{n, m, k}{numeric scalars. \code{n} represents the total population size;
\code{m} represents the infected population size; \code{k} represents the immunized
population size. Note: \eqn{n \geq 3} and \eqn{{m,k} > 0}.}

\item{fixedAllocation}{logical with default set to \code{FALSE}. When \code{FALSE},
the initial sub population sizes are randomly determined using multinomial distribution
given \code{{n,m,k}}. When \code{TRUE}, the initial sub population sizes will be exactly
as stated given \code{{n,m,k}}.}

\item{infect}{numeric value between 0 and 1 with default set at 0.75. This
represents the infection rate.}

\item{infect_radius}{numeric value between 0 and 1 with default set at 0.075.
This represents the radius of a transition-causing subject within which an eligible
subject will undergo a status change when exposed. The smaller the \code{infect_radius}
the smaller the likelihood of transition events by a transition-causing subject.}

\item{nIterations}{positive integer. Represents the number of runs during which each
subject's sub population membership and location will potentlially be updated. This
may be viewed as some unit of time.}

\item{flight}{logical with default set to \code{TRUE}. When \code{TRUE}, the infected are
programmed to move away from the immunized population. When \code{FALSE}, the infected are programmed
to move towards the susceptible. This represents the radius of a transition-causing
subject within which an eligible subject will undergo a status change when exposed.
The smaller the infect_radius the smaller the likelihood of transition events by a transition-causing subject.}

\item{object}{an object of class \code{wd2}}

\item{...}{not in use currently}

\item{x}{an object of class \code{wd2}}
}
\description{
\code{wd2} simulates the possible event of transition from one state to another as
defined in a population of susceptible, infected, and immunized subjects under defined contexts.
Susceptible subjects can transition to be infected; infected subjects can transition to be
immunized; immunized subjects do not undergo any transition. As such, equilibrium in the
sub-population sizes is achieved when there are no infected subjects.
}
\examples{
require(data.table)
temp <- wd2(n=100, m=15, k=5)
summary(temp) |> head()
plot(temp)
print(temp)
snapshot(temp)
}
