#' @title Walking Dead 2
#' @description This package simulates the events of infection and healing
#' in a population of susceptible individuals, infected individuals, and healed
#' individuals. The function is only at the development stage and thus has a
#' limited number of arguments that can be altered for the simulation process.
#' @param n,m,k Numeric scalars. \code{n} represents the total population size;
#' \code{m} represents the infected population size; \code{k} represents the
#' healing population size.
#' @param infect numeric value between 0 and 1. This represents the probability
#' of rate of change from one status to another on contact.
#' @param infect_radius numeric value between 0 and 1. Represents the radius within which an
#' individual can experience a status change when exposed
#' @param nIterations positive integer. Represents the number of trials to run
#' where each individual's status and position may be updated.
#' @param flight logical with default set to TRUE. When TRUE, the infected are
#' programmed to repel from the healers. When FALSE, the infected are programmed
#' to be attracted to the non-infected.
#' infected.
#' @examples
#' require(data.table)
#' require(grDevices)
#' temp <- wd2(n=100, m=15, k=5)
#' temp2 <- summary(temp)
#' temp3 <- plot(temp)
#' # head(temp2)
#' # temp3
#' @export
wd2 <-  function(n,
                 m, # infected group
                 k, # curing group
                 infect = .75, # infection rate
                 infect_radius = .075, # radius of infection
                 nIterations = 10, # Number of iterations
                 flight = TRUE # Do infected flee or not
                 ){
  # Group individuals into three categories given respective probs
  statusGrouping <- stats::rmultinom(n, size= 1, prob = c((1-m/n-k/n),m/n,k/n))

  # Creating the baseline environment
  ans <- list2env(list(
    n      = n,
    pos    = matrix(runif(n*2), ncol=2),
    degree = stats::runif(n, 0, 2*pi),
    status   = apply(statusGrouping, 2, which.max),
    infect = infect,
    radii  = infect_radius,
    speed  = stats::runif(n, 0.005, .025)/2
  ))

  # Computing distances
  ans$D <- as.matrix(stats::dist(ans$pos))

  # Update status of individuals in population
  update_status <- function(x) {
    # Moving contagion
    if (sum(x$status==3) < x$n) {
      newcured <- which(
        (apply(x$D[, x$status==3, drop=FALSE], 1, min) < x$radii) &
          (stats::runif(x$n) < x$infect) & (x$status!=1)
      )
      x$status[newcured] <- 3

      #allcured <- x$status==3

      newsick <- which(
        (apply(x$D[, x$status==2, drop=FALSE], 1, min) < x$radii) &
          (stats::runif(x$n) < x$infect) & (x$status!=3)
      )

      x$status[newsick] <- 2
    }
    invisible()
  }

  # Compute a weighting measure
  weighted_avg <- function(x, W) {

    invW <- 1/(W^2 + 1e-5)
    invW %*% x / rowSums(invW)
  }

  # Update location of individuals in population within defined plane
  update_pos <- function(x){

    # Update angle

    if (sum(x$status==3) < x$n) {

      # attraction <- weighted_avg(
      #   x$pos[!x$sick,,drop=FALSE], x$D[x$sick, !x$sick, drop=FALSE]
      # ) - x$pos[x$sick, ,drop=FALSE]


      attractioncured <- weighted_avg(
        x$pos[x$status==2,,drop=FALSE], x$D[x$status==3, x$status==2, drop=FALSE]
      ) - x$pos[x$status==3, ,drop=FALSE]

      # repulsion  <- weighted_avg(
      #   x$pos[x$sick,,drop=FALSE], x$D[!x$sick, x$sick, drop=FALSE]
      # ) - x$pos[!x$sick, ,drop=FALSE]

      repulsionsick  <- weighted_avg(
        x$pos[x$status==2,,drop=FALSE], x$D[x$status==1, x$status==2, drop=FALSE]
      ) - x$pos[x$status==1, ,drop=FALSE]

      if(flight){
        repulsioncured  <- weighted_avg(
          x$pos[x$status==3,,drop=FALSE], x$D[x$status==2, x$status==3, drop=FALSE]
        ) - x$pos[x$status==2, ,drop=FALSE]
        usethis <- atan2(repulsioncured[,2], repulsioncured[,1]) +
          stats::runif(sum(x$status==2, 0, pi/4)) - pi
      }else{
        attractionsick <- weighted_avg(
          x$pos[x$status==1,,drop=FALSE], x$D[x$status==2, x$status==1, drop=FALSE]
        ) - x$pos[x$status==2, ,drop=FALSE]
        usethis <- atan2(attractionsick[,2], attractionsick[,1]) +
          stats::runif(sum(x$status==2, 0, pi/4))
      }

      # Arctan2. We add pi to the healthy indivuduals' angle so that they go in the
      # opposite direction
      # x$degree[x$sick]  <- atan2(attraction[,2], attraction[,1]) +
      #   runif(sum(x$sick, 0, pi/4))
      x$degree[x$status==2]  <- usethis
      x$degree[x$status==3]  <- atan2(attractioncured[,2], attractioncured[,1]) +
        stats::runif(sum(x$status==3, 0, pi/4))
      x$degree[x$status==1]  <- atan2(repulsionsick[,2], repulsionsick[,1]) +
        stats::runif(sum(x$status==1, 0, pi/4)) + pi

    } else {
      x$degree <- stats::runif(x$n, 0, 2*pi)
    }


    # Update position
    x$pos <- x$pos + x$speed*cbind(cos(x$degree), sin(x$degree))

    # Boundaries
    x$pos[x$pos > 1] <- 1
    x$pos[x$pos < 0] <- 0


    x$D <- as.matrix(stats::dist(x$pos))

    invisible()
  }

  # Plot of each individual's status and position at any time t
  plot_process <- function(x) {

    op <- graphics::par(mar=rep(1, 4), mai = rep(0, 4))
    on.exit(par(op))

    graphics::plot.new()

    graphics::plot.window(c(0,1), c(0,1))
    graphics::rect(0,0,1,1, col = "gray70", border = "transparent")
    graphics::points(
      x$pos,
      bg = c("black", "red","green")[x$status],
      pch = 21, col="white", cex=2
    )
  }
  # ------------------------------------------------------------------------------

  # Initialize
  snapshot <- matrix(0, ncol = n, nrow = nIterations+1)
  snapshot[1,] <- ans$status # Status at baseline

  # Setting animation
  grDevices::graphics.off()
  fig <- magick::image_device(300, 300)

  for (i in 1:nIterations) {

    update_status(ans)
    snapshot[i+1,] <- ans$status
    plot_process(ans)
    update_pos(ans)
  }
  while (!is.null(grDevices::dev.list()))
  grDevices::dev.off()

  structure(list(wd2Summary = snapshot,
                 wd2Plot = fig,
                 n=n),
            class = "wd2")
}
#' @rdname wd2
#' @export
#' @param x an object of class \code{wd2}
#' @param ... additional arguments affecting the summary produced.
plot.wd2 <- function(x,...){
  magick::image_animate(x[["wd2Plot"]], fps = 20)
}

#' @rdname wd2
#' @export
#' @param object an object of class \code{wd2}
#' @param ... additional arguments affecting the summary produced.
summary.wd2 <- function(object,...){
  snapshot <- object[["wd2Summary"]]
  n <- object[["n"]]
  nIteration <- dim(snapshot)[1]
  state <- matrix(NA, ncol = 3, nrow = nIteration)
  for(i in 1:nIteration){
    state[i, ] <- c(mean(snapshot[i,]==1),
                    mean(snapshot[i,]==2),
                    mean(snapshot[i,]==3)
    )
  }
  stateDT <- data.table::setDT(as.data.frame(state))[]
  data.table::setnames(stateDT,1:3,c('Non_Infected','Infected','Healed'))
  stateDT[,run:=1:nrow(stateDT)]
  stateDT[,changeNI:= (Non_Infected - c(NA,Non_Infected[-.N]))*n]
  stateDT[,changeIH:= (Infected-- c(NA,Infected[-.N]))*n]
  stateDT[,transitions:= changeNI+changeIH]
  return(stateDT)
}

