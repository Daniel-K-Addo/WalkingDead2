#' @title wd2 - Walking Dead 2
#' @description `wd2` simulates the possible event of transition from one state to another as
#' defined in a population of susceptible, infected, and immunized subjects under defined contexts.
#' Susceptible subjects can transition to be infected; infected subjects can transition to be
#' immunized; immunized subjects do not undergo any transition. As such, equilibrium in the
#' sub-population sizes is achieved when there are no infected subjects.
#' @param n,m,k numeric scalars. \code{n} represents the total population size;
#' \code{m} represents the infected population size; \code{k} represents the immunized
#' population size. Note: \eqn{n \geq 3} and \eqn{{m,k} > 0}.
#' @param fixedAllocation logical with default set to \code{FALSE}. When \code{FALSE},
#' the initial sub population sizes are randomly determined using multinomial distribution
#' given \code{{n,m,k}}. When \code{TRUE}, the initial sub population sizes will be exactly
#' as stated given \code{{n,m,k}}.
#' @param infect numeric value between 0 and 1 with default set at 0.75. This
#' represents the infection rate.
#' @param infect_radius numeric value between 0 and 1 with default set at 0.075.
#' This represents the radius of a transition-causing subject within which an eligible
#' subject will undergo a status change when exposed. The smaller the \code{infect_radius}
#' the smaller the likelihood of transition events by a transition-causing subject.
#' @param nIterations positive integer. Represents the number of runs during which each
#' subject's sub population membership and location will potentlially be updated. This
#' may be viewed as some unit of time.
#' @param flight logical with default set to \code{TRUE}. When \code{TRUE}, the infected are
#' programmed to move away from the immunized population. When \code{FALSE}, the infected are programmed
#' to move towards the susceptible. This represents the radius of a transition-causing
#' subject within which an eligible subject will undergo a status change when exposed.
#' The smaller the infect_radius the smaller the likelihood of transition events by a transition-causing subject.
#' @examples
#' require(data.table)
#' # temp <- wd2(n=100, m=15, k=5)
#' # summary(temp) |> head()
#' # plot(temp)
#' # print(temp)
#' @export
wd2 <-  function(n, # Population size
                 m, # Initial Infected Population size
                 k, # # Initial Immunized Population size
                 fixedAllocation = FALSE, # Allocation type: randomized, if FALSE
                 infect = .75, # infection rate
                 infect_radius = .075, # radius of infection
                 nIterations = 10, # Number of runs
                 flight = FALSE # Infected will flee from immunized if TRUE
                 ){

  # Group subjects into one of susceptible, infected, OR immunized given fixedAllocation argument
  if(fixedAllocation){
    statusGrouping <- c(rep(1,n-m-k), rep(2,m), rep(3,k))
    statusGrouping <- sample(statusGrouping)
  }else{
    statusGrouping <- stats::rmultinom(n, size= 1, prob = c((1-m/n-k/n),m/n,k/n)) |>
      apply(2, which.max)
    if(all(statusGrouping!=3)){
      statusGrouping[length(statusGrouping)] <- 3
    }
  }

  # Creating the baseline environment
  ans <- list2env(list(
    n      = n,
    pos    = matrix(stats::runif(n*2), ncol=2),
    degree = stats::runif(n, 0, 2*pi),
    status   = statusGrouping,
    infect = infect,
    radii  = infect_radius,
    speed  = stats::runif(n, 0.005, .025)/2
  ))

  # Computing distances
  ans$D <- as.matrix(stats::dist(ans$pos))

  # Update status of subjects in population
  update_status <- function(x) {
    # Moving contagion
    if (sum(x$status==2) != 0) {
      # transitions from infected to immunized
      newcured <- which(
        (apply(x$D[, x$status==3, drop=FALSE], 1, min) < x$radii) &
          (stats::runif(x$n) < x$infect) & (x$status!=1)
      )
      x$status[newcured] <- 3

      # transitions from susceptible to infected
      if (sum(x$status==1) != 0){
        newsick <- which(
          (apply(x$D[, x$status==2, drop=FALSE], 1, min) < x$radii) &
            (stats::runif(x$n) < x$infect) & (x$status!=3)
        )

        x$status[newsick] <- 2
      }
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
    # if (sum(x$status==3) < x$n) {
    if (any(x$status==2)) {

      # Attraction of immunized to infected
      attractioncured <- weighted_avg(
        x$pos[x$status==2,,drop=FALSE], x$D[x$status==3, x$status==2, drop=FALSE]
      ) - x$pos[x$status==3, ,drop=FALSE]

      if (any(x$status==1)) {
      # Repulsion of susceptible from infected
      repulsionsick  <- weighted_avg(
        x$pos[x$status==2,,drop=FALSE], x$D[x$status==1, x$status==2, drop=FALSE]
      ) - x$pos[x$status==1, ,drop=FALSE]
      }

      # Repulsion OR attraction behavior of infected given flight argument
      if(flight){
        repulsioncured  <- weighted_avg(
          x$pos[x$status==3,,drop=FALSE], x$D[x$status==2, x$status==3, drop=FALSE]
        ) - x$pos[x$status==2, ,drop=FALSE]
        usethis <- atan2(repulsioncured[,2], repulsioncured[,1]) +
          stats::runif(sum(x$status==2, 0, pi/4)) - pi
      }else if(any(x$status==1)){
        attractionsick <- weighted_avg(
          x$pos[x$status==1,,drop=FALSE], x$D[x$status==2, x$status==1, drop=FALSE]
        ) - x$pos[x$status==2, ,drop=FALSE]
        usethis <- atan2(attractionsick[,2], attractionsick[,1]) +
          stats::runif(sum(x$status==2, 0, pi/4)) + 2*pi
      }

      # Arctan2. We subtract pi from infected subjects'
      # angle so that they go in the opposite direction of the immunized if flight
      # argument is TRUE.
      x$degree[x$status==2]  <- usethis
      x$degree[x$status==3]  <- atan2(attractioncured[,2], attractioncured[,1]) +
        stats::runif(sum(x$status==3, 0, pi/4))
      if(any(x$status==1)){
        x$degree[x$status==1]  <- atan2(repulsionsick[,2], repulsionsick[,1]) +
          stats::runif(sum(x$status==1, 0, pi/4)) + pi
        }

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

  # Plot of each subject's status and position at any time t
  plot_process <-  function(ans) {
    if(length(unique(ans$status))==3){
      # Dataframe for plotting
      plotData <- data.frame(Status=factor(ans$status, labels = c("Susceptible","Infected","Immune")),
                             posX=ans$pos[,1],
                             posY=ans$pos[,2])
      gp <- ggplot2::ggplot(plotData) +
        ggplot2::geom_point(ggplot2::aes(x=posX,y=posY,col=Status), size=2.5) +
        ggplot2::scale_color_manual(values = c("black", "firebrick","chartreuse3"))#,
                                    # labels = c(paste("Susceptible:", mean(ans$status==1)*ans$n),
                                    #            paste("Infected:", mean(ans$status==2)*ans$n),
                                    #            paste("Immune:", mean(ans$status==3)*ans$n)))
    }else if(all(sort(unique(ans$status))==c(1,2))==TRUE){
      # Dataframe for plotting
      plotData <- data.frame(Status=factor(ans$status, labels = c("Susceptible","Infected")),
                             posX=ans$pos[,1],
                             posY=ans$pos[,2])
      gp <- ggplot2::ggplot(plotData) +
        ggplot2::geom_point(ggplot2::aes(x=posX,y=posY,col=Status), size=2.5) +
        ggplot2::scale_color_manual(values = c("black", "firebrick"))#,
                                    # labels = c(paste("Susceptible:", mean(ans$status==1)*ans$n),
                                    #            paste("Infected:", mean(ans$status==2)*ans$n)))
    }else if(all(sort(unique(ans$status))==c(1,3))==TRUE){
      # Dataframe for plotting
      plotData <- data.frame(Status=factor(ans$status, labels = c("Susceptible","Immune")),
                             posX=ans$pos[,1],
                             posY=ans$pos[,2])
      gp <- ggplot2::ggplot(plotData) +
        ggplot2::geom_point(ggplot2::aes(x=posX,y=posY,col=Status), size=2.5) +
        ggplot2::scale_color_manual(values = c("black","chartreuse3"))#,
                                    # labels = c(paste("Susceptible:", mean(ans$status==1)*ans$n),
                                    #            paste("Immune:", mean(ans$status==3)*ans$n)))
    }else if(all(sort(unique(ans$status))==c(2,3))==TRUE){
      # Dataframe for plotting
      plotData <- data.frame(Status=factor(ans$status, labels = c("Infected","Immune")),
                             posX=ans$pos[,1],
                             posY=ans$pos[,2])
      gp <- ggplot2::ggplot(plotData) +
        ggplot2::geom_point(ggplot2::aes(x=posX,y=posY,col=Status), size=2.5) +
        ggplot2::scale_color_manual(values = c("firebrick","chartreuse3"))#,
                                    # labels = c(paste("Infected:", mean(ans$status==2)*ans$n),
                                    #            paste("Immune:", mean(ans$status==3)*ans$n)))
    }else{
      # Dataframe for plotting
      plotData <- data.frame(Status=factor(ans$status, labels = c("Immune")),
                             posX=ans$pos[,1],
                             posY=ans$pos[,2])
      gp <- ggplot2::ggplot(plotData) +
        ggplot2::geom_point(ggplot2::aes(x=posX,y=posY,col=Status), size=2.5) +
        ggplot2::scale_color_manual(values = c("chartreuse3"))#,
                                    # labels = c(paste("Infected:", mean(ans$status==3)*ans$n)))
    }
    gp + #ggplot2::ggtitle(paste("Simulation Run:", i)) +
      ggplot2::xlim(0,1) + ggplot2::ylim(0,1) +
      ggplot2::theme_minimal() +
      ggplot2::theme(panel.grid.major.x = ggplot2::element_blank(),
                     panel.grid.minor.x = ggplot2::element_blank(),
                     panel.grid.minor.y = ggplot2::element_blank(),
                     panel.grid.major.y = ggplot2::element_blank(),
                     axis.text.x=ggplot2::element_blank(),
                     axis.text.y=ggplot2::element_blank(),
                     axis.ticks.x=ggplot2::element_blank(),
                     axis.ticks.y=ggplot2::element_blank(),
                     axis.title.x=ggplot2::element_blank(),
                     axis.title.y=ggplot2::element_blank(),
                     legend.position="bottom",
                     legend.text = ggplot2::element_text(size = 15),
                     plot.title = ggplot2::element_text(hjust = 0.5),
                     panel.background=ggplot2::element_rect(colour="white", fill = "grey90"))

  }
  # ------------------------------------------------------------------------------

  # Initialize
  snapshot <- matrix(0, ncol = n, nrow = nIterations+1)
  snapshot[1,] <- ans$status # Status at baseline

  # Setting animation
  grDevices::graphics.off()
  fig <- magick::image_graph(400, 400)

  # Begin runs
  for (i in 1:nIterations) {

    update_status(ans)
    snapshot[i+1,] <- ans$status
    p <- plot_process(ans)
    print(p)
    update_pos(ans)
  }

  structure(list(wd2Summary = snapshot,
                 wd2Plot = fig,
                 wd2Object = ans),
            class = "wd2")
}

#' @rdname wd2
#' @export
#' @param object an object of class \code{wd2}
#' @param ... not in use currently
summary.wd2 <- function(object,...){
  # Obtain data for summary
  snapshot <- object[["wd2Summary"]]
  n <- object[["wd2Object"]][["n"]]
  nIteration <- dim(snapshot)[1]
  state <- matrix(NA, ncol = 3, nrow = nIteration)
  for(i in 1:nIteration){
    state[i, ] <- c(mean(snapshot[i,]==1),
                    mean(snapshot[i,]==2),
                    mean(snapshot[i,]==3)
    )
  }
  # Set as data.table and compute estimates of interest
  stateDT <- data.table::setDT(as.data.frame(state))[]
  data.table::setnames(stateDT,1:3,c('Susceptible','Infected','Immune'))
  stateDT[,Run:=0:(nrow(stateDT)-1)]
  stateDT[,ChangeSI:= abs(Susceptible - c(NA,Susceptible[-.N]))*n]
  stateDT[,ChangeII:= (Immune - c(NA,Immune[-.N]))*n]
  stateDT[,Transitions:= ChangeSI+ChangeII]
  data.table::setcolorder(stateDT, c("Run",'Susceptible','Infected','Immune',
                         "ChangeSI","ChangeII","Transitions"))
  return(stateDT)
}

#' @rdname wd2
#' @export
#' @param x an object of class \code{wd2}
#' @param ... not in use currently
plot.wd2 <- function(x,...){
  magick::image_animate(x[["wd2Plot"]], fps = 10)
}

#' @rdname wd2
#' @export
#' @param x an object of class \code{wd2}
#' @param ... not in use currently
print.wd2 <- function(x,...){
  snapshot <- x[["wd2Summary"]]
  n <- x[["wd2Object"]][["n"]]
  nIteration <- dim(snapshot)[1]
  state <- matrix(NA, ncol = 3, nrow = nIteration)
  for(i in 1:nIteration){
    state[i, ] <- c(round(mean(snapshot[i,]==1)*n),
                    round(mean(snapshot[i,]==2)*n),
                    round(mean(snapshot[i,]==3)*n)
    )
  }
  final <- state |> utils::tail()
  if(state[dim(state)[1],2]==0){
    this <- min(which(state[,2]==0)) - 1
    conclusion <- paste("The pathogen was eliminated successfully during run", this)
  }else{
    conclusion <- paste("The pathogen was NOT eliminated. Infections are still possible.")
  }
  print(conclusion)
  output <- list(snapshot = final)
  return(output)
}
