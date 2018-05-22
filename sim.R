library(magrittr)
library(dplyr)
library(ggplot2)
library(lme4)
library(gridExtra)
library(lmerTest)
library(viridis)
library(ggbeeswarm)

# choose a colour-blind friendly pallette 
ggplot <- function(...) ggplot2::ggplot(...) + scale_colour_viridis(discrete = T, begin = 0.05, end = 0.95) + scale_fill_viridis(discrete = T, end = 0.5)

# for reproducibility
set.seed(1)

# generic plot method
plot.sim <- function(d){
  ggplot(d, aes(x = condition, y = y, colour = animal.id, fill = animal.id)) +
    geom_quasirandom(dodge.width = 1, alpha = 1/2, size = 2) +
    geom_point(aes(y = y.mean, group = animal.id), colour = "black", size = 4, position = position_dodge(width = 1), shape = 3)
} 

# `expand.grid.df` works like `expand.grid` but enables the user to define paired columns
# it will return every combination of each argument
expand.grid.df <- function(...) Reduce(function(...) merge(..., by=NULL), list(...))

# `%nin%` "not in", e.g c("a", "b", "c") %nin% c("a", "c") will return "b"
`%nin%` <- Negate(`%in%`)

# if sd.error is specified, `sim.mm` will simulate data with 
# fixed and random effects and their parameters, otherwise
# a data.frame is returned with y-hat but no residuals
# if random.assignment = F, will simulate data from a nested (split-plot) design
# so there will be twice as many observations (need to consider this if estimating
# the power of two different designs)
sim.mm <- function(n.replicates = 4,
                   sd.random = 0.5,
                   fixed.effects = c(0, 0),
                   n.pseudo = 3,
                   random.assignment = T,
                   sd.error = NULL){
  animals <- data.frame(animal.id = as.factor(seq(n.replicates)),
                        random.effect = rnorm(n.replicates, 0, sd.random))
  conditions <- data.frame(condition = c("contr", "exp"),
                           fixed.effect = fixed.effects)
  pseudoreplicate <- data.frame(pseudoreplicate = c(1:n.pseudo))
  d <- expand.grid.df(pseudoreplicate, animals, conditions)
  if(random.assignment){
    i <- sample(unique(d$animal.id), floor(n.replicates/2), replace = FALSE)
    d %<>% filter(xor(animal.id %nin% i, condition == "exp"))
  }
  d %<>% mutate(y.hat = random.effect + fixed.effect)
  if(!is.null(sd.error)){
    d %<>% mutate(y = scale(y.hat + rnorm(nrow(d), 0, sd.error))) %>%
      group_by(animal.id, condition) %>%
      mutate(y.mean = mean(y))
  }
  d
}

# example
plot.sim(sim.mm(n.replicates = 8, n.pseudo = 5, sd.error = 1/2, sd.random = 1, fixed.effects = c(0, 3)))
