library(magrittr)
library(dplyr)
library(ggplot2)
library(lme4)
library(gridExtra)
library(lmerTest)
library(viridis)
library(ggbeeswarm)
library(purrr)

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
sim.mm <- function(...,
                   n.replicates = 4,
                   sd.random = 0.5,
                   fixed.effects = 0,
                   n.pseudo = 3,
                   random.assignment = T,
                   sd.error = NULL){
  animals <- data.frame(animal.id = as.factor(seq(n.replicates)),
                        random.effect = rnorm(n.replicates, 0, sd.random))
  conditions <- data.frame(condition = c("contr", "exp"),
                           fixed.effect = c(0,fixed.effects))
  pseudoreplicate <- data.frame(pseudoreplicate = c(1:n.pseudo))
  d <- expand.grid.df(pseudoreplicate, animals, conditions)
  if(random.assignment){
    i <- sample(unique(d$animal.id), n.replicates/2, replace = FALSE)
    d %<>% filter(xor(animal.id %nin% i, condition == "exp"))
  }
  d %<>% mutate(y.hat = random.effect + fixed.effect)
  if(!is.null(sd.error)){
    d %<>% mutate(y = y.hat + rnorm(nrow(d), 0, sd.error)) %>%
      group_by(animal.id, condition) %>%
      mutate(y.mean = mean(y))
  }
  d
}

# example
plot.sim(sim.mm(n.replicates = 8, n.pseudo = 5, sd.error = 1/2, sd.random = 1, fixed.effects = 2))

# now we have a method of simulating data with theta parameters (and a method of plotting should we need to)
# we can explore the parameter space and see what effect modelling approach has on Type 1 and Type 2 error
# rates

# consider the case of the randomized controlled trial to begin with
# for each trial, choose quasi-random parameters
# estimate Type II error rate over a range of effect sizes

n.sims <- 500
n.replicates <- floor(runif(n.sims, 6, 40)/2)*2 # even number between 6 and 40
n.pseudo <- sample(c(2:18), n.sims, replace = T)
sd.error <- 1
random.assignment = TRUE
sd.random <- sample(c(0.5, 1, 2), n.sims, replace = T)
fixed.effects <- sample(c(0.1, 0.5, 1, 2), n.sims, replace = T)

sims.meta <- data.frame(n.replicates = n.replicates,
                   n.pseudo = n.pseudo,
                   sd.error = sd.error,
                   random.assignment = random.assignment,
                   sd.random = sd.random,
                   fixed.effects = fixed.effects)

sims <- pmap(sims.meta, sim.mm)  # simulates a data-set using parameters from each row of sims.meta

# let's plot the first of these simulations and compare to the parameters as a sanity check
plot.sim(sims[[1L]]) 
sims.meta[1,]

tests <- function(x){
  m1 <- anova(lm(y ~ condition, data = x))
  m2 <- anova(lmer(y ~ condition + (1|animal.id), data = x))
  return(cbind(m1$`Pr(>F)`[1], m2$`Pr(>F)`))
}

x <- t(sapply(sims, tests))
sims.meta$lm <- x[,1]
sims.meta$lmer <- x[,2]
sims.meta$lm_type2 <- sims.meta$lm >= .05

ggplot(sims.meta, aes(x = factor(n.replicates), y = factor(n.pseudo), fill = lm)) +
  geom_raster(interpolate = T, alpha = 4/5) +
  scale_fill_viridis(discrete = F, begin = 0.05, end = 0.95)

ggplot(sims.meta, aes(x = factor(fixed.effects), y = n.pseudo, fill = factor(lm_type2))) +
  geom_boxplot()

p.pseudo.lm <- ggplot(sims.meta, aes(x = n.pseudo, y = lm, colour = factor(fixed.effects))) +
  geom_smooth(method = "lm") 

p.rep.lm <- ggplot(sims.meta, aes(x = n.replicates, y = lm, colour = factor(fixed.effects))) +
  geom_smooth(method = "lm") 

p.pseudo.lmer <- ggplot(sims.meta, aes(x = n.pseudo, y = lmer, colour = factor(fixed.effects))) +
  geom_smooth(method = "lm") 

p.rep.lmer <- ggplot(sims.meta, aes(x = n.replicates, y = lmer, colour = factor(fixed.effects))) +
  geom_smooth(method = "lm") 

grid.arrange(p.pseudo.lm, p.pseudo.lmer, p.rep.lm, p.rep.lmer, ncol = 2)

p.pseudo.lm <- ggplot(sims.meta, aes(x = n.pseudo, y = lm, colour = factor(sd.random))) +
  geom_smooth(method = "lm") 

p.rep.lm <- ggplot(sims.meta, aes(x = n.replicates, y = lm, colour = factor(sd.random))) +
  geom_smooth(method = "lm") 

p.pseudo.lmer <- ggplot(sims.meta, aes(x = n.pseudo, y = lmer, colour = factor(sd.random))) +
  geom_smooth(method = "lm") 

p.rep.lmer <- ggplot(sims.meta, aes(x = n.replicates, y = lmer, colour = factor(sd.random))) +
  geom_smooth(method = "lm") 

grid.arrange(p.pseudo.lm, p.pseudo.lmer, p.rep.lm, p.rep.lmer, ncol = 2)
