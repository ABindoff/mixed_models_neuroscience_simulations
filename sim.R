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

n.sims <- 20000
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
save(sims.meta, file = "sims_meta_20000.RData")

ggplot(sims.meta, aes(x = factor(n.replicates), y = factor(fixed.effects), fill = lmer)) +
  geom_raster(interpolate = F, alpha = 1/10) +
  scale_fill_viridis(discrete = F, begin = 0.05, end = 0.95)

ggplot(sims.meta, aes(x = factor(fixed.effects), y = n.pseudo, fill = factor(lm_type2))) +
  geom_boxplot()

m1 <- gam(as.numeric(lm > .05) ~ s(n.pseudo) + fixed.effects, data = sims.meta)
d <- expand.grid(n.pseudo = seq(2, 18, by = 1),
                 fixed.effects = c(0.1, 0.5, 1, 2))
d$Beta_lm <- predict(m1, d)
p.pseudo.lm <- ggplot(d, aes(x = n.pseudo, y = Beta_lm, colour = factor(fixed.effects))) +
  geom_line() +
  ylim(0,1)

m1 <- gam(as.numeric(lm > .05) ~ s(n.replicates) + fixed.effects, data = sims.meta)
d <- expand.grid(n.replicates = seq(6, 40, by = 1),
                 fixed.effects = c(0.1, 0.5, 1, 2))
d$Beta_lm <- predict(m1, d)
p.rep.lm <- ggplot(d, aes(x = n.replicates, y = Beta_lm, colour = factor(fixed.effects))) +
  geom_line() +
  ylim(0,1)

m1 <- gam(as.numeric(lmer > .05) ~ s(n.pseudo) + fixed.effects, data = sims.meta)
d <- expand.grid(n.pseudo = seq(2, 18, by = 1),
                 fixed.effects = c(0.1, 0.5, 1, 2))
d$Beta_lmer <- predict(m1, d)
p.pseudo.lmer <- ggplot(d, aes(x = n.pseudo, y = Beta_lmer, colour = factor(fixed.effects))) +
  geom_line() +
  ylim(0,1)


m1 <- gam(as.numeric(lmer > .05) ~ s(n.replicates) + fixed.effects, data = sims.meta)
d <- expand.grid(n.replicates = seq(6, 40, by = 1),
                 fixed.effects = c(0.1, 0.5, 1, 2))
d$Beta_lmer <- predict(m1, d)
p.rep.lmer <- ggplot(d, aes(x = n.replicates, y = Beta_lmer, colour = factor(fixed.effects))) +
  geom_line() +
  ylim(0,1)

grid.arrange(p.pseudo.lm, p.pseudo.lmer, p.rep.lm, p.rep.lmer, ncol = 2)

m1 <- gam(as.numeric(lm > .05) ~ s(n.pseudo) + sd.random, data = sims.meta)
d <- expand.grid(n.pseudo = seq(2, 18, by = 1),
                 sd.random = c(0.5, 1, 2))
d$Beta_lm <- predict(m1, d)
p.pseudo.lm <- ggplot(d, aes(x = n.pseudo, y = Beta_lm, colour = factor(sd.random))) +
  geom_line() +
  ylim(0,1)

m1 <- gam(as.numeric(lm > .05) ~ s(n.replicates) + sd.random, data = sims.meta)
d <- expand.grid(n.replicates = seq(6, 40, by = 1),
                 sd.random = c(0.5, 1, 2))
d$Beta_lm <- predict(m1, d)
p.rep.lm <- ggplot(d, aes(x = n.replicates, y = Beta_lm, colour = factor(sd.random))) +
  geom_line() +
  ylim(0,1)

m1 <- gam(as.numeric(lmer > .05) ~ s(n.pseudo) + sd.random, data = sims.meta)
d <- expand.grid(n.pseudo = seq(2, 18, by = 1),
                 sd.random = c(0.5, 1, 2))
d$Beta_lmer <- predict(m1, d)
p.pseudo.lmer <- ggplot(d, aes(x = n.pseudo, y = Beta_lmer, colour = factor(sd.random))) +
  geom_line() +
  ylim(0,1)


m1 <- gam(as.numeric(lmer > .05) ~ s(n.replicates) + sd.random, data = sims.meta)
d <- expand.grid(n.replicates = seq(6, 40, by = 1),
                 sd.random = c(0.5, 1, 2))
d$Beta_lmer <- predict(m1, d)
p.rep.lmer <- ggplot(d, aes(x = n.replicates, y = Beta_lmer, colour = factor(sd.random))) +
  geom_line() +
  ylim(0,1)

grid.arrange(p.pseudo.lm, p.pseudo.lmer, p.rep.lm, p.rep.lmer, ncol = 2)




p.pseudo.lm <- ggplot(sims.meta, aes(x = n.pseudo, y = lm, colour = factor(sd.random))) +
  geom_smooth() +
  ylim(0,1)

p.rep.lm <- ggplot(sims.meta, aes(x = n.replicates, y = lm, colour = factor(sd.random))) +
  geom_smooth() +
  ylim(0,1)

p.pseudo.lmer <- ggplot(sims.meta, aes(x = n.pseudo, y = lmer, colour = factor(sd.random))) +
  geom_smooth() +
  ylim(0,1)

p.rep.lmer <- ggplot(sims.meta, aes(x = n.replicates, y = lmer, colour = factor(sd.random))) +
  geom_smooth() +
  ylim(0,1)

grid.arrange(p.pseudo.lm, p.pseudo.lmer, p.rep.lm, p.rep.lmer, ncol = 2)


# can we model these interactions?
m1 <- glm(as.numeric(lm < .05) ~ n.pseudo*sd.random + fixed.effects*n.replicates,
          family = binomial(link = "logit"), 
          sims.meta)
summary(m1)


sims.meta %>%
  group_by(fixed.effects) %>%
  summarise(type2.lm = mean(lm > .05), type2.lmer = mean(lmer > .05))

sims.meta %>%
  group_by(sd.random) %>%
  summarise(type2.lm = mean(lm > .05), type2.lmer = mean(lmer > .05))


##  Type 1 error rate (no effect)

n.sims <- 10000
n.replicates <- floor(runif(n.sims, 6, 40)/2)*2 # even number between 6 and 40
n.pseudo <- sample(c(2:18), n.sims, replace = TRUE)
sd.error <- 1
random.assignment = TRUE
sd.random <- sample(c(0.5, 1, 2), n.sims, replace = TRUE)
fixed.effects <- 0  # no difference between groups

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


x <- t(sapply(sims, tests))
sims.meta$lm <- x[,1]
sims.meta$lmer <- x[,2]
save(sims.meta, file = "sims_meta_no_effect_10000.RData")

ggplot(sims.meta, aes(x = factor(n.replicates), y = factor(fixed.effects), fill = lmer)) +
  geom_raster(interpolate = F, alpha = 1/10) +
  scale_fill_viridis(discrete = F, begin = 0.05, end = 0.95)

sims.meta$lm_type1 <- sims.meta$lm < .05


p.pseudo.lm <- ggplot(sims.meta, aes(x = n.pseudo, y = lm, colour = factor(fixed.effects))) +
  geom_smooth() +
  ylim(0,1)

p.rep.lm <- ggplot(sims.meta, aes(x = n.replicates, y = lm, colour = factor(fixed.effects))) +
  geom_smooth() +
  ylim(0,1)

p.pseudo.lmer <- ggplot(sims.meta, aes(x = n.pseudo, y = lmer, colour = factor(fixed.effects))) +
  geom_smooth() +
  ylim(0,1)

p.rep.lmer <- ggplot(sims.meta, aes(x = n.replicates, y = lmer, colour = factor(fixed.effects))) +
  geom_smooth() +
  ylim(0,1)

grid.arrange(p.pseudo.lm, p.pseudo.lmer, p.rep.lm, p.rep.lmer, ncol = 2)

p.pseudo.lm <- ggplot(sims.meta, aes(x = n.pseudo, y = lm, colour = factor(sd.random))) +
  geom_smooth() +
  ylim(0,1)

p.rep.lm <- ggplot(sims.meta, aes(x = n.replicates, y = lm, colour = factor(sd.random))) +
  geom_smooth() +
  ylim(0,1)

p.pseudo.lmer <- ggplot(sims.meta, aes(x = n.pseudo, y = lmer, colour = factor(sd.random))) +
  geom_smooth() +
  ylim(0,1)

p.rep.lmer <- ggplot(sims.meta, aes(x = n.replicates, y = lmer, colour = factor(sd.random))) +
  geom_smooth() +
  ylim(0,1)

grid.arrange(p.pseudo.lm, p.pseudo.lmer, p.rep.lm, p.rep.lmer, ncol = 2)

sims.meta %>%
  group_by(fixed.effects) %>%
  summarise(type1.lm = mean(lm < .05), type1.lmer = mean(lmer < .05))

sims.meta %>%
  group_by(sd.random) %>%
  summarise(type1.lm = mean(lm < .05), type1.lmer = mean(lmer < .05))

m1 <- gam(as.numeric(lm < .05) ~ s(n.pseudo) + sd.random, data = sims.meta)
d <- expand.grid(n.pseudo = seq(2, 18, by = 1),
                 sd.random = c(0.5, 1, 2))
d$fit.lm <- predict(m1, d)
p1 <- ggplot(d, aes(x = n.pseudo, y = fit.lm, colour = factor(sd.random))) +
  geom_line() +
  ylim(0,0.6) +
  ylab("Type 1 error rate, lm")

m1 <- gam(as.numeric(lmer < .05) ~ s(n.pseudo) + sd.random, data = sims.meta)
d$fit.lmer <- predict(m1, d)
p2 <- ggplot(d, aes(x = n.pseudo, y = fit.lmer, colour = factor(sd.random))) +
  geom_line() +
  ylim(0,0.6) +
  ylab("Type 1 error rate, lmer")

m1 <- gam(as.numeric(lmer < .05) ~ s(n.replicates) + sd.random, data = sims.meta)
d <- expand.grid(n.replicates = seq(6, 40, by = 1),
                 sd.random = c(0.5, 1, 2))
d$fit.lmer.rep <- predict(m1, d)
p3 <- ggplot(d, aes(x = n.replicates, y = fit.lmer.rep, colour = factor(sd.random))) +
  geom_line() +
  ylim(0,0.6) +
  ylab("Type 1 error rate, lmer")


m1 <- gam(as.numeric(lm < .05) ~ s(n.replicates) + sd.random, data = sims.meta)
d <- expand.grid(n.replicates = seq(6, 40, by = 1),
                 sd.random = c(0.5, 1, 2))
d$fit.lm.rep <- predict(m1, d)
p4 <- ggplot(d, aes(x = n.replicates, y = fit.lm.rep, colour = factor(sd.random))) +
  geom_line() +
  ylim(0,0.6) +
  ylab("Type 1 error rate, lm")

grid.arrange(p1, p2, p4, p3, ncol = 2)
