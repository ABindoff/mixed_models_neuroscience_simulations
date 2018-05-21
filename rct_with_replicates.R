library(magrittr)
library(dplyr)
library(ggplot2)
library(lme4)
library(gridExtra)
library(lmerTest)
library(viridis)
library(ggbeeswarm)
ggplot <- function(...) ggplot2::ggplot(...) + scale_colour_viridis(discrete = T, begin = 0.05, end = 0.8) + scale_fill_viridis(discrete = T, end = 0.5)
options(digits = 5)
set.seed(1)


plot.sim <- function(d){
  ggplot(d, aes(x = condition, y = y, colour = animal.id, fill = animal.id)) +
    geom_quasirandom(dodge.width = 1, alpha = 1/2) +
    geom_point(aes(y = y.mean, group = animal.id), colour = "black", size = 4, position = position_dodge(width = 1), shape = 3) +
    ylim(-3, 4)
} 


expand.grid.df <- function(...) Reduce(function(...) merge(..., by=NULL), list(...))
`%nin%` <- Negate(`%in%`)

animals <- data.frame(animal.id = c("01", "02", "03", "04", "05", "06"),
                      random.effect = rnorm(6, 0, 0.5))
conditions <- data.frame(condition = c("exp", "contr"),
                         fixed.effect = c(0, 0))


replicates <- data.frame(replicate = c(1:20))
d <- expand.grid.df(replicates, animals, conditions)
d %<>% filter(xor(animal.id %nin% c("01", "02", "03"), condition == "exp")) %>% # balanced RCT
  mutate(y.hat = random.effect + fixed.effect)

d1 <- d %>% mutate(y = y.hat + rnorm(nrow(d), 0, 1)) %>% group_by(animal.id, condition) %>% mutate(y.mean = mean(y))
plot.sim(d1)

p1 <- list(NULL)
p2 <- list(NULL)
p3 <- list(NULL)
p4 <- list(NULL)

for(j in 1:10){
  replicates <- data.frame(replicate = c(1:(j*9-7)))
  d <- expand.grid.df(replicates, animals, conditions)
  d %<>% filter(xor(animal.id %nin% c("01", "02", "03"), condition == "exp")) %>% # balanced RCT
    mutate(y.hat = random.effect + fixed.effect)
  p1[[j]] <- list()
  p2[[j]] <- list()
  p3[[j]] <- list()
  p4[[j]] <- list()
  
  for(i in 1:100){
    d1 <- d %>% mutate(y = y.hat + rnorm(nrow(d), 0, 1))
    a1 <- anova(lm(y ~ condition, data = d1))
    #a2 <- anova(lm(y ~ condition/animal.id, data = d1))
    #a3 <- anova(lmer(y ~ condition + (1|animal.id), data = d1), ddf = "Kenward-Roger")
    a4 <- anova(lmer(y ~ condition + (1|animal.id), data = d1), ddf = "Satterthwaite")
    
    p1[[j]] <- c(p1[[j]], a1$`Pr(>F)`[1])
    #p2[[j]] <- c(p2[[j]], a2$`Pr(>F)`[1])
    #p3[[j]] <- c(p3[[j]], a3$`Pr(>F)`)
    p4[[j]] <- c(p4[[j]], a4$`Pr(>F)`)
  }
}

lapply(p1, function(x) mean(x < .05))
lapply(p4, function(x) mean(x < .05))

hist(unlist(p1), breaks = 20)
hist(unlist(p4), breaks = 20)
