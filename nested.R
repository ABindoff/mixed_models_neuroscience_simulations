library(magrittr)
library(dplyr)
library(ggplot2)
library(lme4)
library(gridExtra)
library(lmerTest)
library(viridis)
library(ggbeeswarm)
ggplot <- function(...) ggplot2::ggplot(...) + scale_colour_viridis(discrete = T, end = 0.5) + scale_fill_viridis(discrete = T, end = 0.5)
options(digits = 3)
set.seed(1)

sim.mm <- function(replicates = 6, treatment = 0.5, sd1 = 0.5, sd2 = 1){
  replicates <- 2*(floor(replicates/2))
  replicate1 <- c(1:floor(replicates/2))
  replicate2 <- c((floor(replicates/2)+1):replicates) 
  d1 <- expand.grid(animalid  = c(-0.5, 0.1, 0.4),
                    treatment = c(0, treatment),
                    replicate = replicate1)
  d2 <- expand.grid(animalid  = c(-0.5, 0.1, 0.4),
                    treatment = c(0, treatment),
                    replicate = replicate2)
  # smaller effect sizes (random and fixed)
  small.noise <- abs(rnorm(nrow(d1), 0, sd1))
  d1$y1 <- d1$treatment + d1$animalid + small.noise
  d2$y1 <- d2$treatment + d2$animalid - small.noise
  # larger effect sizes (random and fixed)
  large.noise <- abs(rnorm(nrow(d1), 0, sd2))
  d1$y2 <- d1$treatment + d1$animalid + large.noise
  d2$y2 <- d2$treatment + d2$animalid - large.noise
  
  d <- bind_rows(d1, d2)
  
  
  d$animalid %<>% factor(labels = c("01", "02", "03"))
  d$treatment %<>% factor(labels = c("a", "b"))
  d %<>% group_by(treatment, animalid) %>% mutate(y1mean = mean(y1), y2mean = mean(y2))
  d
}




plot.sim <- function(d){
  p1 <- ggplot(d, aes(x = treatment, y = y1, colour = animalid, fill = animalid)) +
    geom_quasirandom(dodge.width = 1, alpha = 1/3) +
    geom_point(aes(y = y1mean), size = 4, position = position_dodge(width = 1), shape = 3) +
    ylim(-3, 3)
  p2 <- ggplot(d, aes(x = treatment, y = y2, colour = animalid)) +
    geom_quasirandom(dodge.width = 1, alpha = 1/3) +
    geom_point(aes(y = y2mean), size = 4, position = position_dodge(width = 1), shape = 3) +
    ylim(-3, 3)

  grid.arrange(p1, p2, ncol = 2)
}

d <- sim.mm()
plot.sim(d)

# which one would you trust more?


# recover fixed and random parameters using `lmer`
# y1
m1 <- lmer(y1 ~ treatment + (treatment|animalid), data = d)
ranef(m1)
fixef(m1)  # so far, so good

# y2
m2 <- lmer(y2 ~ treatment + (treatment|animalid), data = d)
ranef(m2)
fixef(m2)  # obviously not as good as y1, but not terrible

anova(m1, ddf = "Kenward-Roger")  # justifiably reject null
anova(m2, ddf = "Kenward-Roger")  

# what happens if we don't account for random effect?
# y1
m3 <- lm(y1 ~ treatment, data = d)
coef(m3)

# y2
m4 <- lm(y2 ~ treatment, data = d)
coef(m4)

anova(m3)
anova(m4)


# what happens if we increase the number of replicates?
set.seed(1)

d <- sim.mm(40)

plot.sim(d)



# recover fixed and random parameters using `lmer`
# y1
m1 <- lmer(y1 ~ treatment + (treatment|animalid), data = d)
ranef(m1)
fixef(m1)  # so far, so good

# y2
m2 <- lmer(y2 ~ treatment + (treatment|animalid), data = d)
ranef(m2)
fixef(m2)  

anova(m1, ddf = "Kenward-Roger")  # justifiably reject null
anova(m2, ddf = "Kenward-Roger")  


m3 <- lm(y1 ~ treatment, data = d)
anova(m3)
m4 <- lm(y2 ~ treatment, data = d)
anova(m4)


# we can see that m3 and m4 overestimate degrees of freedom
# not such a problem when we KNOW that there is an effect, right?
# but what if we KNOW that the effect is really, really tiny,
# like the sort of effect that is more likely to do with protocol
# than actual treatment?

d <- sim.mm(replicates = 10, treatment = 0.01)
plot.sim(d)

# recover fixed and random parameters using `lmer`
# y1
m1 <- lmer(y1 ~ treatment + (treatment|animalid), data = d)
ranef(m1)
fixef(m1)  # so far, so good

# y2
m2 <- lmer(y2 ~ treatment + (treatment|animalid), data = d)
ranef(m2)
fixef(m2)  

anova(m1, ddf = "Kenward-Roger")  
anova(m2, ddf = "Kenward-Roger")  

# our estimates are very good, and arguably, we correctly reject the null

# repeat the experiment with a very large number of replicates
d <- sim.mm(replicates = 500, treatment = 0.05)
plot.sim(d)

# recover fixed and random parameters using `lmer`
# y1
m1 <- lmer(y1 ~ treatment + (treatment|animalid), data = d)
ranef(m1)
fixef(m1)  # so far, so good


anova(m1, ddf = "Kenward-Roger")  # fail to reject null
anova(lm(y1 ~ treatment, data = d))  # confidently reject null

# Depending on which way you look at it, we should reject the null - it's just a tiny, tiny effect
# but it's 3 animals...
# what are the degrees of freedom in the ANOVA table for the single-level regression? 
# Does that look right?

# Most of the time our data suffer from errors that aren't perfectly centred around the mean
# what happens if there is absolutely NO EFFECT at all?







