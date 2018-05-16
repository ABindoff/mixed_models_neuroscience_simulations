library(magrittr)
library(dplyr)
library(ggplot2)
library(lme4)
library(gridExtra)
library(lmerTest)
library(viridis)
library(ggbeeswarm)
ggplot <- function(...) ggplot2::ggplot(...) + scale_colour_viridis(discrete = T, begin = 0.1, end = 0.5) + scale_fill_viridis(discrete = T, end = 0.5)
options(digits = 5)
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
    geom_point(aes(y = y1mean, group = animalid), colour = "black", size = 4, position = position_dodge(width = 1), shape = 3) +
    ylim(-3, 3)
  p2 <- ggplot(d, aes(x = treatment, y = y2, colour = animalid)) +
    geom_quasirandom(dodge.width = 1, alpha = 1/3) +
    geom_point(aes(y = y2mean, group = animalid), colour = "black", size = 4, position = position_dodge(width = 1), shape = 3) +
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
anova(lm(y1 ~ treatment, data = d))  # confidently reject null?

# Depending on which way you look at it, we should reject the null - it's just a tiny, tiny effect
# but it's 3 animals...
# what are the degrees of freedom in the ANOVA table for the single-level regression? 
# Does that look right?
# (put aside the fact that we have perfectly symmetric noise)

# Most of the time our data suffer from errors that aren't perfectly centred around the mean difference
# what happens if there is absolutely NO EFFECT at all, and there really is a very good chance
# that differences in mean are caused entirely by noise?


sim.mm2 <- function(replicates = 6, treatment = 0.5, sd1 = 0.5, sd2 = 1){
  replicates <- 2*(floor(replicates/2))  # make sure it's an even number of replicates
  
  d <- expand.grid(animalid  = c(-0.5, 0.1, 0.4),
                   treatment = c("a", "b"),
                   replicate = c(1:replicates))
  d %<>% mutate(effect = 0)
  d$effect[d$treatment == "b"] <- treatment
  
  d$y1 <- d$effect + d$animalid + rnorm(nrow(d), 0, sd1)
  d$y2 <- d$effect + d$animalid + rnorm(nrow(d), 0, sd2)
 
  d$animalid %<>% factor(labels = c("01", "02", "03"))
  d %<>% group_by(treatment, animalid) %>% mutate(y1mean = mean(y1), y2mean = mean(y2))
  d
}



d <- sim.mm2(replicates = 10, treatment = 0.0)
plot.sim(d)

p1 <- NULL
p2 <- NULL

for(i in 1:100){
  d <- sim.mm2(replicates = 100, treatment = 0.0, sd1 = 2)
  a1 <- anova(lmer(y1 ~ treatment + (treatment|animalid), data = d), ddf = "Kenward-Roger")
  a2 <- anova(lm(y1 ~ treatment, data = d))
  p1 <- c(p1, a1$`Pr(>F)`)
  p2 <- c(p2, a2$`Pr(>F)`[1])
}
hist(p1, breaks = 20)
hist(p2, breaks = 20)
mean(p1 < .05)
mean(p2 < .05)

# we can see that both methods produce acceptable Type 1 error rates
# (actually might be more concerned about the conservatism of the Kenward-Roger
# approximation in this example)

# what happens when we reduce the number of replicates?

p1 <- NULL
p2 <- NULL
p3 <- NULL
p4 <- NULL

for(i in 1:1000){
  d <- sim.mm2(replicates = 5, treatment = 0.0, sd1 = 2)
  a1 <- anova(lmer(y1 ~ treatment + (treatment|animalid), data = d), ddf = "Kenward-Roger")
  a2 <- anova(lm(y1 ~ treatment, data = d))
  a3 <- anova(lmer(y1 ~ treatment + (1|animalid), data = d), ddf = "Kenward-Roger")
  a4 <- anova(lmer(y1 ~ treatment + (1|animalid), data = d), ddf = "Satterthwaite")
  p1 <- c(p1, a1$`Pr(>F)`)
  p2 <- c(p2, a2$`Pr(>F)`[1])
  p3 <- c(p3, a3$`Pr(>F)`)
  p4 <- c(p4, a4$`Pr(>F)`)
}
hist(p1, breaks = 20)
hist(p2, breaks = 20)
hist(p3, breaks = 20)
hist(p4, breaks = 20)

mean(p1 < .05)
mean(p2 < .05)
mean(p3 < .05)
mean(p4 < .05)

# these observations are unreservedly dependent, however the Type 1 error rate remains exquisitely controlled
# no matter how you analyse the data



