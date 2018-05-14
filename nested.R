library(magrittr)
library(dplyr)
library(ggplot2)
library(lme4)
library(gridExtra)
library(lmerTest)
library(viridis)
ggplot <- function(...) ggplot2::ggplot(...) + scale_colour_viridis(discrete = T) + scale_fill_viridis(discrete = T)

set.seed(1)


d <- expand.grid(animalid  = c(-0.5, 0, 0.5),
                 treatment = c(0, 0.5),
                 replicate = c(1:10))
# large effect sizes (random and fixed)
d$y1 <- d$treatment + d$animalid + rnorm(nrow(d), 0, 0.2)
# small effect sizes (random and fixed)
d$y2 <- d$treatment + d$animalid + rnorm(nrow(d), 0, 2)
d$animalid %<>% factor(labels = c("01", "02", "03"))
d$treatment %<>% factor(labels = c("a", "b"))
d %<>% group_by(treatment, animalid) %>% mutate(y1mean = mean(y1), y2mean = mean(y2))

p1 <- ggplot(d, aes(x = treatment, y = y1, colour = animalid, fill = animalid)) +
  geom_point(position = position_dodge(width = 1/2), alpha = 3/4) +
  geom_point(aes(y = y1mean), size = 4, position = position_dodge(width = 1/2), shape = 3) +
  ylim(-3, 3)
p2 <- ggplot(d, aes(x = treatment, y = y2, colour = animalid)) +
  geom_point(position = position_dodge(width = 1/2), alpha = 3/4) +
  geom_point(aes(y = y2mean), size = 4, position = position_dodge(width = 1/2), shape = 3) +
  ylim(-3, 3)

grid.arrange(p1, p2, ncol = 2)



# recover fixed and random parameters using `lmer`
# y1
m1 <- lmer(y1 ~ treatment + (treatment|animalid), data = d)
ranef(m1)
fixef(m1)  # so far, so good

# y2
m2 <- lmer(y2 ~ treatment + (treatment|animalid), data = d)
ranef(m2)
fixef(m2)  # obviously not as good as y1, but not terrible

anova(m1)  # justifiably reject null
anova(m2)  

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


d <- expand.grid(animalid  = c(-0.5, 0, 0.5),
                 treatment = c(0, 0.5),
                 replicate = c(1:500))
# large effect sizes (random and fixed)
d$y1 <- d$treatment + d$animalid + rnorm(nrow(d), 0, 0.2)
# small effect sizes (random and fixed)
d$y2 <- d$treatment + d$animalid + rnorm(nrow(d), 0, 2)
d$animalid %<>% factor(labels = c("01", "02", "03"))
d$treatment %<>% factor(labels = c("a", "b"))
d %<>% group_by(treatment, animalid) %>% mutate(y1mean = mean(y1), y2mean = mean(y2))

p1 <- ggplot(d, aes(x = treatment, y = y1, colour = animalid, fill = animalid)) +
  geom_point(position = position_dodge(width = 1/2), alpha = 3/4) +
  geom_point(aes(y = y1mean), size = 4, position = position_dodge(width = 1/2), shape = 3) +
  ylim(-3, 3)
p2 <- ggplot(d, aes(x = treatment, y = y2, colour = animalid)) +
  geom_point(position = position_dodge(width = 1/2), alpha = 3/4) +
  geom_point(aes(y = y2mean), size = 4, position = position_dodge(width = 1/2), shape = 3) +
  ylim(-3, 3)

grid.arrange(p1, p2, ncol = 2)



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
