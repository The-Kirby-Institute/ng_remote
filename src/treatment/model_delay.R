library(tidyverse)


# data
data = 
  tibble(delay = c(0, 2, 7, 4*30),
         prob = c(0.566, 0.783, 0.874, 0.968)/0.968)   # intervention
         # prob = c(0.274, 0.301, 0.472, 0.857)/0.857) # baseline


# Fit a model
logistic = function(x, t = data$delay){ 1/(1 + exp(-x[1]*(t - x[2])))}
rss = function(x){ 
  fit = logistic(x)
  sum((fit - data$prob)^2)
}
x0 = c(0.001, 21)
x = optim(x0, rss)
fit = 
  tibble(delay = seq(0, 125)) %>%
  mutate(prob = logistic(x$par, delay))


# Check
data %>%
  ggplot(aes(x = delay, y = prob)) +
  geom_point() +
  geom_line(data = fit) +
  labs(x = 'Delay from test to Treatment',
       y = 'Probability',
       title = 'Probability of Delay to Treatment Conditioned on Treatment')

