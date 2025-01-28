# GROUND CLASS DISTRIBUTIONS
# package imports
box::use(
  dplyr[...],
  stats[...],
  ggplot2[ggplot, aes, geom_point, geom_path, geom_density]
)

## =============================================================================
## deriving GSI index distribution from tunnel data

# lengths of tunnels in different GSI categories (from Paraskevoupoulos et al, 2012)
cat_A_m <- 676
cat_B_m <- 570 + 61 + 310 +1606 + 1715
cat_C_m <- 600 + 741 + 237 + 8 + 475 + 401 + 65 + 515
cat_D_m <- 34 + 153 + 150 + 72 + 140 + 153

# definition of the categories
# we are taking into consideration GSI indices from 0 to 70 
# this will lead to an underestimation of good ground conditions
cat_A_GSI <- 55:70
cat_B_GSI <- 35:55
cat_C_GSI <- 15:35
cat_D_GSI <- 0:15

x = seq(1, 70, 0.1)

# creating a "fake" step-wise distribution of GSIs, on which to fit a continuous distribution
cat_A <- runif(cat_A_m, min = min(cat_A_GSI), max = max(cat_A_GSI))
cat_B <- runif(cat_B_m, min = min(cat_B_GSI), max = max(cat_B_GSI))
cat_C <- runif(cat_C_m, min = min(cat_C_GSI), max = max(cat_C_GSI))
cat_D <- runif(cat_D_m, min = min(cat_D_GSI), max = max(cat_D_GSI))

all_cat <- c(cat_A, cat_B, cat_C, cat_D)

# fit various distributions to see which one fits the data
dist_data_weibull <- MASS::fitdistr(all_cat, "weibull")
dist_prob_weibull <- dweibull(x, shape = dist_data_weibull$estimate[["shape"]], scale = dist_data_weibull$estimate[["scale"]])

dist_data_lnorm <- MASS::fitdistr(all_cat, "lognormal")
dist_prob_lnorm <- dlnorm(x, meanlog = dist_data_lnorm$estimate[["meanlog"]], sdlog = dist_data_lnorm$estimate[["sdlog"]])

dist_data_gamma <- MASS::fitdistr(all_cat, "gamma")
dist_prob_gamma <- dgamma(x, shape = dist_data_gamma$estimate[["shape"]], rate = dist_data_gamma$estimate[["rate"]])

dist_data_norm <- MASS::fitdistr(all_cat, "normal")
dist_prob_norm <- dnorm(x, mean = dist_data_norm$estimate[["mean"]], sd = dist_data_norm$estimate[["sd"]])

# we found (empirically) that the best way to fit the distribution to the data
# was to fit it to the reverse of the data and then invert it
all_cat_data <- tibble::tibble(data = all_cat, dist = "orig") |> 
    mutate(data = abs(75 - data))

weibull_dist_data <- tibble::tibble(data = dist_prob_weibull, dist = "weibull")
lognorm_dist_data <- tibble::tibble(data = dist_prob_lnorm, dist = "lognormal")
gamma_dist_data <- tibble::tibble(data = dist_prob_gamma, dist = "gamma")
norm_dist_data <- tibble::tibble(data = dist_prob_norm, dist = "norm")

cat_db <- all_cat_data
database <- rbind(weibull_dist_data, lognorm_dist_data, gamma_dist_data, norm_dist_data) |> 
    mutate(x = rep(x, 4))

# as can be seen from the inverted graphs, the Weibull distribution is the one 
# that best approximates the data
distr_plots <- ggplot() + 
    geom_density(data = cat_db, aes(x = data, color = dist)) +
    geom_path(data = database, aes(x, data, color = dist))

# now we can plot the actual distributions
all_cat_data_rev <- all_cat_data |> mutate(data = abs(75-data))
dist_prob_weibull_rev <- dweibull(x, shape = dist_data_weibull$estimate[["shape"]], scale = dist_data_weibull$estimate[["scale"]]) |> rev() 
weibull_dist_data_rev <- tibble::tibble(data = dist_prob_weibull_rev, dist = "weibull")

real_distr_plots <- ggplot() + 
    geom_density(data = all_cat_data_rev, aes(data, color = dist)) +
    geom_path(data = weibull_dist_data_rev, aes(x, data, color = dist))



## =============================================================================
## discretizing GSI index distribution

# definition of the discretization interval
classes <- seq(0, 70, 5)

# computation of probabilities for each class
class_prob <- pweibull(classes, shape = dist_data_weibull$estimate[["shape"]], scale = dist_data_weibull$estimate[["scale"]]) |> rev()
class_prob <- class_prob[-length(class_prob)] - class_prob[-1]

# "normalization" of the results, to have the sum of all probs equal to 1
class_prob <- class_prob / sum(class_prob)



## =============================================================================
## definition of conditional probability tables for the ground classes 
## transition between unit length intervals

# we'll define the cpts based on the normal distribution, assuming that it is 
# more likely for adjacent unit length intervals to have similar ground classes
p <- seq(-3.75, 0, 0.75)
cond_prob <- pnorm(p[-1]) - pnorm(p[-length(p)]) |> abs()

# double the middle probability and mirror the distribution
cond_prob[length(cond_prob)] <- cond_prob[length(cond_prob)]*2
cond_prob <- c(cond_prob, rev(cond_prob[-length(cond_prob)]))

# normalize the probabilities so that their sum is equal to 1
cond_prob <- cond_prob / sum(cond_prob)

# build conditional probability table (matrix)
cpt_ground_class <- lapply(1:14, function(n) {

    if (n > 5) leading_zeros <- n - 5
    else leading_zeros <- 0
    
    if (n < 5) prob <- cond_prob[(6-n):length(cond_prob)]
    else if (n > 9) prob <- cond_prob[1:(length(cond_prob) - (n-10))]
    else prob <- cond_prob
    
    # normalize the conditional probabilities
    prob <- prob / sum(prob)

    trailing_zeros <- 14 - leading_zeros - length(prob)
    if (trailing_zeros < 1) trailing_zeros <- 0

    c(rep(0, leading_zeros), prob, rep(0, trailing_zeros))

}) |> 
    do.call(cbind, args = _)

# # we transform the initial probabilities and the cpt into matrices to perform matrix multiplication
# class_prob_m <- matrix(class_prob, nrow = 1)
# 
# # the real conditional probability table is cpt_ground_class
# # to perform the matrix multiplication and obtain the final probabilities, we have to transpose it
# final_prob <- class_prob_m %*% (cpt_ground_class |> t())
# x <- seq(2.5, 67.5, 5)
# 
# plot(x, class_prob)
# for (i in 1:1000) {
#     if (i == 1) final_probs <- class_prob_m
#     final_probs <- final_probs %*% (cpt_ground_class |> t())
#     if (!(i %% 10)) lines(x, final_probs)
# }
# 
# print(final_probs)
# print(class_prob)



## =============================================================================
## random sampling
## define a ground class random sampling function based on our probability 
## tables and distributions
##
## define a function to normalize and re-sample the results on our interval
## 
## define a function to smoothen the sample (via moving averages)

#' @export
ground_class_sample <- function(n, prob = class_prob, cpt = cpt_ground_class) {
  
  # take the first ground class from our distribution
  prev_val <- sample.int(length(prob), 1, prob = class_prob)
  
  # pre-assign the final vector to avoid expensive memory management
  final_vals <- numeric(n)
  final_vals[1] <- prev_val
  
  # fill the vector with values computed from our CPTs 
  for (i in 2:n) {
    
    # the sampling probability is given by the CPT and the original distribution
    next_val <- sample.int(
      length(prob), 
      1, 
      prob = cpt[, prev_val]*prob
    )
    
    final_vals[i] <- prev_val <- next_val
    
  }
  
  final_vals
  
}



#' @export
resample <- function(sample, smin = 1, smax = 14, rmax) {
  
  (sample - smin)/(smax - smin)*rmax
  
}



#' @export
smooth_sample <- function(sample, mav_filter = 3) {
  
  sample %>% stats::filter(rep(1/mav_filter, mav_filter)) %>% coalesce(., sample)
  
}
