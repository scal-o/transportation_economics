## COST FUNCTIONS
# package imports
box::use(
  dplyr[...],
  stats[rnorm, rbeta],
  tidyr[pivot_longer],
  utils[read.csv2]
)


## =============================================================================
## cost data from Swiss tunnels

min_cost_2001 <- 110 # €/m^3 2001€
max_cost_2001 <- 1100 # €/m^3 2001€


# switzerland eurostat construction cost index
cci_eurostat <- read.csv2("Eurostat_CCI.csv") %>% 
  as_tibble() %>% 
  mutate(across(starts_with("x"), \(x) if (is.character(x)) sub(",", ".", x) else x)) %>% 
  mutate(across(starts_with("x"), as.numeric)) %>%
  filter(State == "Switzerland") %>% 
  select(-State, -X2000) %>% 
  pivot_longer(cols = everything(), values_to = "CCI", names_to = "year", names_prefix = "X") %>% 
  mutate(
    CCI_percent = 1 + CCI/100,
    CCI_cumulative = cumprod(CCI_percent)
    )

cci_increase <- cci_eurostat$CCI_cumulative[nrow(cci_eurostat)]

min_cost <- min_cost_2001*cci_increase
max_cost <- max_cost_2001*cci_increase



## =============================================================================
## functions to compute various cost data

# function to compute cost/m^3, based on data from Greek tunnels
greek_cost_fun <- function(x) {
  
  x[x < 1] <- 1
  -55.99 * log(x) + 259.07
  
}


x <- 1:70
greek_costs <- greek_cost_fun(x)

# function to compute normalized cost/m^3
normalized_cost_fun <- function(x) {
    (greek_cost_fun(x) - min(greek_costs))/(max(greek_costs) - min(greek_costs))
}


# function to compute small disturbances
noise_fun <- function(x) {
  
  noise_l <- sample(c(0, 1), length(x), replace = TRUE, prob = c(0.9, 0.1))
  noise <- rbeta(length(x), shape1 = 3, shape2 = 3) * 4 + 1
  
  noise[!noise_l] <- 1
  noise
}


#' @export
#' function to compute actual cost for tunnels in Switzerland, based on cost
#' data from 2001
cost_fun <- function(x) {
    normalized_cost_fun(x) * noise_fun(x) * max_cost
}
