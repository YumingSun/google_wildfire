library(MASS)
library(tidyverse)
library(TestIndVars)
library(spdep)
library(pracma)
library(INLA)

# -------- helpers --------
simulate_icar = function(nrow = 3, ncol = 4, sigma_u = 1){
  nb = cell2nb(nrow, ncol)
  W  = nb2mat(nb, style = "B")
  Q  = diag(rowSums(W)) - W
  Q_pinv = pinv(Q)
  u = as.vector(mvtnorm::rmvnorm(1, sigma = (sigma_u^2) * Q_pinv))
  u = u - mean(u)
  return(list(u = u, W = W))
}

simulate_rw2 = function(n, sigma_v){
  v = numeric(n)
  v[1:2] = rnorm(2, 0 , sigma_v)
  for (t in 3:n) v[t] = 2*v[t-1] - v[t-2] + rnorm(1, 0, sigma_v)
  v = v - mean(v)
  v = v / sd(v) * 0.5
  return(v)
}

# -------- parameters --------
parameters = as.numeric(commandArgs(trailingOnly = TRUE))
expid      = parameters[1]
n_spatial  = parameters[2] # 7, 14, 21
n_time     = parameters[3] # 180, 260, 320
gamma_shape = if (length(parameters) >= 4) parameters[4] else 5

n        = n_spatial * n_time
sigma_u  = 1
sigma_v  = 0.05
sigma_w  = 0.1
p        = 10

ncol = n_spatial / 7
icar_lst = simulate_icar(nrow = 7, ncol = ncol, sigma_u = sigma_u)
u = icar_lst$u
g = icar_lst$W
v = simulate_rw2(n = n_time, sigma_v = sigma_v)

space_id = rep(1:n_spatial, times = n_time)
time_id  = rep(1:n_time, each  = n_spatial)
st_id    = as.numeric(interaction(space_id, time_id, drop = TRUE))

w = rnorm(n, mean = 0, sd = sigma_w)
w = w - mean(w)

X = mvrnorm(n, mu = rep(0,p), Sigma = diag(rep(1,p)))

dat <- data.frame(
  X,
  space = space_id,
  time  = time_id,
  st_id = st_id
)

dat_lag = dat %>%
  group_by(space) %>%
  mutate(
    X1_lag1  = lag(X1,  n = 1, order_by = time, default = 0),
    X1_lag2  = lag(X1,  n = 2, order_by = time, default = 0),
    X1_lag3  = lag(X1,  n = 3, order_by = time, default = 0),
    
    X2_lag1  = lag(X2,  n = 1, order_by = time, default = 0),
    X2_lag2  = lag(X2,  n = 2, order_by = time, default = 0),
    X2_lag3  = lag(X2,  n = 3, order_by = time, default = 0),
    
    X3_lag1  = lag(X3,  n = 1, order_by = time, default = 0),
    X3_lag2  = lag(X3,  n = 2, order_by = time, default = 0),
    X3_lag3  = lag(X3,  n = 3, order_by = time, default = 0),
    
    X4_lag1  = lag(X4,  n = 1, order_by = time, default = 0),
    X4_lag2  = lag(X4,  n = 2, order_by = time, default = 0),
    X4_lag3  = lag(X4,  n = 3, order_by = time, default = 0),
    
    X5_lag1  = lag(X5,  n = 1, order_by = time, default = 0),
    X5_lag2  = lag(X5,  n = 2, order_by = time, default = 0),
    X5_lag3  = lag(X5,  n = 3, order_by = time, default = 0),
    
    X6_lag1  = lag(X6,  n = 1, order_by = time, default = 0),
    X6_lag2  = lag(X6,  n = 2, order_by = time, default = 0),
    X6_lag3  = lag(X6,  n = 3, order_by = time, default = 0),
    
    X7_lag1  = lag(X7,  n = 1, order_by = time, default = 0),
    X7_lag2  = lag(X7,  n = 2, order_by = time, default = 0),
    X7_lag3  = lag(X7,  n = 3, order_by = time, default = 0),
    
    X8_lag1  = lag(X8,  n = 1, order_by = time, default = 0),
    X8_lag2  = lag(X8,  n = 2, order_by = time, default = 0),
    X8_lag3  = lag(X8,  n = 3, order_by = time, default = 0),
    
    X9_lag1  = lag(X9,  n = 1, order_by = time, default = 0),
    X9_lag2  = lag(X9,  n = 2, order_by = time, default = 0),
    X9_lag3  = lag(X9,  n = 3, order_by = time, default = 0),
    
    X10_lag1 = lag(X10, n = 1, order_by = time, default = 0),
    X10_lag2 = lag(X10, n = 2, order_by = time, default = 0),
    X10_lag3 = lag(X10, n = 3, order_by = time, default = 0)
  ) %>% ungroup()

eta =
  1 +  0.8 * dat_lag$X1 - 0.5 * dat_lag$X1 * dat_lag$X2 + 0.75 * dat_lag$X1 * dat_lag$X3 +
  0.13 * dat_lag$X1 * dat_lag$X4 - 0.25 * dat_lag$X1 * dat_lag$X5 - 0.66 * dat_lag$X1 * dat_lag$X6 +
  0.4  * dat_lag$X1 * dat_lag$X7 - 0.7  * dat_lag$X1 * dat_lag$X8 + 0.9  * dat_lag$X1 * dat_lag$X9 -
  0.2  * dat_lag$X1 * dat_lag$X10 +
  
  0.8 * dat_lag$X1_lag1 - 0.5 * dat_lag$X1_lag1 * dat_lag$X2_lag1 + 0.75 * dat_lag$X1_lag1 * dat_lag$X3_lag1 +
  0.13 * dat_lag$X1_lag1 * dat_lag$X4_lag1 - 0.25 * dat_lag$X1_lag1 * dat_lag$X5_lag1 - 0.66 * dat_lag$X1_lag1 * dat_lag$X6_lag1 +
  0.4  * dat_lag$X1_lag1 * dat_lag$X7_lag1 - 0.7  * dat_lag$X1_lag1 * dat_lag$X8_lag1 + 0.9  * dat_lag$X1_lag1 * dat_lag$X9_lag1 -
  0.2  * dat_lag$X1_lag1 * dat_lag$X10_lag1 +
  
  0.8 * dat_lag$X1_lag2 - 0.5 * dat_lag$X1_lag2 * dat_lag$X2_lag2 + 0.75 * dat_lag$X1_lag2 * dat_lag$X3_lag2 +
  0.13 * dat_lag$X1_lag2 * dat_lag$X4_lag2 - 0.25 * dat_lag$X1_lag2 * dat_lag$X5_lag2 - 0.66 * dat_lag$X1_lag2 * dat_lag$X6_lag2 +
  0.4  * dat_lag$X1_lag2 * dat_lag$X7_lag2 - 0.7  * dat_lag$X1_lag2 * dat_lag$X8_lag2 + 0.9  * dat_lag$X1_lag2 * dat_lag$X9_lag2 -
  0.2  * dat_lag$X1_lag2 * dat_lag$X10_lag2 +
  
  0.8 * dat_lag$X1_lag3 - 0.5 * dat_lag$X1_lag3 * dat_lag$X2_lag3 + 0.75 * dat_lag$X1_lag3 * dat_lag$X3_lag3 +
  0.13 * dat_lag$X1_lag3 * dat_lag$X4_lag3 - 0.25 * dat_lag$X1_lag3 * dat_lag$X5_lag3 - 0.66 * dat_lag$X1_lag3 * dat_lag$X6_lag3 +
  0.4  * dat_lag$X1_lag3 * dat_lag$X7_lag3 - 0.7  * dat_lag$X1_lag3 * dat_lag$X8_lag3 + 0.9  * dat_lag$X1_lag3 * dat_lag$X9_lag3 -
  0.2  * dat_lag$X1_lag3 * dat_lag$X10_lag3 +
  
  u[space_id] + v[time_id] + w[st_id]

mu = exp(eta) 
y  = rgamma(n, shape = gamma_shape, rate = gamma_shape / mu)

dat_lag$y = y

formula = y ~ X1 + X1:X2 + X1:X3 + X1:X4 + X1:X5 + X1:X6 + X1:X7 + X1:X8 +
  X1:X9 + X1:X10 +
  
  X1_lag1 + X1_lag1:X2_lag1 + X1_lag1:X3_lag1 + X1_lag1:X4_lag1 + X1_lag1:X5_lag1 +
  X1_lag1:X6_lag1 + X1_lag1:X7_lag1 + X1_lag1:X8_lag1 + X1_lag1:X9_lag1 + X1_lag1:X10_lag1 +
  
  X1_lag2 + X1_lag2:X2_lag2 + X1_lag2:X3_lag2 + X1_lag2:X4_lag2 + X1_lag2:X5_lag2 +
  X1_lag2:X6_lag2 + X1_lag2:X7_lag2 + X1_lag2:X8_lag2 + X1_lag2:X9_lag2 + X1_lag2:X10_lag2 +
  
  X1_lag3 + X1_lag3:X2_lag3 + X1_lag3:X3_lag3 + X1_lag3:X4_lag3 + X1_lag3:X5_lag3 +
  X1_lag3:X6_lag3 + X1_lag3:X7_lag3 + X1_lag3:X8_lag3 + X1_lag3:X9_lag3 + X1_lag3:X10_lag3 +
  
  f(space, model = "besag", graph = g, constr = TRUE) +
  f(time,  model = "rw2",   constr = TRUE) +
  f(st_id, model = "iid",   constr = TRUE)

start = proc.time()  
result = inla(
  formula,
  family = "gamma",
  data = dat_lag,
  control.family   = list(link = "log"),
  control.predictor= list(compute = TRUE),
  control.compute  = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE),
  control.fixed    = list(
    mean = list(default = 0, intercept = 0),
    prec = list(default = 1e-4, intercept = 1e-4)
  )
)
end = proc.time() 
time_spent = (end - start)[["elapsed"]]

fixed_effects   = result$summary.fixed[, 1:5]
random_spatial  = result$summary.random$space[, 2:6]
random_temporal = result$summary.random$time[,  2:6]

res_loc = paste0("/results/",
                 sprintf("sp_%02d_tp_%03d/", n_spatial, n_time))
fixed_res_name = sprintf("fixed_coef_%03d.txt", expid)
sp_res_name    = sprintf("random_sp_%03d.txt", expid)
tp_res_name    = sprintf("random_tp_%03d.txt", expid)
run_time_name  = sprintf("time_%03d.txt", expid)

write.table(fixed_effects,    paste0(res_loc, fixed_res_name), sep = "\t", row.names = TRUE)
write.table(random_spatial,   paste0(res_loc, sp_res_name),    sep = "\t", row.names = TRUE)
write.table(random_temporal,  paste0(res_loc, tp_res_name),    sep = "\t", row.names = TRUE)
write(time_spent, file = paste0(res_loc, run_time_name))
