# METADATA ====
# Description: Line search methods homework (Deterministic models and optimization)
# Created: 2018-11-15 (Reid Falconer)
# Updated: 2018-11-20 (Reid Falconer)
# Reviewed: 

# INITIALISE ====

rm(list=ls())

#> Libraries ----
library(tidyverse) # always load tidyverse last to avoid namespace conflicts

#> Set options ----

# disable scientific notation
options(scipen = 999)


# FUNCTIONS ====

#> gradient descent ----

# definition of the gradient descent
gradient_descent <- function(p, start, n = 2, step_1=1, 
                             tol=1e-8, max_it = 100000, 
                             display = T ) {
  if(n==2) { # our 2 two-dimensional case
    m <- c(0.5, 0) # define m vector
    Amat <- function(p) {
    Amat <- matrix(c(1,p,p,1),2,2) 
    }
    A <- Amat(p) # define A matrix
    }
    else{  # mutiple dimensions
      m <- runif(n, -1, 1)# define m vector
      Amat <- matrix(rnorm(n*n, mean = 0, sd = 0.1), ncol=n, nrow = n)
      A <- Amat%*%(t(Amat)) # define A matrix (making sure that it is a fixed positive definite matrix)
  }
  # define function
  defined_function <- function(x, m, A) {
    func = 0.5*t(x - m)%*%A%*%(x - m) - sum(log((x)^2)) 
    return(c(func))
  }
  # define the derivative 
  derivative <- function(x, m, A) {
    grad <- A%*%(x-m) - 2/(x)
    return(c(grad))
  }
  i <- 0 # Counter for number of iterations
  point_1 <- start
  hi <- Sys.time() # before convergence 
  grdnt <- derivative(point_1, m, A)
  step = step_1
  point_2 <- point_1 - step_1*grdnt
  path_mat <- matrix(NA, nrow = max_it, ncol = n)
  path_mat[i+1,] <- start
  while (abs(defined_function(point_1, m, A)-defined_function(point_2, m, A)) > tol) { 
    while(defined_function(point_2, m, A) > defined_function(point_1, m, A)) {
      step = step*0.5
      point_2 <- point_1 - step*grdnt
    } 
    i <- i + 1
    path_mat[i+1,] <- point_2
    step = step_1
    point_1 <- point_2
    grdnt <- derivative(point_1, m, A)
    point_2 <- point_1 - step*grdnt
    #print(defined_function(point_2, m, A)) # optional of one wants to print progress
  }
  paths = path_mat[0:i+1,] %>% as.data.frame()
  point = point_2
  iterations = i
  point_text = paste0(round(point[1],2)," ; ", round(point[2],2))
  
  if(display == T) {
    x <- seq(-5, 5, length.out = 100)
    y <- seq(-5, 5, length.out = 100)
    z <- apply(expand.grid(x, y), 1, "defined_function", m = m, A = A)
    
    df <- cbind(expand.grid(x, y),z) %>% 
      dplyr::rename(x = Var1,
                    y = Var2) %>% 
      as.data.frame() 
    
    plot <- ggplot() + geom_raster(data = df, aes(x = x, y = y, fill = z), 
                                   interpolate = TRUE, show.legend = FALSE) +
      geom_contour(data = df, colour = "white", aes(x = x, y = y, z = z), 
                   binwidth = 2, alpha = 0.4) + 
      geom_path(aes(x = paths$V1, y = paths$V2, colour = "Gradient Descent"), 
                size = 2,  linetype = 1) +
      geom_point(size = 3, color = "white", aes(x = point[1], point[2])) +
      geom_point(size = 2, color = "black", aes(x = paths$V1, y = paths$V2),shape = 1) +
      theme_grey() +
      labs(x = "", y = "") +
      geom_label(color = "black", aes(x = point[1], point[2]), 
                label = point_text, nudge_x = 0, nudge_y = 1, size = 4)+
      scale_colour_manual("", values = c("Gradient Descent"="#5addd0")) +
      theme(axis.text=element_text(size=12), plot.title = element_text(size=22))
  }
  else {
    na <- Sys.time() # after convergence 
    diff <- na-hi # calculate convergence rate
    return(list(point = point_2, iterations = i, point_text = point_text, 
                paths = path_mat[0:i+1,] %>%
                  as.data.frame() %>%
                    rename(x = V1,
                          y = V2),
                    time = diff))
  }
  print(plot) # return the plot
  na <- Sys.time() # after convergence
  diff <- na-hi # calculate convergence rate
  return(list(point = point_2, iterations = i, point_text = point_text, 
              paths = path_mat[0:i+1,] %>%
                as.data.frame() %>%
                rename(x = V1,
                      y = V2),
              time = diff)) 
  }

# locate the minimum of the function using the Gradient Descent method
gd_results <- gradient_descent(
  p = 0.5, # define rho (only applicable in 2 dimensional case)
  start = runif(2, -5, 5), # start point of the search. Initialised point in paper c(-3.5,1)
  step_1 = 1, # step size (alpha)
  tol = 1e-5, # tolerance
  display = T, # display plot (only applicable in 2 dimensional case)
  n = 2, # dimensions
  max_it = 100000) # max iterations 


#> newton's method ----

# definition of newton's method
newtons_method <- function(p, start, n = 2, step_1=1, 
                           tol=1e-8, max_it = 100000, 
                           display = T ) {
  if(n==2) {
    m <- c(0.5, 0) # define m vector
    Amat <- function(p) {
      Amat <- matrix(c(1,p,p,1),2,2) # define A matrix
    }
    A <- Amat(p)
  }
  else{
    m <- runif(n, -1, 1) # define m vector
    Amat <- matrix(rnorm(n*n, mean = 0, sd = 0.1), ncol=n, nrow = n) # define A matrix
    A <- Amat%*%(t(Amat)) # making sure that it is a fixed positive definite matrix
  }
  # define function
  defined_function <- function(x, m, A) {
    func = 0.5*t(x - m)%*%A%*%(x - m) - sum(log((x)^2)) 
    return(c(func))
  }
  # define the derivative 
  derivative <- function(x, m, A) {
    grad <- A%*%(x-m) - 2/(x)
    return(c(grad))
  }
  #define the hessian 
  hessian <- function(x, m, A) {
     A + diag((diag(2/x%*%t(x))))
  }
  i <- 0 # Counter for number of iterations
  point_1 <- start
  hi <- Sys.time() # before convergence
  grdnt <- derivative(point_1, m, A)
  hess <- hessian(point_1, m, A)
  step = step_1
  point_2 <- point_1 - step_1*solve(hess)%*%grdnt
  path_mat <- matrix(NA, nrow = max_it, ncol = n)
  path_mat[i+1,] <- start
  while (abs(defined_function(point_1, m, A)-defined_function(point_2, m, A)) > tol) {
    while(defined_function(point_2, m, A) > defined_function(point_1, m, A)) {
      step <- step*0.5
      point_2 <- point_1 - step*solve(hess)%*%grdnt
    } 
    i <- i + 1
    path_mat[i+1,] <- point_2
    step = step_1
    point_1 <- point_2
    grdnt <- derivative(point_1, m, A)
    hess <- hessian(point_1, m, A)
    point_2 <- point_1 - step_1*solve(hess)%*%grdnt
    #print(defined_function(point_2, m, A)) # print progress
  }
  paths = path_mat[0:i+1,] %>% as.data.frame()
  point = point_2
  iterations = i
  point_text = paste0(round(point[1],2)," ; ", round(point[2],2))
  
  if(display == T) {
    x <- seq(-5, 5, length.out = 100)
    y <- seq(-5, 5, length.out = 100)
    z <- apply(expand.grid(x, y), 1, "defined_function", m = m, A = A)
    
    df <- cbind(expand.grid(x, y),z) %>% 
      dplyr::rename(x = Var1,
                    y = Var2) %>% 
      as.data.frame() 
    
    plot <- ggplot() + geom_raster(data = df, aes(x = x, y = y, fill = z),
                                   interpolate = TRUE, show.legend = FALSE) +
      geom_contour(data = df, colour = "white", aes(x = x, y = y, z = z), 
                   binwidth = 2, alpha = 0.4) + 
      geom_path(aes(x = paths$V1, y = paths$V2, colour = "Newton's Method"), 
                size = 2,  linetype = 1) +
      geom_point(size = 3, color = "white", aes(x = point[1], point[2])) +
      theme_grey() +
      geom_point(size = 2, color = "black", aes(x = paths$V1, y = paths$V2),shape = 1) +
      labs(x = "", y = "") +
      geom_label(color = "black", aes(x = point[1], point[2]), 
                 label = point_text, nudge_x = 0, nudge_y = 1, size = 4)+
      scale_colour_manual("", values = c("Newton's Method"="#ff6666")) +
      theme(axis.text=element_text(size=12), plot.title = element_text(size=22))
  }
  else {
    na <- Sys.time() # after convergence
    diff <- na-hi
    return(list(point = point_2, iterations = i, point_text = point_text, 
                paths = path_mat[0:i+1,] %>%
                  as.data.frame() %>%
                  rename(x = V1,
                         y = V2),
                time = diff))
  }
  print(plot) # return the plot
  na <- Sys.time()
  diff <- na-hi
  return(list(point = point_2, iterations = i, point_text = point_text, 
              paths = path_mat[0:i+1,] %>%
                as.data.frame() %>%
                rename(x = V1,
                       y = V2),
              time = diff)) 
}

# locate the minimum of the function using the Gradient Descent method
newtons_result <- newtons_method(
  p = -0.7, # define rho (only applicable in 2 dimensional case)
  start = runif(2, -5, 5), # start point of the search. Initialised point in paper c(-3.5,1)
  step_1 = 1, # initialize step size (alpha)
  tol = 1e-5, # tolerance
  display = T, # display plot (only applicable in 2 dimensional case)
  n = 2, # dimensions
  max_it = 100000) # max iterations 


# Fucntion with both newton and gradient.
# This fuction was purely built to  facilitate the plotting of 
# both algoithms on the same plot. 

combination<- function(p, start, n = 2, step_1_k = 1, step_1 = 1, 
                       tol=1e-8, max_it = 100000, display = T, 
                       Gradient = T, Newton = T) {
  if(n==2) {
    m <- c(0.5, 0) 
    Amat <- function(p) {
      Amat <- matrix(c(1,p,p,1),2,2)
    }
    A <- Amat(p) # define A matrix
  }
  else{
    m <- runif(n, -1, 1)
    Amat <- matrix(rnorm(n*n, mean = 0, sd = 0.1), ncol=n, nrow = n)
    A <- Amat%*%(t(Amat)) # define A matrix (making sure that it is a fixed positive definite matrix)
  }
  # define function
  defined_function <- function(x, m, A) {
    func = 0.5*t(x - m)%*%A%*%(x - m) - sum(log((x)^2)) 
    return(c(func))
  }
  # define the derivative 
  derivative <- function(x, m, A) {
    grad <- A%*%(x-m) - 2/(x)
    return(c(grad))
  }
  #define the hessian 
  hessian <- function(x, m, A) {
    A + diag((diag(2/x%*%t(x))))
  }

  if(Newton == T) {
  i <- 0 # Counter for number of iterations
  point_1 <- start
  hi <- Sys.time() # before convergence
  grdnt <- derivative(point_1, m, A)
  hess <- hessian(point_1, m, A)
  step = step_1
  point_2 <- point_1 - step_1*solve(hess)%*%grdnt
  path_mat <- matrix(NA, nrow = max_it, ncol = n)
  path_mat[i+1,] <- start
  while (abs(defined_function(point_1, m, A)-defined_function(point_2, m, A)) > tol) {
    while(defined_function(point_2, m, A) > defined_function(point_1, m, A)) {
      step = step*0.5
      point_2 <- point_1 - step*solve(hess)%*%grdnt
    } 
    i <- i + 1
    path_mat[i+1,] <- point_2
    step = step_1
    point_1 <- point_2
    grdnt <- derivative(point_1, m, A)
    hess <- hessian(point_1, m, A)
    point_2 <- point_1 - step_1*solve(hess)%*%grdnt
    #print(defined_function(point_2, m, A)) # print progress
  }
  paths = path_mat[0:i+1,] %>% as.data.frame()
  point = point_2
  iterations = i
  point_text = paste0(round(point[1],2)," ; ", round(point[2],2))
  }
  
  if(Gradient == T) {
  k <- 0 # Counter for number of iterations
  point_1_k <- start
  hi_k <- Sys.time() # before convergence
  grdnt_k <- derivative(point_1_k, m, A)
  step_k <- step_1_k
  point_2_k <- point_1_k - step_1_k*grdnt_k
  path_mat_k <- matrix(NA, nrow = max_it, ncol = n)
  path_mat_k[k+1,] <- start
  while (abs(defined_function(point_1_k, m, A)-defined_function(point_2_k, m, A)) > tol) {
    while(defined_function(point_2_k, m, A) > defined_function(point_1_k, m, A)) {
      step_k = step_k*0.5
      point_2_k <- point_1_k - step_k*grdnt_k
    } 
    k <- k + 1
    path_mat_k[k+1,] <- point_2_k
    step_k = step_1_k
    point_1_k <- point_2_k
    grdnt_k <- derivative(point_1_k, m, A)
    point_2_k <- point_1_k - step_k*grdnt_k
    #print(defined_function(point_2_k, m, A)) # print progress
  }
  paths_k = path_mat_k[0:k+1,] %>% as.data.frame()
  point_k = point_2_k
  iterations_k = k
  point_text_k = paste0(round(point_k[1],2)," ; ", round(point_k[2],2))
  }
  
  if(display == T) {
    x <- seq(-5, 5, length.out = 100)
    y <- seq(-5, 5, length.out = 100)
    z <- apply(expand.grid(x, y), 1, "defined_function", m = m, A = A)
    df <- cbind(expand.grid(x, y),z) %>% 
      dplyr::rename(x = Var1,
                    y = Var2) %>% 
      as.data.frame() 
    title_p<- paste("p =", p)
    plot <- ggplot() + geom_raster(data = df, aes(x = x, y = y, fill = z), 
                                   interpolate = TRUE, show.legend = FALSE) +
      geom_contour(data = df, colour = "white", aes(x = x, y = y, z = z),
                   binwidth = 2, alpha = 0.4) + 
      geom_path(aes(x = paths$V1, y = paths$V2, colour = "Newton's Method") , 
                size = 2,  linetype = 1, alpha = 0.8, show.legend = T) +
      geom_path(aes(x = paths_k$V1, y = paths_k$V2, colour = "Gradient Descent"), 
                size = 2,  linetype = 1, alpha = 0.8, show.legend = T) +
      geom_point(size = 3, color = "white", aes(x = point[1], point[2])) +
      geom_point(size = 3, color = "white", aes(x = point_k[1], point_k[2])) +
      geom_point(size = 2, color = "black", aes(x = paths_k$V1, y = paths_k$V2),shape = 1) +
      geom_point(size = 2, color = "black", aes(x = paths$V1, y = paths$V2),shape = 1) +
      theme_grey() +
      scale_colour_manual("", values = c("Newton's Method"="#ff6666", 
                                         "Gradient Descent"="#5addd0")) +
      theme(legend.text=element_text(size=12), 
            plot.title = element_text(hjust = 0.5, size = 20)) +
      labs(x = "", y = "") +
      ggtitle(paste0(title_p)) 
  }
  else {
    na <- Sys.time() # after convergence
    diff <- na-hi # calculate convergence rate
    return(list(point = point_2, iterations = i, point_text = point_text, 
                paths = path_mat[0:i+1,] %>%
                  as.data.frame() %>%
                  rename(x = V1,
                         y = V2),
                time = diff))
  }
  print(plot) # return the plot
  na <- Sys.time()
  diff <- na-hi
  return(list(point = point_2, iterations = i, point_text = point_text, 
              paths = path_mat[0:i+1,] %>%
                as.data.frame() %>%
                rename(x = V1,
                       y = V2),
              time = diff)) 
}

# locate the minimum of the function using the Gradient Descent and Newtons method
newton_gd_result <- combination(
  p = -0.7, # define rho (only applicable in 2 dimensional case)
  start = runif(2, -5, 5), # start point of the search. Initialised point in paper c(-3.5,1)
  step_1_k = 1,# step size (alpha)
  step_1 = 1, # initialize step size (alpha)
  tol = 1e-5, # tolerance
  display = T, # display plot (only applicable in 2 dimensional case)
  n = 2, # dimensions
  max_it = 100000) # max iterations 

# SIMULATIONS ====

#> Study convergence plots ----

name <- 0
for (i in seq(-0.9, 0.9, 0.1)) {
  name <- name +1
  combination(
    p = i, # loops over all p (-0.9 - 0.9)
    start = c(-3.5,1), #runif(2, -5, 5), # start point of the search 
    step_1_k = 1, # step size (alpha)
    step_1 = 1,
    tol = 1e-5, # tolerance
    display = T, # display plot (only applicable in 2 dimensional case)
    n = 2, # dimensions
    max_it = 100000) # max iterations
ggsave(paste0("reports/images/fig_",name,".png"), plot = last_plot(), #save each plot to disk 
       scale = 1, dpi = 300, limitsize = TRUE)
}

#> Study higher dimensions ----

# where `number` is the number of iterations, `n` is the amount of 
# dimensions and rho is the p parameter (only applicable in 2 dimesional case)
iterations_gd <- function(number, n, rho) {
  iter <- matrix(0, nrow = number, ncol = 1)
  time <- matrix(0, nrow = number, ncol = 1)
  for(i in 1:number) {
    results <- gradient_descent(
      p = rho, # rho value (only applicable in 2 dimensional case)
      start = runif(n, -5, 5), # start point of the search 
      step_1 = 1, # step size (alpha)
      tol = 1e-5, # tolerance
      display = F, # display plot (only applicable in 2 dimensional case)
      n = n, # dimensions
      max_it = 1000000)   # max iterations
    iter[i,] <- results$iterations
    time[i,]<- results$time
  }
  return(list(iterations = mean(iter), time = mean(time), number = n))
}

# where `number` is the number of iterations, `n` is the amount of dimensions 
# and rho is the p parameter (only applicable in 2 dimesional case)
iterations_new <- function(number, n, rho) {
  iter <- matrix(0, nrow = number, ncol = 1)
  time <- matrix(0, nrow = number, ncol = 1)
  for(i in 1:number) {
    results <- newtons_method(
      p = rho, # rho value (only applicable in 2 dimensional case)
      start = runif(n, -5, 5), # start point of the search 
      tol = 1e-5,  # tolerance
      display = F,# display plot (only applicable in 2 dimensional case)
      n = n, # dimensions
      max_it = 1000000) # max iterations
    iter[i,] <- results$iterations
    time[i,]<- results$time
    #  print(iter)
  }
  return(list(iterations = mean(iter), time = mean(time), number = n))
}

# create matrix vector of correct dimensions to populate
rho_seq <- length(seq(-0.9,0.9,0.1))
col_int <- length(c(1,2,3,4,5))
output = matrix(NA, nrow = rho_seq , ncol = col_int)

index = 0
for (i in seq(-0.9,0.9,0.1)) {
  index = index + 1
interation_gd <- iterations_gd(number = 500, n = 2, rho = i) 
interation_new <- iterations_new(number = 500, n = 2, rho = i) 
print(index)
output[index,2] <- interation_gd$iterations
output[index,3] <-interation_new$iterations
output[index,1] <- i
output[index,4] <- interation_gd$time
output[index,5] <-interation_new$time
}

output %>% 
  as.data.frame() %>% 
  ggplot() +
  geom_line(aes(x = V3, y = V1, colour = "Gradient Descent"), size = 1.5) +
  geom_line(aes(x = V3, y = V2, colour = "Newton's Method"), size = 1.5) +
  theme_grey() +
  scale_colour_manual("", values = c("Newton's Method"="#ff6666", 
                                     "Gradient Descent"="#5addd0")) +
  theme(legend.text=element_text(size=12), 
        axis.text = element_text(hjust = 0.5, size = 12), 
        axis.title=element_text(size=12)) +
  labs(x = "p", y = "Iterations")

# save plot
ggsave("reports/images/iterations.png", plot = last_plot(),
       scale = 1, dpi = 300, limitsize = TRUE)

output %>% 
  as.data.frame() %>% 
  ggplot() +
  geom_line(aes(x = V1, y = V4, colour = "Gradient Descent"), size = 1.5) +
  geom_line(aes(x = V1, y = V5, colour = "Newton's Method"), size = 1.5) +
  theme_grey() +
  scale_colour_manual("", values = c("Newton's Method"="#ff6666", 
                                     "Gradient Descent"="#5addd0")) +
  theme(legend.text=element_text(size=12), 
        axis.text = element_text(hjust = 0.5, size = 12), 
        axis.title=element_text(size=12)) +
  labs(x = "p", y = "Runtime (Seconds)")
  
# save plot
ggsave("reports/images/convergence.png", plot = last_plot(),
       scale = 1, dpi = 300, limitsize = TRUE)

# create matrix vector of correct dimensions to populate
dimensions <- length(c(10,50,100,200,500,1000))
int <- length(c(1,2,3))
output_dimensions = matrix(NA, nrow = dimensions , ncol = int)

# WARNING - This takes approximalty 50-90 minutes
index = 0
for (i in c(50,100,200,500,1000)) {
  index = index + 1
  interation_gd <- iterations_gd(number = 100, n = i, rho = 0.5) 
  interation_new <- iterations_new(number = 100, n = i, rho = 0.5) 
  print(index)
  output_dimensions[index,1] <- i
  output_dimensions[index,2] <- interation_gd$iterations
  output_dimensions[index,3] <-interation_new$iterations
  output_dimensions[index,4] <- interation_gd$time
  output_dimensions[index,5] <-interation_new$time
}

output_dimensions %>% 
  as.data.frame() %>% 
  ggplot() +
  geom_line(aes(x = V1, y = V2, colour = "Gradient Descent"), size = 1.5) +
  geom_line(aes(x = V1, y = V3, colour = "Newton's Method"), size = 1.5) +
  theme_grey() +
  scale_colour_manual("", values = c("Newton's Method"="#ff6666", 
                                     "Gradient Descent"="#5addd0")) +
  theme(legend.text=element_text(size=12), 
        axis.text = element_text(hjust = 0.5, size = 12), 
        axis.title=element_text(size=12)) +
  labs(x = "Dimensions (n)", y = "Iterations")

# save plot
ggsave("reports/images/dimensions.png", plot = last_plot(),
       scale = 1, dpi = 300, limitsize = TRUE)

output_dimensions %>% 
  as.data.frame() %>% 
  ggplot() +
  geom_line(aes(x = V1, y = V4, colour = "Gradient Descent"), size = 1.5) +
  geom_line(aes(x = V1, y = V5, colour = "Newton's Method"), size = 1.5) +
  theme_grey() +
  scale_colour_manual("", values = c("Newton's Method"="#ff6666", 
                                     "Gradient Descent"="#5addd0")) +
  theme(legend.text=element_text(size=12), 
        axis.text = element_text(hjust = 0.5, size = 12), 
        axis.title=element_text(size=12)) +
  labs(x = "Dimensions (n)", y = "Runtime (Seconds)")

# save plot
ggsave("reports/images/time.png", plot = last_plot(),
       scale = 1, dpi = 300, limitsize = TRUE)









