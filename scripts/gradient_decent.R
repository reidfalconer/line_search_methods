# define function
defined_function <- function(x) {
  (1/2)*t((x - m))%*%A%*%(x - m) - log((x)^2) 
}
  