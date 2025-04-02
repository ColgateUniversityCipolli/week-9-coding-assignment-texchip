
library(tidyverse)

dat.precip <- read_csv(file = "agacis.csv")

dat.precip.long <- dat.precip |>    
  dplyr::select(-Annual) |>                   
  pivot_longer(cols = c(Jan, Feb, Mar, Apr,   
                        May, Jun, Jul, Aug, 
                        Sep, Oct, Nov, Dec), 
               values_to = "Precipitation",
               names_to = "Month") |>         
  mutate(Precipitation = case_when(Precipitation == "M" ~ NA_character_,
                                   TRUE                 ~ Precipitation))|>
  mutate(Precipitation = as.numeric(Precipitation))

#####################################
# MLEs (Weibull)
#####################################
llweibull <- function(par, data, neg=F){
  a <- exp(par[1]) 
  sigma <- exp(par[2]) 
  
  ll <- sum(log(dweibull(x=data, shape=a, scale=sigma)), na.rm=T)
  
  return(ifelse(neg, -ll, ll))
}

MLE.weibull <- optim(fn = llweibull,
              par = c(1,1),
              data = dat.precip.long$Precipitation,
              neg=T)

(MLE.weibull$par <- exp(MLE.weibull$par)) 
ll.weibull <- MLE.weibull$value

#####################################
# a.) MLEs (Gamma)
#####################################
llgamma <- function(par, data, neg=F){
  alpha <- exp(par[1])
  beta <- exp(par[2])
  
  ll <- sum(log(dgamma(x=data, shape=alpha, scale=beta)), na.rm=T)
  
  return(ifelse(neg, -ll, ll))
}

MLE.gamma <- optim(fn = llgamma,
                    par = c(1,1),
                    data = dat.precip.long$Precipitation,
                    neg=T)

(MLE.gamma$par <- exp(MLE.gamma$par))
ll.gamma <- MLE.gamma$value

#####################################
# b.) MLEs (Lognormal)
#####################################
lllognormal <- function(par, data, neg=F){
  mu <- par[1] 
  sigma <- exp(par[2]) 
  
  ll <- sum(log(dlnorm(x=data, meanlog=mu, sdlog=sigma)), na.rm=T)
  
  return(ifelse(neg, -ll, ll))
}

MLE.lognormal <- optim(fn = lllognormal,
                        par = c(0,1),
                        data = dat.precip.long$Precipitation,
                        neg=T)

MLE.lognormal$par[2] <- exp(MLE.lognormal$par[2]) 
(MLE.lognormal$par)
ll.lognormal <- MLE.lognormal$value

#####################################
# c.) Weibull / Gamma
#####################################

(Q.wg <- exp(ll.weibull - ll.gamma))
# > 1 => Weibull > Gamma

#####################################
# d.) Weibull / Log-Normal
#####################################

(Q.wl <- exp(ll.weibull - ll.lognormal))
# < 1 => Weibull < Log-Normal

#####################################
# e.) Gamma / Log-Normal
#####################################

(Q.gl <- exp(ll.gamma - ll.lognormal))
# < 1 => Gamma < Log-Normal
