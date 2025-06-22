rm(list=ls())
library(dplyr, warn.conflicts = FALSE)
library(boot)

# Ex1

# Declaration des variables et fonctions
x <- c(1274,4860,2376,10567,526,1323,1890,2307,3042,1999,5106,833)
z <- c(913,510,1224,63,326,305,486,213,1813,843,1826,1826)
y <- c(2187,5370,3600,10630,852,1628,2376,2520,4855,2842,NA,NA)

lambda.k <- function(k){
  if (k == 1){
    return(0.0001)
  } else{
    return(12 / (sum(x+z) + (12 - 2)/lambda.k(k-1))) # sum(y[1:10]) + sum(x[11]+x[12]+z[11]+z[12])
  }
}

Qfunction <- function(lambda,k){
  return (12*log(lambda)-lambda*(sum(x+z) + (12 - 2)/lambda.k(k)))
}

# Etape Maximisation
lambda.vector <- c()
lambda.vector[1] <- lambda.k(1)
lambda.vector[2] <- lambda.k(2)
k <- 2 
while (abs(lambda.k(k) - lambda.k(k-1))/lambda.k(k-1)>0.000015){
  lambda.vector[k+1] = lambda.k(k+1)
  k <- k + 1
}
# Return the last element of the vector
print("The parameter lambda for this exponential law is:")
print(tail(lambda.vector, n=1))

# Ex 2
set.seed(10)
N = 1000
n = 20
# Ex 2
shape <- 5 
scale <- 1/2 

#creation d'une liste de 1000 echantillons de taille "20 valeurs"
# l <- list()
# for (i in (1:N)){
#   l <- c(l,list(rgamma(n, shape = shape, scale = scale)))
# }
df <- data.frame(values=rgamma(n, shape = shape, scale = scale))
# echantillon boostrapé de l'estimateur theta.one (l'esperance des Xi)
function.mean <- function(df,i){
  d <- df[i,]
  return (mean(d))
}

# function.mean.test <- function(l,i){
#   return (mean(unlist(l[i])))
# }
# sample.bootstrap.one.test <- boot::boot(data = l,statistic = function.mean.test,R=1000)
#https://dimension.usherbrooke.ca/dimension/ssrBootstrapPropriete.html
# cf exemple mais ici notre daframe n'a qu'une seule colonne. IL ne faut pas la spécifier
sample.bootstrap.one <- boot::boot(data = df,statistic = function.mean,R=1000)
ci.one.perc <- boot::boot.ci(sample.bootstrap.one,type=c("perc"),conf = 0.9)
ci.one.bc <- boot::boot.ci(sample.bootstrap.one,type=c("basic"),conf = 0.9)
ci.one.bca <- boot::boot.ci(sample.bootstrap.one,type=c("bca"),conf = 0.9)

# echantillon boostrapé de l'estimateur theta.two (variance des Xi)
function.variance <- function(df,i){
  d <- df[i,]
  return (var(d))
}
sample.bootstrap.two <- boot::boot(data = df,statistic = function.variance,R=1000)
ci.two.perc <- boot::boot.ci(sample.bootstrap.two,type=c("perc"),conf = 0.9)
ci.two.bc <- boot::boot.ci(sample.bootstrap.two,type=c("basic"),conf = 0.9)
ci.two.bca <- boot::boot.ci(sample.bootstrap.two,type=c("bca"),conf = 0.9)

# echantillon boostrapé de l'estimateur theta.three (standard deviation des Xi)
function.sd <- function(df,i){
  d <- df[i,]
  return (sd(d))
}
sample.bootstrap.three <- boot::boot(data = df,statistic = function.sd,R=1000)
ci.three.perc <- boot::boot.ci(sample.bootstrap.three,type=c("perc"),conf = 0.9)
ci.three.bc <- boot::boot.ci(sample.bootstrap.three,type=c("basic"),conf = 0.9)
ci.three.bca <- boot::boot.ci(sample.bootstrap.three,type=c("bca"),conf = 0.9)

# ---------------------------------------------------------------------------
# first.sample <- unlist(l[1])
# # echantillon boostrapé de l'estimateur theta.one (l'esperance des Xi)
# function.mean.bis <- function(data,i){
#   return (mean(data[i]))
# }
# sample.bootstrap.one.bis <- boot::boot(data = first.sample,statistic = function.mean.bis,R=1000)
# ci.one.perc.bis <- boot::boot.ci(sample.bootstrap.one.bis,type=c("perc"),conf = 0.9)
# ci.one.bca.bis <- boot::boot.ci(sample.bootstrap.one.bis,type=c("bca"),conf = 0.9)
# borneInferieure=ci.one.bca.bis$bca[4]
# borneSuperieure=ci.one.bca.bis$bca[5]

# echantillon boostrapé de l'estimateur theta.four (cf formula question 5. ex2)
# function.formula.four <- function(l,i){
#   print(i)
#   list.elements <- unlist(l[i])
#   print('------------------')
#   print(list.elements)
#   print(length(list.elements))
#   # 26992 2.0588
#   list.elements.mean <- mean(list.elements)
#   numerator <- sum(sapply(list.elements, function(x){(x - list.elements.mean)^3})) / length(list.elements)
#   return (numerator/sd(list.elements)^3)
# }

function.formula.four <- function(df,i){
  d <- df[i,]
  df.mean <- mean(d)
  numerator <- sum(sapply(d, function(x){(x - df.mean)^3})) / length(d)
  return (numerator/sd(d)^3)
}
sample.bootstrap.four <- boot::boot(data = df,statistic = function.formula.four,R=1000)
ci.three.perc <- boot::boot.ci(sample.bootstrap.four,type=c("perc"),conf = 0.9)
ci.three.bc <- boot::boot.ci(sample.bootstrap.four,type=c("basic"),conf = 0.9)
ci.three.bca <- boot::boot.ci(sample.bootstrap.four,type=c("bca"),conf = 0.9)

# Optimiser l'intervalle perc,bc,bca en variant la taille de l'echantillon dans le boot::boot()
# Finding the R, number of boostrap replicates to fit the best CI
IC.values <- function(estimator,R){
  # inputs: estimator (function object) and R the number of sample
  # output: IC lists for our study (perc,bc,bca)
  
  # boostraped sample for a theta i with i in [1:4] and i integer. 
  sample.boostrap.generic <- boot::boot(data = df,statistic = estimator,R=R)
  # IC for estimator theta i with i in [1:4] and i integer. 
  ci.perc.generic <- boot::boot.ci(sample.boostrap.generic,type=c("perc"),conf = 0.9)
  ci.bc.generic <- boot::boot.ci(sample.boostrap.generic,type=c("basic"),conf = 0.9)
  ci.bca.generic <- boot::boot.ci(sample.boostrap.generic,type=c("bca"),conf = 0.9)
  return(list(ci.perc.generic$percent,ci.bc.generic$basic,ci.bca.generic$bca))
}

delta <- function(df){
  # input: data.frame containing the method 'conf' from the returned object boot.ci()
  # output: data.frame with a new column which is the interval.
  
  # so we take the V5 and V4 elements which are the bounds.
  df$IC <- df$V5-df$V4
  return(df)
}

best.IC.calculate <- function(estimator,R){
  IC.types <- c("perc", "bc", "bca")
  res.IC.matrix <- sapply(R,function(x){IC.values(estimator=estimator,x)})
  # delta entre borne sup et borne inf
  res.IC.matrix <- apply(res.IC.matrix,c(1,2),function(x){delta(data.frame(x))})
  
  ll <- list()
  i <- 1
  # on itere sur les lignes (les trois types d'intervalles)
  # ici notre etude est sur 3 types d'intervalles (on boucle sur une condition d'arrêt pour qu'il sorte a i<-4)
  while(i < 4) {
    # res.IC.matrix[i,]) binding the n elements on a row with n = length(R)
    # output data.frame: aggregation des elements par ligne puis ajout d'une colone type d'intervale IC et index, 
    # creation d'une liste de data.frame
    ll[[i]] <- cbind(bind_rows(res.IC.matrix[i,]),typeIC=IC.types[i],index=seq(1:length(R)))
    i <- i + 1
  }
  # https://github.com/tidyverse/dplyr/issues/4489
  ll <- suppressWarnings(dplyr::bind_rows(ll))
  min.ll <- ll %>% group_by(typeIC) %>% summarise(IC=min(IC))
  res <- merge(min.ll,ll,by=c("typeIC","IC"),all.x=TRUE) %>% mutate(index.value=sapply(index,function(x){R[x]}))
  return(res)
}


main <- function(){
  R <- seq(from=999,to=3000,by = 1000) # mettre b=10
  print("Calcul du  R adequat par IC pour trouver les meilleurs intervales de confiance pour theta 1")
  res.theta1 <- best.IC.calculate(estimator = function.mean,R=R)
  print(res.theta1)
  print("Calcul du R adequat par IC pour trouver les meilleurs intervales de confiance pour theta 2")
  res.theta2 <- best.IC.calculate(estimator = function.variance,R=R)
  print(res.theta2)
  print("Calcul du R adequat par IC pour trouver les meilleurs intervales de confiance pour theta 3")
  res.theta3 <- best.IC.calculate(estimator = function.sd,R=R)
  print(res.theta3)
  print("Calcul du R adequat par IC pour trouver les meilleurs intervales de confiance pour theta 4")
  res.theta4 <- best.IC.calculate(estimator = function.formula.four,R=R)
  print(res.theta4)
}
# Question 6
main()

# piste ameliorations
# plots et mclapply()





