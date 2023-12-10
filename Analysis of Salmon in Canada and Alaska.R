library("numDeriv")



#Datos--------------------------------------------------------------------------
Can_Fresh = c(129,
      148,
      179,
      156,
      140,
      108,
      135,
      152,
      153,
      152,
      148,
      145,
      123,
      117,
      118,
      153,
      154,
      155,
      128,
      144,
      133,
      128,
      144,
      125,
      153,
      108)

Can_Marine = c(420,
371,
407,
419,
362,
330,
355,
301,
397,
301,
383,
337,
364,
355,
379,
403,
390,
349,
400,
403,
375,
383,
373,
346,
352,
339)

Alas_Fresh = c(131,
               105,
               99,
               94,
               99,
               114,
               123,
               104,
               119,
               114,
               109,
               82,
               105,
               121,
               85,
               83,
               53,
               95,
               76,
               95,
               70,
               74,
               80,
               95,
               99,
               87)
Alas_marine = c(355,
                469,
                402,
                440,
                403,
                428,
                372,
                407,
                474,
                396,
                397,
                431,
                388,
                403,
                451,
                453,
                427,
                411,
                442,
                426,
                397,
                451,
                398,
                433,
                481,
                480)
x = Can_Marine-Can_Fresh
x = sort(x)


rpmu=function(mu,t1,t2,n){
  HH=t2-(t1^2)/n
  muhat=t1/n
  Hmu=t2-(2*t1*mu)+(n*(mu^2))
  rp=-(n/2)*(log(Hmu)-log(HH))
  return(rp)} 

sinnombre=function(object1){
  object2=object1
  names(object2)=NULL
  return(object2)}

lvnor=function(vec2,t1,t2,n){
  mu=vec2[1]
  sigma=vec2[2]
  if (sigma >0){
    s2=sigma^2
    lv=-n*log(sigma)-(t2/(2*s2))+(mu*t1/s2)-((n*mu^2)/(2*s2))
  }else{
    lv=-99999999999999999}
  return(lv)}

lpmu=function(vec1,t1,t2,n,mufix){
  mu=mufix
  sigma=vec1
  if (sigma >0){
    s2=sigma^2
    lv=-n*log(sigma) - (t2-(2*t1*mu) + (n*(mu^2)))/(2*s2)
  }else{
    lv=-99999999999999999
  }
  return(lv)
} 




#Estadísticas Descriptivas------------------------------------------------------

n = length(x)
med = mean(x)
std = sqrt(var(x)*(n-1)/n)


#Distribuciones-----------------------------------------------------------------

#Estadísticas suficientes-------------------------------------------------------

t1 = sum(x)
t2 = sum(x^2)
t3 = sum(log(x))

#Gamma--------------------------------------------------------------------------

#Función de Log-Verosimilitud con censura, a es el parámetro de forma, b el de escala

log_likeg = function(par) {
  a = par[1]
  b = par[2]
  vec_dist = (1:n) 
  for(i in (1:n)) {
     vec_dist[i] = log(pgamma(x[i]+0.5,shape = a, scale = b)-pgamma(x[i]-0.5,shape=a, scale=b))
   }
  return(sum(vec_dist))  
}

#Aproximaciones dadas por el método de momentos
a_prox = med^2/std^2
b_prox = std^2/med

#Se optimiza la función para encontrar los parámetros a y b
ab = optim(c(a_prox,b_prox), log_likeg, control = list(fnscale = -1, reltol = 1e-15))

#Se verifica si hay un problema con la convergencia (si no se imprime no hay problema)
if(ab$convergence != 0 && ab$convergence != 1) {
  print("Warning")
}
a = ab$par[1]
b = ab$par[2]

#mediana de la gamma
median_g = qgamma(0.5,shape = a, scale = b)

#Información de Akaike

AIC_g = 4 - 2*ab$value

#Promedio de distancias verticales ponderadas

u_g = pgamma(x, shape = a, scale = b)

k = c(1:n)

V_ug = k*(n+1-k)/((n+1)^2*(n+2))
PDV_g = (1/n)*sum(abs(u_g-(k/(n+1)))/sqrt(V_ug))

#Exponencial-------------------------------------------------------------------

log_likee = function(par) {
  vec_dist = (1:n) 
  for(i in (1:n)) {
    vec_dist[i] = log(pexp(x[i]+0.5, rate = 1/par)-pexp(x[i]-0.5,rate = 1/par))
  }
  return(sum(vec_dist))  
}

#Aproximación dada por el método de momentos
theta_aprox = med

#Lista de optimize
theta_l = optimize(log_likee, c(theta_aprox-2,theta_aprox+2), 
                 lower = theta_aprox-2, upper = theta_aprox+2, tol = 1e-8
                 ,maximum = TRUE)

theta = theta_l$maximum
mediam_e = qexp(0.5, rate = 1/theta)

#Información de Akaike

AIC_e = 2-2*theta_l$objective

#

#Promedio de distancias verticales ponderadas

u_e = pexp(x, rate = 1/theta)
u_e = sort(u_e)

k = c(1:n)

V_ue = k*(n+1-k)/((n+1)^2*(n+2))
PDV_e = (1/n)*sum(abs(u_e-k/(n+1))/sqrt(V_ue))

#Normal------------------------------------------------------------------------

#Logverosimilitud normal con censura
log_liken = function(par) {
  mul = par[1]
  sigl = par[2]
  vec_dist = (1:n) 
  for(i in (1:n)) {
    vec_dist[i] = log(pnorm(x[i]+0.5, mul, sigl)-pnorm(x[i]-0.5, mul, sigl))
  }
  return(sum(vec_dist))  
}

musig = optim(c(med,std), log_liken, control = list(fnscale = -1, reltol = 1e-15))

#Se verifica si hay un problema con la convergencia (si no se imprime no hay problema)
if(musig$convergence != 0 && musig$convergence != 1) {
  print("Warning")
}

mu = musig$par[1]
sigma = musig$par[2]

median_n = qnorm(0.5, mu, sigma)

#Información de Akaike

AIC_n = 4 - 2*musig$value

#Promedio de distancias verticales ponderadas

u_n = pnorm(x, mu, sigma)
u_n = sort(u_n)

k = c(1:n)

V_un = k*(n+1-k)/((n+1)^2*(n+2))
PDV_n = (1/n)*sum(abs(u_n-k/(n+1))/sqrt(V_un))

#Lognormal----------------------------------------------------------------------

#Logverosimilitud lognormal con censura

log_likeln = function(par) {
  mul = par[1]
  sigl = par[2]
  vec_dist = (1:n) 
  for(i in (1:n)) {
    vec_dist[i] = log(plnorm(x[i]+0.5, mul, sigl)-plnorm(x[i]-0.5, mul, sigl))
  }
  return(sum(vec_dist))  
}

#Aproximaciones dada por el método de momentos
mu_l = t3/n
sigma_l = sqrt(var(log(x))*(n-1)/n)

musigl = optim(c(mu_l,sigma_l), log_likeln, control = list(fnscale = -1, reltol = 1e-15))

#Se verifica si hay un problema con la convergencia (si no se imprime no hay problema)
if(musigl$convergence != 0 && musigl$convergence != 1) {
  print("Warning")
}

mu_l = musigl$par[1]
sig_l = musigl$par[2]

median_log = qlnorm(0.5, mu_l, sigma_l)

#Información de Akaike
T_1 = sum(log(x))
T_2 = 2*sum(log(x))
AIC_ln = 4-2*musigl$value

#Promedio de distancias verticales ponderadas

u_ln = plnorm(x, mu_l, sigma_l)
u_ln = sort(u_ln)

k = c(1:n)

V_uln = k*(n+1-k)/((n+1)^2*(n+2))
PDV_ln = (1/n)*sum(abs(u_ln-k/(n+1))/sqrt(V_uln))


#Weibull-------------------------------------------------------------------








