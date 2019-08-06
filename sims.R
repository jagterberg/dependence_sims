library(Matrix)
library(irlba)
library(igraph)
library(Rcpp)

simulate_Erdos_renyi <- function(n,p) {
  A <- Matrix(0,n,n)
  A[upper.tri(A,diag = TRUE)] <- rbinom(choose(n,2) + n,1,p)
  A <- forceSymmetric(A,"U")
  return(A)
}

Rcpp::sourceCpp("generate_corr_sbm.cpp")


ns <- seq(1000,10000,500)
#n <- 4000
p <- .4
a <- .7
b <- .1
corr = .1

sims <- 1

vech1 <- rep(0,length(ns))
vech2 <- rep(0,length(ns))
j <- 1
for (n in ns) {
  omega <- log(n)
  
  for (i in 1:sims) {
    A1 <- simulate_Erdos_renyi(n,p)
    A2 <- Matrix(simulate_corr_SBM(a,b,as.matrix(A1),corr,p))
    Mhat_naive <- bdiag(A1,A2)
    uhat_naive <- irlba(Mhat_naive,1,1)
    uhat_naive <- uhat_naive$u
    
    #Mhat_better <- Mhat_naive
    diag(Mhat_naive[c( (n+1):(2*n)), c((n+1):(2*n))]) <- rep(omega,n)
    uhat_better <- irlba(Mhat_naive,1,1)
    uhat_better <- uhat_better$u
    
    vech1[j] <- vech1[j] + min(min(abs(sqrt(n)*uhat_better - 1),min(abs(sqrt(n)*uhat_better + 1))))
    vech2[j] <- vech2[j] + min(min(abs(sqrt(n)*uhat_naive - 1),min(abs(sqrt(n)*uhat_naive + 1))))
  }
  vech1[j] <- vech1[j]/sims
  vech2[j] <- vech2[j]/sims
  j <- j+1
  
}

library(ggplot2)

jpeg('rplot.jpg')
dat <- data.frame(naive = vech2,off_diag =vech2,n = ns)
g <- ggplot(dat, aes(x = n))
g +   geom_line(aes(y = naive, col = 'Variable Name A'),lwd = 1) +
  geom_line(aes(y = off_diag, col = 'Variable Name B'),lwd = 1)  +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(plot.title = element_text(size = 10,hjust = 0.5))+
  ggtitle( 'Estimated Eigenvectors Minus true eigenvectors error (up to sign)') +
  ylab('Error')  + xlab('n') +
  theme(axis.title.y = element_text(size = 10)) +
  theme(axis.title.x = element_text(size =10)) +
  scale_linetype_manual(values = c('dashed','solid'),
                        labels= c('naive estimage','padded off-diagonal')
                        ,name = "") 

dev.off()

  


#cor(as.numeric(A1[upper.tri(A1,diag=TRUE)]),as.numeric(A2[upper.tri(A2,diag=TRUE)]))











