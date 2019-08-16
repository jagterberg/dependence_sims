library(Matrix)
library(irlba)
if (!require(ggplot2)) {
  install.packages("ggplot2")
  library(ggplot2)
}
#if(!require(pushoverr)) {
#  install.packages("pushoverr")
#  library(pushoverr)
#}

library(Rcpp)

simulate_Erdos_renyi <- function(n,p) {
  A <- Matrix(0,n,n)
  A[upper.tri(A,diag = TRUE)] <- rbinom(choose(n,2) + n,1,p)
  A <- forceSymmetric(A,"U")
  return(A)
}

title <- "omega = I(both are 1)"

Rcpp::sourceCpp("generate_corr_sbm.cpp")
print("Sourced C++ Code, beginning simulations...")

set.seed(1234)

ns <- seq(500,8000,500)
#n <- 4000
p <- .05
a <- .03
b <- .07
corr = .4

sims <- 10

vech1 <- rep(0,length(ns))
vech2 <- rep(0,length(ns))
j <- 1
for (n in ns) {
  
  print(paste0("simulations for n = ",n))
  for (i in 1:sims) {
    
    A1 <- simulate_Erdos_renyi(n,p)
    A2 <- Matrix(simulate_corr_SBM(a,b,as.matrix(A1),corr,p))
    Mhat_naive <- bdiag(A1,A2)
    uhat_naive <- irlba(Mhat_naive,1,1)
    uhat_naive <- uhat_naive$u
    
    #Mhat_better <- Mhat_naive
    omega <- rowSums(A1 == 1 & A2==1)
    diag(Mhat_naive[c( (n+1):(2*n)), c((n+1):(2*n))]) <- omega
    uhat_better <- irlba(Mhat_naive,1,1)
    uhat_better <- uhat_better$u
    
    vech1[j] <- vech1[j] + min(min(abs(sqrt(2*n)*uhat_better - 1),min(abs(sqrt(2*n)*uhat_better + 1))))
    vech2[j] <- vech2[j] + min(min(abs(sqrt(2*n)*uhat_naive - 1),min(abs(sqrt(2*n)*uhat_naive + 1))))
  }
  vech1[j] <- vech1[j]/sims 
  vech2[j] <- vech2[j]/sims
  j <- j+1
  
}




dat <- data.frame(naive = vech2,off_diag =vech1,n = ns)

save(dat,file = paste0(title,"_output.RData"))
  


#cor(as.numeric(A1[upper.tri(A1,diag=TRUE)]),as.numeric(A2[upper.tri(A2,diag=TRUE)]))



#load(paste0(title,"_output.RData"))
jpeg(paste0(title,'.jpg'))
g <- ggplot(dat, aes(x = n))
g +   geom_line(aes(y = naive, linetype = 'Variable Name A'),lwd = 1) +
  geom_line(aes(y = off_diag, linetype = 'Variable Name B'),lwd = 1)  +
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





#message <- paste0(title,' simulations finished.')
#pushover(message=message, user='u89ff6gbp38dg2key3p5q8x5qa5ndu', app='ajg6age36i9jnqxwbaw7izjy845hd3')
  
  
  
  
  