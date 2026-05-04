
library(MASS)
library(Matrix)
library(doParallel)
setwd("~/1_doc/Research/AbbVie/Safety/1_code/v2/4_sim_4_alt/")

#########################################################################
alpha = 0.05
n.v = 15  ## number of variables per region or area
n.itt = 10^5  ## number of simulation iterations
n.cluster = 8  ## number of clusters for parallel computing 

###################################################################
for (n.r in c(5)){ ## number of regions or areas

final.out.mat = matrix(NA, nrow = 28, ncol = 3*n.r+7)

for (scen.ind in c(1:28)){

  print(scen.ind)
  
  if (scen.ind==1){mean.ind=1; h1.ind=1; cor.r.ind = 1; cor.v.ind = 1}
  if (scen.ind==2){mean.ind=1; h1.ind=1; cor.r.ind = 1; cor.v.ind = 2}
  if (scen.ind==3){mean.ind=1; h1.ind=1; cor.r.ind = 1; cor.v.ind = 3}
  if (scen.ind==4){mean.ind=1; h1.ind=1; cor.r.ind = 2; cor.v.ind = 1}
  if (scen.ind==5){mean.ind=1; h1.ind=1; cor.r.ind = 2; cor.v.ind = 2}
  if (scen.ind==6){mean.ind=1; h1.ind=1; cor.r.ind = 2; cor.v.ind = 3}
  if (scen.ind==7){mean.ind=1; h1.ind=1; cor.r.ind = 3; cor.v.ind = 3}
  
  if (scen.ind==8){mean.ind=2; h1.ind=1; cor.r.ind = 1; cor.v.ind = 1}
  if (scen.ind==9){mean.ind=2; h1.ind=1; cor.r.ind = 1; cor.v.ind = 2}
  if (scen.ind==10){mean.ind=2; h1.ind=1; cor.r.ind = 1; cor.v.ind = 3}
  if (scen.ind==11){mean.ind=2; h1.ind=1; cor.r.ind = 2; cor.v.ind = 1}
  if (scen.ind==12){mean.ind=2; h1.ind=1; cor.r.ind = 2; cor.v.ind = 2}
  if (scen.ind==13){mean.ind=2; h1.ind=1; cor.r.ind = 2; cor.v.ind = 3}
  if (scen.ind==14){mean.ind=2; h1.ind=1; cor.r.ind = 3; cor.v.ind = 3}

  if (scen.ind==15){mean.ind=2; h1.ind=2; cor.r.ind = 1; cor.v.ind = 1}
  if (scen.ind==16){mean.ind=2; h1.ind=2; cor.r.ind = 1; cor.v.ind = 2}
  if (scen.ind==17){mean.ind=2; h1.ind=2; cor.r.ind = 1; cor.v.ind = 3}
  if (scen.ind==18){mean.ind=2; h1.ind=2; cor.r.ind = 2; cor.v.ind = 1}
  if (scen.ind==19){mean.ind=2; h1.ind=2; cor.r.ind = 2; cor.v.ind = 2}
  if (scen.ind==20){mean.ind=2; h1.ind=2; cor.r.ind = 2; cor.v.ind = 3}
  if (scen.ind==21){mean.ind=2; h1.ind=2; cor.r.ind = 3; cor.v.ind = 3}
  
  if (scen.ind==22){mean.ind=2; h1.ind=3; cor.r.ind = 1; cor.v.ind = 1}
  if (scen.ind==23){mean.ind=2; h1.ind=3; cor.r.ind = 1; cor.v.ind = 2}
  if (scen.ind==24){mean.ind=2; h1.ind=3; cor.r.ind = 1; cor.v.ind = 3}
  if (scen.ind==25){mean.ind=2; h1.ind=3; cor.r.ind = 2; cor.v.ind = 1}
  if (scen.ind==26){mean.ind=2; h1.ind=3; cor.r.ind = 2; cor.v.ind = 2}
  if (scen.ind==27){mean.ind=2; h1.ind=3; cor.r.ind = 2; cor.v.ind = 3}
  if (scen.ind==28){mean.ind=2; h1.ind=3; cor.r.ind = 3; cor.v.ind = 3}
  
  
  ## mean vector for n.v variables per region or area
  if (mean.ind==1){
    mean.stats.vec = rep(0, n.v)    
  } else{
    mean.stats.vec = c(2, rep(0, n.v-1))  
  }
  
  ## number of regions under alternative hypothesis
  if (h1.ind==1) n.r.false = 0 
  if (h1.ind==2) n.r.false = 1 
  if (h1.ind==3) n.r.false = 2 
  
  mean.stats.false = 1.6
  mean.stats.all.vec = rep(mean.stats.vec, n.r) + 
    c(rep(mean.stats.false, n.r.false*n.v), rep(0, n.r*n.v - n.r.false*n.v))
  
  # print(mean.stats.all.vec)
  ## corr coef across the region
  if (cor.r.ind == 1) cor.stats.num.r = 0  
  if (cor.r.ind == 2) cor.stats.num.r = -1/(n.r*n.v-1)
  if (cor.r.ind == 3) cor.stats.num.r = 0.7  
  
  ## corr coef for variables within a region
  if (cor.v.ind == 1) cor.stats.num.v = 0  
  if (cor.v.ind == 2) cor.stats.num.v = -1/(n.r*n.v-1)
  if (cor.v.ind == 3) cor.stats.num.v = 0.7  
  
  ## cov matrix for the n.v variables
  stats.sigma.mat = matrix(cor.stats.num.v, nrow = n.v, ncol = n.v)
  diag(stats.sigma.mat) = 1
  
  ## cov matrix for n.v*n.r variables
  big.stats.sigma.mat = matrix(cor.stats.num.r, nrow = n.v*n.r, ncol = n.v*n.r)
  for (n.r.ind in 1:n.r){
    big.stats.sigma.mat[((1:n.v)+(n.r.ind-1)*n.v),((1:n.v)+(n.r.ind-1)*n.v)] = 
      stats.sigma.mat
  }

  ## output matrix
  # out.bon.mat = matrix(NA, nrow = n.itt, ncol = n.r)
  # out.FDR.overall.vec = rep(NA, n.itt)
  
  #####################################################################
  cl = makeCluster(n.cluster)
  registerDoParallel(cl)
  output.vec.temp = foreach(itt = 1:n.itt) %dopar% {
  
  # for (itt in 1:n.itt){
   
    library(MASS)
    library(Matrix)
    
    set.seed(itt+n.itt*scen.ind)
    
    # print(itt)
    
    data.stats.vec = mvrnorm(1, mu = mean.stats.all.vec, 
                             Sigma = big.stats.sigma.mat)
    
    ## matrix of test statistics; row is region, column is variable
    data.stats.temp = matrix(data.stats.vec,
                                nrow = n.r, ncol = n.v, byrow = TRUE)
    
    ## one-sided p-value
    data.pvalue.temp = 1-pnorm(data.stats.temp)
    
    ### proposed method
    ## the 2nd smallest p-values
    p.2nd.bon.vec = sapply(1:n.r, function(x.r){
      sort.temp = as.numeric(data.pvalue.temp[x.r, ])
      index.2nd = (order(sort.temp))[2]
      p.adj.temp = p.adjust(sort.temp, method = "holm")
      return(p.adj.temp[index.2nd])
    })
    
    ## FDR control across n.r regions
    p.FDR.vec = p.adjust(p.2nd.bon.vec, method = "BH")
  
    if (sum(p.FDR.vec<=alpha)>0){
      FDR.out.temp = sum(p.FDR.vec[(n.r.false+1):n.r]<=alpha)/sum(p.FDR.vec<=alpha)
    } else {
      FDR.out.temp = 0
    }
    
    ### proposed method with at least 3
    ## the 3rd smallest p-values
    p.3rd.bon.vec = sapply(1:n.r, function(x.r){
      sort.temp = as.numeric(data.pvalue.temp[x.r, ])
      index.2nd = (order(sort.temp))[3]
      p.adj.temp = p.adjust(sort.temp, method = "holm")
      return(p.adj.temp[index.2nd])
    })
    
    ## FDR control across n.r regions
    p.FDR.3rd.vec = p.adjust(p.3rd.bon.vec, method = "BH")
    
    if (sum(p.FDR.3rd.vec<=alpha)>0){
      FDR.out.3rd.temp = 
        sum(p.FDR.3rd.vec[(n.r.false+1):n.r]<=alpha)/sum(p.FDR.3rd.vec<=alpha)
    } else {
      FDR.out.3rd.temp = 0
    }
  
    ### proposed method with at least 4
    ## the 4th smallest p-values
    p.4th.bon.vec = sapply(1:n.r, function(x.r){
      sort.temp = as.numeric(data.pvalue.temp[x.r, ])
      index.2nd = (order(sort.temp))[4]
      p.adj.temp = p.adjust(sort.temp, method = "holm")
      return(p.adj.temp[index.2nd])
    })
    
    ## FDR control across n.r regions
    p.FDR.4th.vec = p.adjust(p.4th.bon.vec, method = "BH")
    
    if (sum(p.FDR.4th.vec<=alpha)>0){
      FDR.out.4th.temp = 
        sum(p.FDR.4th.vec[(n.r.false+1):n.r]<=alpha)/sum(p.FDR.4th.vec<=alpha)
    } else {
      FDR.out.4th.temp = 0
    }
    
    ## summary
    return(c(
      p.2nd.bon.vec<=alpha,
      FDR.out.temp,
      p.3rd.bon.vec<=alpha,
      FDR.out.3rd.temp,
      p.4th.bon.vec<=alpha,
      FDR.out.4th.temp
    ))
  
  }
  stopCluster(cl)
  
  output.mat.temp = matrix(unlist(output.vec.temp),
                           nrow = n.itt, ncol = 3*(n.r+1), byrow = TRUE)
  
  final.out.mat[scen.ind, ] = c(mean.ind, n.r.false, cor.stats.num.r, 
                                cor.stats.num.v,
                              apply(output.mat.temp, 2, mean))
  
  colnames(final.out.mat) = c("mean.ind", "n.r.false", "cor.stats.num.r", 
                              "cor.stats.num.v",
                              "v1", "v2", "v3", "v4", "v5", "FDR",
                              "v1.3rd", "v2.3rd", "v3.3rd", 
                              "v4.3rd", "v5.3rd", "FDR.3rd",
                              "v1.4th", "v2.4th", "v3.4th", 
                              "v4.4th", "v5.4th", "FDR.4th"
                              )

}

print(final.out.mat)

write.csv(final.out.mat, paste0("safe_alt_nr_", n.r, ".csv"))
}


##############################################################
## latex table, m = 5
library(xtable)

data.output.1 = data.frame(read.csv("safe_alt_nr_5.csv")[,-1])
data.output = data.frame("mu" = data.output.1$mean.ind,
                         "w" = data.output.1$n.r.false)

data.output$r = sprintf('%.3f', data.output.1$cor.stats.num.r)
data.output$v = sprintf('%.3f', data.output.1$cor.stats.num.v)

data.output$h1.m1 = paste0(sprintf('%.1f', data.output.1$v1*100), "%")
data.output$FDR.m1 = paste0(sprintf('%.1f', data.output.1$FDR*100), "%")

data.output$h1.m2 = paste0(sprintf('%.1f', data.output.1$v1.3rd*100), "%")
data.output$FDR.m2 = paste0(sprintf('%.1f', data.output.1$FDR.3rd*100), "%")

data.output$h1.m3 = paste0(sprintf('%.1f', data.output.1$v1.4th*100), "%")
data.output$FDR.m3 = paste0(sprintf('%.1f', data.output.1$FDR.4th*100), "%")

data.output[data.output=="0.0%"] = "<0.1%"
data.output$r[data.output$r=="0.000"] = "0"
data.output$r[data.output$r=="0.700"] = "0.7"
data.output$r[data.output$r=="-0.014"] = "-1/74"

data.output$v[data.output$v=="0.000"] = "0"
data.output$v[data.output$v=="0.700"] = "0.7"
data.output$v[data.output$v=="-0.014"] = "-1/74"

### additional modifications
data.output$mu = c("M1", rep("", 6),
                   "M2", rep("", 6),
                   "M2", rep("", 6),
                   "M2", rep("", 6))

data.output$w = c("S1", rep("", 6),
                   "S1", rep("", 6),
                   "S2", rep("", 6),
                   "S3", rep("", 6))

print(xtable(data.output),include.rownames = FALSE)


