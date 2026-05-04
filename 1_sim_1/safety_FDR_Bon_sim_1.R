
library(MASS)
library(Matrix)
library(doParallel)
setwd("~/1_doc/Research/AbbVie/Safety/1_code/v2/1_sim_1/")

#########################################################################
alpha = 0.05
n.v = 15  ## number of variables per region or area
n.itt = 10^5  ## number of simulation iterations
n.cluster = 8  ## number of clusters for parallel computing 

###################################################################
for (n.r in c(5, 10)){ ## number of regions or areas

final.out.mat = matrix(NA, nrow = 28, ncol = n.r+5)

for (scen.ind in 1:28){

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
    mean.stats.vec = c(6, rep(0, n.v-1))  
  }
  
  ## number of regions under alternative hypothesis
  if (h1.ind==1) n.r.false = 0 
  if (h1.ind==2) n.r.false = 1 
  if (h1.ind==3) n.r.false = 2 
  
  mean.stats.false = 3
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
    
    ## the 2nd smallest p-values
    p.2nd.bon.vec = sapply(1:n.r, function(x.r){
      sort.temp = as.numeric(data.pvalue.temp[x.r, ])
      index.2nd = (order(sort.temp))[2]
      p.adj.temp = p.adjust(sort.temp, method = "holm")
      return(p.adj.temp[index.2nd])
    })
    
    # p.2nd.vec = sapply(1:n.r, function(x.r){
    #   # index.2nd = rank(data.pvalue.temp[x.r, ])[2]
    #   # p.adj.temp = p.adjust(data.pvalue.temp[x.r, ], method = "bonferroni")
    #   # return(p.adj.temp[index.2nd])
    #   (sort(data.pvalue.temp[x.r, ]))[2]
    # })
    # 
    # ## Bonferroni adjusted pvalues across n.v variables
    # p.2nd.bon.vec = pmin(1, p.2nd.vec*n.v)
    
    ## FDR control across n.r regions
    p.FDR.vec = p.adjust(p.2nd.bon.vec, method = "BH")
    
    ## summary
    # out.bon.mat[itt, ] = p.2nd.bon.vec<=alpha
  
    if (sum(p.FDR.vec<=alpha)>0){
      FDR.out.temp = sum(p.FDR.vec[(n.r.false+1):n.r]<=alpha)/sum(p.FDR.vec<=alpha)
      # out.FDR.overall.vec[itt] = sum(p.FDR.vec[(n.r.false+1):n.r]<=alpha)/sum(p.FDR.vec<=alpha)
    } else {
      FDR.out.temp = 0
      # out.FDR.overall.vec[itt] = 0
    }
  
    return(c(
      p.2nd.bon.vec<=alpha,
      FDR.out.temp
    ))
  
  }
  stopCluster(cl)
  
  output.mat.temp = matrix(unlist(output.vec.temp),
                           nrow = n.itt, ncol = n.r+1, byrow = TRUE)
  
  final.out.mat[scen.ind, ] = c(mean.ind, n.r.false, cor.stats.num.r, cor.stats.num.v,
                              apply(output.mat.temp, 2, mean))

}

print(final.out.mat)

write.csv(final.out.mat, paste0("safe_nr_", n.r, ".csv"))
}

stop("a")

##############################################################
## latex table, m = 5
library(xtable)

data.output = data.frame(read.csv("safe_nr_5.csv")[,-1])

data.output$V3 = sprintf('%.3f', data.output$V3)
data.output$V4 = sprintf('%.3f', data.output$V4)

data.output$V5 = paste0(sprintf('%.1f', data.output$V5*100), "%")
data.output$V6 = paste0(sprintf('%.1f', data.output$V6*100), "%")
data.output$V7 = paste0(sprintf('%.1f', data.output$V7*100), "%")
data.output$V8 = paste0(sprintf('%.1f', data.output$V8*100), "%")
data.output$V9 = paste0(sprintf('%.1f', data.output$V9*100), "%")
data.output$V10 = paste0(sprintf('%.1f', data.output$V10*100), "%")

data.output[data.output=="0.0%"] = "<0.1%"
data.output$V3[data.output$V3=="0.000"] = "0"
data.output$V3[data.output$V3=="0.700"] = "0.7"
data.output$V3[data.output$V3=="-0.014"] = "-1/74"

data.output$V4[data.output$V4=="0.000"] = "0"
data.output$V4[data.output$V4=="0.700"] = "0.7"
data.output$V4[data.output$V4=="-0.014"] = "-1/74"

print(xtable(data.output),include.rownames = FALSE)

##############################################################
## latex table, m = 10
library(xtable)

data.output = data.frame(read.csv("safe_nr_10.csv")[,-c(1, 6:15)])

data.output$V3 = sprintf('%.3f', data.output$V3)
data.output$V4 = sprintf('%.3f', data.output$V4)

data.output$V15 = paste0(sprintf('%.1f', data.output$V15*100), "%")

data.output[data.output=="0.0%"] = "<0.1%"
data.output$V3[data.output$V3=="0.000"] = "0"
data.output$V3[data.output$V3=="0.700"] = "0.7"
data.output$V3[data.output$V3=="-0.007"] = "-1/149"

data.output$V4[data.output$V4=="0.000"] = "0"
data.output$V4[data.output$V4=="0.700"] = "0.7"
data.output$V4[data.output$V4=="-0.007"] = "-1/149"

data.output$V1 = c("M1", rep("", 6),
                   "M2", rep("", 6),
                   "M2", rep("", 6),
                   "M2", rep("", 6))
data.output$V2 = c("S1", rep("", 6),
                   "S1", rep("", 6),
                   "S2", rep("", 6),
                   "S3", rep("", 6))

print(xtable(data.output),include.rownames = FALSE)











