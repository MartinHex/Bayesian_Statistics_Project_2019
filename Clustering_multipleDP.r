path = getwd()
DataSetPath = paste(path, "/MyCSV/DataSets")
WorkSpacePath = paste(path,"MCMC_output")
setwd(DataSetPath)

rm(list=ls())

input = read.table("Combined_Lab_Well_Age_Sex_Onset_v2.csv",header=T, sep = ',')
names(input)

summary(input)

index = input$index + 1
nbrCovariates = dim(input)[2] - 3
n.data=dim(input)[1] # number of patients.
J=n.data
M=35 # Truncation for DP

library(rjags)      # JAGS interface for R
library(plotrix)    # To plot the CIs
set.seed(1)

nbrBins = max(input$Death_Bin)

data <- list(Y=input$Death_Bin, Other_Onset=input$Other_Onset, Bulbar_Onset=input$Bulbar_Onset,
             Limb_Onset=input$Limb_Onset, ALSFRS=input$AlsfrsScore, FVC = input$FVCScore,
             Age = input$Age, Sex = input$Sex, Creatine=input$CREATININE, Hemoglobin = input$HEMOGLOBIN,
             n.data=n.data, nbrCovariates=nbrCovariates, nbrBins = nbrBins, M = M, J = J, index = index)


inits = list(beta=matrix(0,nrow=nbrCovariates,ncol=nbrBins),
             .RNG.seed = 2,
             .RNG.name = 'base::Wichmann-Hill')

#-------------------  Create the JAGS model.  --------------------------
setwd(path)

modelGLMM_DP=jags.model("GLMM_multipleDP.bug",data=data,inits=inits,n.adapt=2000,n.chains=1)
update(modelGLMM_DP,2000)

variable.names = c("bb", "beta", "alpha")
n.iter=5000
thin=2

outGLMM_DP=coda.samples(model=modelGLMM_DP,variable.names=variable.names,n.iter=n.iter,thin=thin)

setwd(WorkSpacePath)
save(outGLMM_DP,file='Clustering_Combined_MultiDP_SmallerSet.rdata')


# -------------------------------- Potentially load data -----------------------------------
rm(list=ls())
# NOTE: file not on Github due to filesize limits, need to run locally to obtain the workspace.
load('Clustering_Combined_MultiDP_SmallerSet.Rdata')

# -------------------------------- Visualize data -------------------------------------

data=as.matrix(outGLMM_DP)
data=data.frame(data)
attach(data)
n.chain=dim(data)[1]

# -------------------------------- From Henrik, automated visualization ------------------------------------
# --------------------------------- Open a lot of windows for plotting -------------------------------------
B = nbrBins
nrCovars = nbrCovariates
beta_str = matrix(nrow = nrCovars, ncol = B)
for (i in 1:nrCovars) {
  beta_str[i, ] = paste('beta.', i, '.', 1:B, '.', sep = '')
}

cov_names = c('Limb_Onset', 'FVC_Score', 'AGE', 'SEX', 'CREATINE', 'HEMOGLOBIN')
y_strN = paste(rep('beta: ', 6), cov_names, sep='')


##################
### TRACEPLOTS ###
##################

# plot all traceplots for each covariate and death bin

# for (c in 1:nrCovars) {
#   x11()
#   par(mfrow=c(2,5))
#   for (b in 11:20) {
#     plot(out_data[, beta_str[c, b]], type='l', ylab = beta_str[c, b], xlab = "Iterations")
#   }
# }

# plot all traceplots for each bin and all covariates

for(b in 1:B){
  t_name = paste('Traceplot_bin ', b, sep = '')
  x11(title = t_name)
  par(mfrow=c(2, 3))
  for (c in 1:nrCovars) {
    plot(data[, beta_str[c, b]], type='l', ylab = y_strN[c], xlab = "Iterations")
  }
  mtext(paste('Bin ', b, sep = ''), side = 3, line = -3, outer = TRUE, font=2)
}

graphics.off()

############################
## Autocorrelation plots ###
############################

# plot all ACF for each covariate and death bin

# for (c in 1:nrCovars) {
#   x11()
#   par(mfrow=c(2,3))
#   for (b in 1:1) {
#     acf(out_data[, beta_str[c, 1]], lwd=3,col="red3",main=beta_str[c, b])
#   }
# }

# plot all ACF for each bin and all covariates

for(b in 1:B){
  x11(title = paste('Autocorrelation_bin ', b, sep = ''))
  par(mfrow=c(2, 3))
  for (c in 1:nrCovars) {
    acf(data[, beta_str[c, b]],lwd=3, col="red3", main=y_strN[c])
  }
  mtext(paste('Bin ', b, sep = ''), side = 3, line = -1.3, outer = TRUE, font=2)
}

graphics.off()

#################
# Density plots#
################

# plot all densities for each covariate and death bin

# for (c in 1:nrCovars) {
#   x11()
#   par(mfrow=c(2,5))
#   for (b in 1:10) {
#     plot(density(out_data[, beta_str[c, b]]), main = toString(c))#beta_str[c, b], xlab = NA, ylab = NA)
#   }
# }

# plot all densities for each bin and all covariates

for(b in 1:B){
  x11(title = paste('Density_function_bin ', b, sep = ''))
  par(mfrow=c(2, 3))
  for (c in 1:nrCovars) {
    plot(density(data[, beta_str[c, b]]), main=y_strN[c])
  }
  mtext(paste('Bin ', b, sep = ''), side = 3, line = -1.5, outer = TRUE, font=2)
}

graphics.off()


# Mean
for (b in 1:B) {
  for (c in 1:nrCovars) {
    print(beta_str[c, b])
    print(mean(data[, beta_str[c, b]]<0))
  }
}

#Plots the traceplot of bb for eachpatient for bin 14.

#nrPat = J
#for (i in 1:nrPat) {
#  if(i%%12==0) {
#    x11()
#    par(mfrow=c(2,6)) # some of the traceplots from b - values from the DP
#  }
#  bb_str = paste('bb.', i, '.','14','.', sep = '')
#  b_str = paste('b', i, sep = '')
#  plot(data[, bb_str], main = b_str, type = 'l')
#}

graphics.off()



# ------------------------------ Cluster wrt Binder's loss -------------------------------------


label.mat = list()
for (i in 0:19){
  start = 21 + J*i
  end = 490 + J*i
  label.mat[[i+1]] = as.matrix(data[,start:end]) # extract cluster labels for each bin
}

#First Binder's loss
sest.mat = matrix(0,nrow=nbrBins,ncol=J)
for (k in 1:nbrBins) {
  m=J
  G=n.chain
  pihat <- matrix(0,ncol=m,nrow=m)
  for(i in 1:G){
    ss <- label.mat[[k]][i,]
    cij <- outer(ss,ss,'==')
    pihat <- pihat+cij
  }

  pihat <- pihat/G

  #Binder loss function
  FF <- vector("numeric")
  K <- 0.9
  for(i in 1:G){
    ss <- label.mat[[k]][i,]
    cij <- outer(ss,ss,'==')
    pluto <- (pihat-K)*as.matrix(cij)
    pluto <-  pluto[upper.tri(pluto)]
    FF[i] <- sum(pluto)
  }

  plot(FF)

  ind.bind <- which.max(FF)[1]
  label.mat[[k]][ind.bind,]
  ll.bind <- label.mat[[k]][ind.bind,]
  unici <- unique(ll.bind)
  l.uni <- length(unici)

  ncl=l.uni

  #Find partition for specified bin
  Sest=rep(0,J)
  for(i in 1:ncl){
    Sest[as.numeric(which(ll.bind==unici[i]))] = i
  }

  #Combine all partitions for all bins to a matrix
  for (i in 1:J){
    sest.mat[k,i] = Sest[i]
  }
}

#Second Binder's loss
m=J
G=nbrBins
pihat <- matrix(0,ncol=m,nrow=m)
for(i in 1:G){
  ss <- sest.mat[i,]
  cij <- outer(ss,ss,'==')
  pihat <- pihat+cij
}

pihat <- pihat/G

#Binder loss function
FF <- vector("numeric")
K <- 0.1
for(i in 1:G){
  ss <- sest.mat[i,]
  cij <- outer(ss,ss,'==')
  pluto <- (pihat-K)*as.matrix(cij)
  pluto <-  pluto[upper.tri(pluto)]
  FF[i] <- sum(pluto)
}

ind.bind <- which.max(FF)[1]
ll.bind <- sest.mat[ind.bind,]
unici <- unique(ll.bind)
l.uni <- length(unici)

ncl=l.uni
Sest2=rep(0,J)
for(i in 1:ncl){
  Sest2[as.numeric(which(ll.bind==unici[i]))] = i
}

table(Sest2)

setwd(DataSetPath)
write.csv(Sest2, file ="Partition_multipleDP.csv", row.names=FALSE)








nclusLG=length(unique(Sest2))

bins = as.numeric(names(table(Sest2)))
freqs = as.vector(table(Sest2))

label=rep(0,J)

mylist=list()

for(i in 1:nclusLG)
{
  in_gr=which(Sest2==bins[i])
  label[in_gr]=i
  mylist[[i]]=in_gr
}
label
mylist
