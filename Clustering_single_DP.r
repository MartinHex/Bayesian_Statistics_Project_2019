setwd(".../your/path/MyCSV/DataSets")

rm(list=ls())

# Read input data, here structured as: 
# index, Limb_Onset, Death_Date, Death_Bin, FVCScore, Age, Sex, CREATINE and HEMOGLOBIN
input = read.table("Combined_Lab_Well_Age_Sex_Onset_v2.csv",header=T, sep = ',')
names(input)

summary(input)

nbrCovariates = dim(input)[2] - 3 # one column of indices and death_date+death_bin
n.data=dim(input)[1] # number of patients.
J=n.data  # number of units, our case the same as number of data points/patients
M=30 # truncation parameter in dirichlet process

library(rjags)      # JAGS interface for R
library(plotrix)    # To plot the CIs
set.seed(1)

#If the indexing starts for 0 we need to add 1 to allign, to port between python and R
input$index = input$index + 1
nbrBins = max(input$Death_Bin) # Number of discrete time steps
# Data for JAGS model
data <- list(Y=input$Death_Bin, Limb_Onset=input$Limb_Onset, FVC = input$FVCScore, age=input$Age,
             sex=input$Sex, creatine=input$CREATININE, hemoglobin=input$HEMOGLOBIN, index=input$index, n.data=n.data,
             J=J,M=M, nbrCovs=nbrCovariates, nbrBins=nbrBins)

#Initiate variables for JAGS
r=rep(0.5,M)
theta=rep(0,M)
S=rep(1,n.data)
Y = rep(0, n.data)

inits = list(
  beta=matrix(0,nrow=nbrCovariates,ncol=nbrBins),  # Regresstion parameter for all covariates and each bin
  lambda.bb = 1,
  r=r,
  theta= theta,
  S = S,
  .RNG.seed=2,
  .RNG.name = 'base::Wichmann-Hill'
)

#-------------------  Create the JAGS model.  --------------------------
setwd(".../your/path")

modelGLMM_DP=jags.model("GLMM_single_DP.bug",data=data,inits=inits,n.adapt=10000,n.chains=1)
update(modelGLMM_DP,10000)

variable.names=c("bb", "beta", "alpha")
n.iter=10000
thin=10

outGLMM_DP=coda.samples(model=modelGLMM_DP,variable.names=variable.names,n.iter=n.iter,thin=thin)

setwd(".../your/path/MCMC_output")
save(outGLMM_DP,file='GLMMDP_output_20bins_15_02_10000.rdata')



#####################
##### Results #######
#####################
rm(list=ls())

setwd(".../your/path/MyCSV/DataSets")
ALS_data = read.csv("Combined_Lab_Well_Age_Sex_Onset_v2.csv", header = T, sep=",") # load initial covariates
attach(ALS_data)
head(ALS_data)

nrPat = nrow(ALS_data)
J = nrow(ALS_data)
nrCovars = ncol(ALS_data)-3 # index, death_bin and death_date
B = max(ALS_data$Death_Bin) # number of discrete time steps

library(coda)
library(plotrix)
library(latex2exp)

setwd(".../your/path/MCMC_output")
load('GLMMDP_output_20bins_15_02_10000.rdata') # Load MCMC output workspace

out_data=as.matrix(outGLMM_DP)
out_data=data.frame(out_data)
attach(out_data)
n.chain=dim(out_data)[1]   # final sample size of the MCMC chain

beta_str = matrix(nrow = nrCovars, ncol = B) # Names of regression parameters in MCMC output
for (i in 1:nrCovars) {
  beta_str[i, ] = paste('beta.', i, '.', 1:B, '.', sep = '')
}

cov_names = c('Limb_Onset', 'FVC_Score', 'AGE', 'SEX', 'CREATINE', 'HEMOGLOBIN')
y_strN = paste(rep('beta: ', 6), cov_names, sep='') # For axis labels in plots


##################
### TRACEPLOTS ###
##################

# Plot all traceplots of beta for each bin and all covariates
for(b in 1:B){
  t_name = paste('Traceplot_bin ', b, sep = '')
  x11(title = t_name)
  par(mfrow=c(2, 3))
  for (c in 1:nrCovars) {
    plot(out_data[, beta_str[c, b]], type='l', ylab = y_strN[c], xlab = "Iterations")
  }
  mtext(paste('Bin ', b, sep = ''), side = 3, line = -3, outer = TRUE, font=2)
}

graphics.off()

# Traceplots of b - values from the DP
for (i in 1:nrPat) {
  if(i%%12==0) {
    x11()
    par(mfrow=c(2,6))
  }
  bb_str = paste('bb.', i, '.', sep = '')
  b_str = paste('b', i, sep = '')
  plot(out_data[, bb_str], main = b_str, type = 'l')
}

graphics.off()

############################
## Autocorrelation plots ###
############################

# Plot autocorrelation of beta for each bin and all covariates
for(b in 1:B){
  x11(title = paste('Autocorrelation_bin ', b, sep = ''))
  par(mfrow=c(2, 3))
  for (c in 1:nrCovars) {
    acf(out_data[, beta_str[c, b]],lwd=3, col="red3", main=y_strN[c])
  }
  mtext(paste('Bin ', b, sep = ''), side = 3, line = -1.3, outer = TRUE, font=2)
}

graphics.off()

#################
# Density plots#
################

# Plot all density distribution of beta for each bin and all covariates
for(b in 1:B){
  x11(title = paste('Density_function_bin ', b, sep = ''))
  par(mfrow=c(2, 3))
  for (c in 1:nrCovars) {
    plot(density(out_data[, beta_str[c, b]]), main=y_strN[c])
  }
  mtext(paste('Bin ', b, sep = ''), side = 3, line = -1.5, outer = TRUE, font=2)
}

graphics.off()


# Mean of beta in each bin
for (b in 1:B) {
  for (c in 1:nrCovars) {
    print(beta_str[c, b])
    print(mean(out_data[, beta_str[c, b]]<0))
  }
}


#################
# alpha and K_n #
#################
# Traceplots
x11()
par(mfrow=c(1,2))
plot(out_data[,'alpha'],main='total mass',type='l', ylab = NA) #if a prior is assigned for alpha
plot(out_data[,'nbrClusters'],main='K_n',type='l')

# Density
x11()
par(mfrow=c(1,2))
plot(density(out_data[,'alpha']),main='total mass', ylab = NA)
plot(table(out_data[,'nbrClusters']),main='K_n')

# Mean and Var
mean(alpha);var(alpha)
mean(out_data[,'nbrClusters']);var(out_data[,'nbrClusters'])



###################################
# Cluster estimation by Lau&Green #
###################################
label.mat = as.matrix(out_data[, 2:(nrPat+1)]) # Extract cluster labels from position in MCMC output

pihat <- matrix(0,ncol=J,nrow=J)
for(i in 1:n.chain){
  ss <- label.mat[i,]
  cij <- outer(ss,ss,FUN = '==')
  pihat <- pihat+cij
}
pihat <- pihat/n.chain

### Binder's loss function ###
FF <- vector("numeric")
K <- 0.5
for(i in 1:n.chain){
  ss <- label.mat[i,]
  cij <- outer(ss,ss,'==')
  pluto <- (pihat-K)*as.matrix(cij)
  pluto <-  pluto[upper.tri(pluto)]
  FF[i] <- sum(pluto)
}

plot(FF)

ind.bind <- which.max(FF)[1]
label.mat[ind.bind,]
plot(FF)
ll.bind <- label.mat[ind.bind,]
unici <- unique(ll.bind)
l.uni <- length(unici) # estimated number of groups

ncl=l.uni
for(i in 1:ncl){
  print(as.numeric(which(ll.bind==unici[i])))
}

Sest=rep(0,J)
for(i in 1:ncl){
  Sest[as.numeric(which(ll.bind==unici[i]))] = i
}

Sest
table(Sest)

# if partition requires to be saved as csv:
#setwd(".../your/path/MyCSV/Datasets")
# write.csv(Sest,'partition_sing_DP_10000.csv')

nclusLG=length(unique(Sest))
nclusLG

bins = as.numeric(names(table(Sest)))
freqs = as.vector(table(Sest))

label=rep(0,J)

cluster_list <- list()

for(i in 1:nclusLG)
{
  in_gr=which(Sest==bins[i])
  label[in_gr]=i
  cluster_list[[i]]=in_gr
}
label
cluster_list

# Save patient partition as R workspace
save(cluster_list, file='Cluster_10000.rdata')

#################################
# Patient data for each cluster #
#################################

# Extract covariates and information for each cluster

rm(list=ls())

ALS_data = read.csv("Combined_Lab_Well_Age_Sex_Onset_v2.csv", header = T, sep=",")
load('Cluster_10000.Rdata')

ALS_data$index = ALS_data$index+1

nrPat = nrow(ALS_data)
nrCovars = ncol(ALS_data)-3 # index, death_bin and death_date not covariates

als_tmp = as.matrix(ALS_data[, ])
als_tmp = apply(als_tmp,2, as.numeric)
patient_data_list = list()

for (cluster_nr in 1:length(cluster_list)) {
  pat_in_cl = length(cluster_list[[cluster_nr]])
  patient_data = matrix(0, nrow = pat_in_cl, ncol = nrCovars+3) # covariates + ID and death_date
  index = 1
  for (pat in 1:pat_in_cl) {
    pat_index = cluster_list[[cluster_nr]][pat]
    patient_data[index, ] <- als_tmp[pat_index, ] # where the ID and covariates is located in the als_data
    index <- index+1
  }
  patient_data_list[[cluster_nr]] <- patient_data
}

# Save patient data according to partition in R workspace list
save(patient_data_list, file='patient_cluster_10000.rdata')
load('patient_cluster_10000.rdata')

