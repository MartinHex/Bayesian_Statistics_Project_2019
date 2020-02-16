path = getwd()
DataSetPath = paste(path, "/MyCSV/DataSets")
WorkSpacePath = paste(path,"MCMC_output")

setwd(DataSetPath)
rm(list=ls())

data = read.table("Combined_Lab_Well_Age_Sex_Onset_v2.csv",header=T, sep = ',')
data = as.matrix(data)

setwd("C:/Users/marti/OneDrive/Programmering/R/School/BayesianStatistics")
partition = read.table("Partition_multipleDP.csv",header=T, sep = ',')
partition = partition['x']

J = dim(data)[1]
nbrCovariates = dim(data)[2] -3


# ------------------------- Seperates data into one dataframe for each cluster -------------------------------

bins = as.numeric(names(table(partition)))
freqs = as.vector(table(partition))
nbrClusters = max(partition)

cluster_list=list()

for(i in 1:nbrClusters)
{
  in_gr=which(partition==bins[i])
  cluster_list[[i]]=in_gr
}

patient_data_list = list()
for (cluster_nr in 1:length(cluster_list)) {
  pat_in_cl = length(cluster_list[[cluster_nr]])
  patient_data = matrix(0, nrow = pat_in_cl, ncol = nbrCovariates+3) # covariates + ID and death_date
  index = 1
  for (pat in 1:pat_in_cl) {
    pat_index = cluster_list[[cluster_nr]][pat]
    patient_data[index, ] <- data[pat_index, ] # where the ID and covariates is located in the als_data
    index <- index+1
  }
  patient_data_list[[cluster_nr]] <- patient_data
}




# ----------------------------- Create and run the model with given cluster --------------------

library(rjags)
myCluster = 4
CurrentData = patient_data_list[myCluster][[1]]
n.data = dim(CurrentData)[1]
dat <- list( Age=CurrentData[,6], t=CurrentData[,3], Sex=CurrentData[,7], Is=CurrentData[,2],
             FVC=CurrentData[,5], Cre=CurrentData[,8], Hemoglobin = CurrentData[,9],
             nbrCovariates=nbrCovariates, n.data=n.data)

inits = list(beta=rep(0,7),alpha=1)

modelAFT <- jags.model("Survival.bug",data=dat,inits=inits,n.adapt=50000,n.chains=1)



update(modelAFT,50000)


variable.names=c("mu", "alpha","beta") # monitoring#
n.iter=50000
thin=20

outputAFT <- coda.samples(model=modelAFT,variable.names=variable.names,n.iter=n.iter,thin=thin)

setwd(WorkSpacePath)
save(outputAFT,file='Survival_Cluster4_AFT_multipleDP.rdata')


data.out <- as.matrix(outputAFT)
data.out <- data.frame(data.out)
attach(data.out)
n.chain <- dim(data.out)[1]


# --------------------------------- Visualize traceplot and autocorrelation ------------------------------------

x11()
par(mfrow=c(2,3))
acf(data.out[,'beta.2.'],lwd=3,col="red3",main="beta: Limb_Onset")
acf(data.out[,'beta.3.'],lwd=3,col="red3",main="beta: SEX")
acf(data.out[,'beta.4.'],lwd=3,col="red3",main="beta: AGE")
acf(data.out[,'beta.5.'],lwd=3,col="red3",main="beta: FVC_Score")
acf(data.out[,'beta.6.'],lwd=3,col="red3",main="beta: CREATINE")
acf(data.out[,'beta.7.'],lwd=3,col="red3",main="beta: HEMOGLOBIN")


x11()
par(mfrow=c(2,3))
plot(ts(data.out[,'beta.2.']),xlab="t",ylab="beta: Limb_Onset")
plot(ts(data.out[,'beta.3.']),xlab="t",ylab="beta: SEX")
plot(ts(data.out[,'beta.4.']),xlab="t",ylab="beta: AGE")
plot(ts(data.out[,'beta.5.']),xlab="t",ylab="beta: FVC_Score")
plot(ts(data.out[,'beta.6.']),xlab="t",ylab="beta: CREATINE")
plot(ts(data.out[,'beta.7.']),xlab="t",ylab="beta: HEMOGLOBIN")

# --------------------------------- Plot survival function with Henrik's script ----------------------------------

timeSteps = 1000
Surv = matrix(0, nrow = n.chain, ncol=1)
SurvQ = matrix(0, nrow = timeSteps, ncol = 3)
pat_mu = data.out$mu.110. #<--------------------- Change patient here!!!!
pat_alpha = data.out$alpha

for (t in 1:timeSteps) {
  for (iter in 1:n.chain) {
    Surv[iter, ] =  exp(-exp(log(t)-pat_mu[iter])*pat_alpha[iter])
  }
  SurvQ[t, ]=quantile(Surv, c(0.025,0.5,0.975))
}

library(ggplot2)

Surv_dataframe = data.frame(SQ1 = SurvQ[, 1], SQ2 = SurvQ[, 2], SQ3 = SurvQ[, 3], time=1:timeSteps)
p <- ggplot(Surv_dataframe, aes(x=time, y=SQ2, ymin=SQ1, ymax=SQ3)) +
  geom_line() +
  geom_ribbon(alpha=0.2)
p + labs(x="time", y="S", title = "Survival quantiles for patient 110")
