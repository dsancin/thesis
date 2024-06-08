library(copula)
library(gridExtra)
library(cowplot)
library(gridGraphics)
#?par per il problema dei margini cambia parametri mai, mar, oma, (vedi figura 1.6 e conversazione con chatgpt)
#ora assi vanno bene, ma dimensione immagini è ancora meh
setwd('/home/davide/università/tesi magistrale/tesi magistrale/figure')

#CAPITOLO 1  
#figure 1.1 (independence copula)
d<-2
ic<-indepCopula(dim = d)
set.seed(2008)
#grid.arrange(plot1, plot2, ..., ncol=3, nrow = 3)
pdf('01_independence_copula.pdf', width = 10, height = 10) #imposta width heigth uguali
grid.arrange(wireframe2(ic, FUN=pCopula, col.4=adjustcolor("black", alpha.f=0.25)), contourplot2(ic, FUN=pCopula),ncol=2, nrow = 1)
dev.off()

#figure 1.2 (c-volume of ind cop)
n<-1000
U<-rCopula(n, copula=ic)
pdf('02_c_volume_ind.pdf', width = 10, height = 10) #imposta width heigth uguali
par(mfrow=c(1,1),  mar=c(5,4,4,2)+0.6)
plot(U, xlab=quote(U[1]), ylab=quote(U[2]),
     cex.lab=1.5, cex.axis=1.5)
segments(1/4, 1/2, 1/4, 1, lwd = 3, col=2)
segments(1/3, 1/2, 1/3, 1, lwd = 3, col=2)
segments(1/4, 1/2, 1/3, 1/2, lwd = 3, col=2)
segments(1/4, 1, 1/3, 1, lwd = 3, col=2)
dev.off()

#figure 1.3 (W M copula)
d<-2
set.seed(2008)
U<-runif(100)
pdf('03_W_M.pdf', width = 10, height = 5) #imposta width heigth uguali
par(mfrow=c(1,2), mar=c(5,4,4,2)+0.6, pty="s")  #pty=s to make the plots square shaped
plot(cbind(U,1-U), xlab=quote(U[1]), ylab=quote(U[2]))
plot(cbind(U,U), xlab=quote(U[1]), ylab=quote(U[2]))
dev.off()

#figure 1.4 (Survival copula)
#cop<-claytonCopula(2)
#set.seed(332)
#U<-rCopula(1000, copula=cop)
#V<-1-U # sample from the survival Clayton copula
#par(mfrow=c(1,2))
#plot(U, xlab = quote(U[1]),  ylab = quote(U[2]))
#plot(V, xlab = quote(V[1]),  ylab = quote(V[2]))
n <- 1000 # sample size
d <- 2 # dimension
th <- 1 # Pareto parameter
set.seed(271) # set seed (for reproducibility)>
R <- runif(n)^(-1/th) # sample radial part (here: Pareto on [1,Inf)) 
## Sample from a so-called Pareto-simplex copula
E <- matrix(rexp(n * d), nrow = n, ncol = d) # unit exponentials
S <- E / matrix(rowSums(E), nrow = n, ncol = d) # S uniformly on unit simplex
incBeta <- function(x, a, b) pbeta(x, a, b) * beta(a, b) # incomplete beta
psi <- function(t, th) t^(-1/th) * incBeta(pmin(1,t), a = 1/th, b = d) / th
U <- psi(R * S, th = th)
V<-1-U
pdf('04_survival.pdf', width = 10, height = 5) #imposta width heigth uguali
par(mfrow=c(1,2), pty='s')
plot(U, xlab = quote(U[1]), ylab = quote(U[2])) # Pareto-simplex sample
plot(V, xlab = quote(V[1]),  ylab = quote(V[2]))
dev.off() 

#figure 1.5 (radial and exchangeable copulas)
#PER FARE GRID NON FUNZIONA SU PLOT NORMALI
contourplot2(tCopula(0.5, df=2.5),
             FUN=dCopula, n.grid=64, lwd=1/2)
p2 <- grid.grab()
contourplot2(gumbelCopula(3),
             FUN=dCopula, n.grid=64, lwd=1/4, 
             pretty=FALSE, cuts=42, 
             col.regions=gray(seq(0.5,1, length.out=128)))
p4 <- grid.grab()
pdf('05_symmetry.pdf', width = 10, height = 5) #imposta width heigth uguali
plot_grid(p2,p4)
dev.off()

pdf('05_symmetry_obs.pdf', width = 10, height = 5) #imposta width heigth uguali
par(mfrow=c(1,2), pty='s')
plot(rCopula(1000, copula=tCopula(0.5, df=2.5)),
     xlab = quote(U[1]),  ylab = quote(U[2]))

plot(rCopula(1000, copula=gumbelCopula(3)),
     xlab = quote(U[1]),  ylab = quote(U[2]))
dev.off()


#figure 1.6 (normal copulas)
set.seed(999)
n<-1000
nc<- normalCopula(iTau(normalCopula(), tau = 0.5))
U<-rCopula(n, copula=setTheta(nc, value=0.1))
U2<-rCopula(n, copula = setTheta(nc, value=0.5))
U3<-rCopula(n, copula=setTheta(nc, value=0.95))
pdf('06_normal_copulas.pdf', width = 10, height = 5) #imposta width heigth uguali
par(mfrow=c(1,3), mar=c(5,4,4,2)+0.6, pty='s')
plot(U, xlab=quote(U[1]), ylab=quote(U[2]),
     cex.lab=1.5, cex.axis=1.5)
plot(U2, xlab=quote(U[1]), ylab=quote(U[2]),
     cex.lab=1.5, cex.axis=1.5)
plot(U3, xlab=quote(U[1]), ylab=quote(U[2]),
     cex.lab=1.5, cex.axis=1.5)
dev.off()


#figure 1.7 (t copulas)
pdf('07_t_copulas_fixed_df.pdf', width = 10, height = 5) #imposta width heigth uguali
par(mfrow=c(1,3), mar=c(5,4,4,2)+0.6, pty='s')
plot(rCopula(1000, copula=tCopula(-0.99, df=2.5)),
     xlab = quote(U[1]),  ylab = quote(U[2]),
     cex.lab=1.5, cex.axis=1.5)

plot(rCopula(1000, copula=tCopula(0.5, df=2.5)),
     xlab = quote(U[1]),  ylab = quote(U[2]),
     cex.lab=1.5, cex.axis=1.5)

plot(rCopula(1000, copula=tCopula(0.99, df=2.5)),
     xlab = quote(U[1]),  ylab = quote(U[2]),
     cex.lab=1.5, cex.axis=1.5)
dev.off()

pdf('07_t_copulas_fixed_P.pdf', width = 10, height = 5) #imposta width heigth uguali
par(mfrow=c(1,3), mar=c(5,4,4,2)+0.6, pty='s')
plot(rCopula(1000, copula=tCopula(0.5, df=1)),
     xlab = quote(U[1]),  ylab = quote(U[2]),
     cex.lab=1.5, cex.axis=1.5)

plot(rCopula(1000, copula=tCopula(0.5, df=3)),
     xlab = quote(U[1]),  ylab = quote(U[2]),
     cex.lab=1.5, cex.axis=1.5)

plot(rCopula(1000, copula=tCopula(0.5, df=5)),
     xlab = quote(U[1]),  ylab = quote(U[2]),
     cex.lab=1.5, cex.axis=1.5)
dev.off()

#figure 1.8 (archimedean copulas)
n <- 1000 # sample size
set.seed(27)
# Define the copula parameters chosen such 
# that Pearson's correlation coefficient,
# when computed with N(0,1) margins, is 0.7 
# tc <- frankCopula(iRho(frankCopula(), rho = 0.7))
th.g <- 3 # Gumbel copula parameter
th.c <- 3 # Clayton copula parameter
th.f<- 9 # Frank copula parameter

## Define the copulas
gc <- gumbelCopula(th.g) # Gumbel copula
cc <- claytonCopula(th.c) # Clayton copula
nf <- frankCopula(th.f) # Frank copula

## Generate copula data
U.nf <- rCopula(n, copula = nf)
U.gc <- rCopula(n, copula = gc)
U.cc <- rCopula(n, copula = cc)

pdf('08_archimedean_copulas.pdf', width = 10, height = 5) #imposta width heigth uguali
par(mfrow=c(1,3), mar=c(5,4,4,2)+0.6, pty='s')
plot(U.gc, xlab = expression(U[1]), ylab = expression(U[2]), # Gumbel copula
     cex = 0.4,
     cex.lab=1.5, cex.axis=1.5)
plot(U.cc, xlab = expression(U[1]), ylab = expression(U[2]), # Clayton copula
     cex = 0.4,
     cex.lab=1.5, cex.axis=1.5)
plot(U.nf, xlab = expression(U[1]), ylab = expression(U[2]), # Frank copula
     cex = 0.4,
     cex.lab=1.5, cex.axis=1.5)
dev.off()


#===========================================================================
#CAPITOLO 3
library(xts)
library(copula)
library(zoo)
library(tseries) 
load("/home/davide/università/tesi magistrale/dati/arpa/only_rain.RData")
clean_rain<-as.matrix(clean_rain_xts)
colnames(clean_rain)<-colnames(clean_rain_xts)
clean_rain<-pobs(clean_rain)


#UNIVARITE TIME SERIES
pdf('/home/davide/università/tesi magistrale/tesi magistrale/figure/cap3/rain_ts.pdf', width = 10, height = 10) #imposta width heigth uguali
trace(plot.zoo, 
      quote(mtext <- function(...) graphics::mtext(..., cex = 0.5, las = 1)))
plot.zoo(clean_rain_xts, oma = c(6, 5, 5, 0), main=" ")
untrace(plot.zoo)
dev.off()

#ACF PLOTS
pdf('/home/davide/università/tesi magistrale/tesi magistrale/figure/cap3/rain_acf1.pdf', width = 10, height = 10) #imposta width heigth uguali
par(mfrow=c(3,3))  
for (x in colnames(clean_rain_xts[,1:9])){
  acf(ts(clean_rain_xts[,x]), main=x)    #plot of the ACF
}
dev.off()
pdf('/home/davide/università/tesi magistrale/tesi magistrale/figure/cap3/rain_acf2.pdf', width = 10, height = 10) #imposta width heigth uguali
par(mfrow=c(3,3))  
for (x in colnames(clean_rain_xts[,10:18])){
  acf(ts(clean_rain_xts[,x]), main=x)    #plot of the ACF
}
dev.off()

#PAIRWISE SCATTERPLOT
pdf('/home/davide/università/tesi magistrale/tesi magistrale/figure/cap3/pairwise_rain.pdf', width = 10, height = 10)
pairs2(clean_rain, cex=0.4)
dev.off()

#HEATMAPS
library(heatmaply)
#pdf('/home/davide/università/tesi magistrale/tesi magistrale/figure/cap3/tau_heat.pdf', width = 10, height = 10)
heatmaply_cor(corKendall(clean_rain), xlab = "Stations", 
              ylab = "Stations", dendrogram = "none", scale="none")
#dev.off()
heatmaply_cor(lower_tdc, xlab = "Stations", 
              ylab = "Stations", dendrogram = "none", scale="none")

heatmaply_cor(upper_tdc, xlab = "Stations", 
              ylab = "Stations", dendrogram = "none", scale="none")
dev.off()
#K-PLOTS
library(VineCopula)
pdf('/home/davide/università/tesi magistrale/tesi magistrale/figure/cap3/Kplot.pdf', width = 10, height = 4)
par(mfrow=c(1,3), pty='s')
BiCopKPlot(clean_rain[,"Codroipo"] ,clean_rain[,"Tarvisio"])
BiCopKPlot(clean_rain[,"Codroipo"] ,clean_rain[,"Musi"])
BiCopKPlot(clean_rain[,"Codroipo"] ,clean_rain[,"Palazzolo.dello.Stella"])
dev.off()

#eventually plot it only for the clusters, but due to reason above does not work anyway



#CLUSTERS
dissimilarity_diff_upper<- ((1-upper_tdc))^0.5
upper_clustering_diff_average<-hclust(d=as.dist(dissimilarity_diff_upper), method="average")
upper_clustering_diff_complete <- hclust(d=as.dist(dissimilarity_diff_upper), method = "complete")

pdf('/home/davide/università/tesi magistrale/tesi magistrale/figure/cap3/den.pdf', width = 10, height = 10) #imposta width heigth uguali
par(mfrow=c(1,2))  
plot(upper_clustering_diff_average, main="Average", xlab="", sub="")
abline(h=0.91, col="blue")
plot(upper_clustering_diff_complete, main="Complete", xlab="", sub="")
abline(h=0.91, col="blue")
dev.off()




#COMPLETE CLUSTER DIFF
cluster11<-clean_rain[,c("Enemonzo", "Tarvisio")]
cluster12<-clean_rain[,c("Brugnera", "Vivaro")]
cluster13<-clean_rain[,c("Fagagna", "Gemona.del.Friuli", "Musi")]  
cluster14<-clean_rain[,c("Talmassons", "Cividale.del.Friuli", "Udine.S.O.", "Palazzolo.dello.Stella","San.Vito.al.Tgl.", "Codroipo")]
cluster15<-clean_rain[,c("Fossalon", "Sgonico.Zgonik", "Gradisca.d.Is.", "Capriva.del.Friuli", "Cervignano.del.Friuli")]

#CLUSTER VALIDATION
library(heatmaply)
library(VineCopula)

lambda_clust14<-fitLambda(cluster14, lower.tail  = FALSE)
tau_clust14<-corKendall(cluster14)
colnames(lambda_clust14)=colnames(cluster14)
rownames(lambda_clust14)=colnames(cluster14)
colnames(tau_clust14)=colnames(cluster14)
rownames(lambda_clust14)=colnames(cluster14)


#heatmaply_cor(corKendall(cluster14), xlab = "Stations", 
#              ylab = "Stations", dendrogram = "none", scale="none", main="")
heatmaply_cor(tau_clust14, xlab = "Stations", 
              ylab = "Stations", dendrogram = "none", scale="none", main="")
dev.off()

pairwise14<-combn(colnames(cluster14), 2) #this create all the combinations of the stations (no repetitions)
pdf('/home/davide/università/tesi magistrale/tesi magistrale/figure/cap3/Kplots_cluster14_1.pdf', width = 10, height = 10) #imposta width heigth uguali
par(mfrow=c(3,3), pty='s')  #soluzione temporanea con più plot
for (i in 1:(length(pairwise14)/4)){ 
  BiCopKPlot(clean_rain[,pairwise14[1,i]],clean_rain[,pairwise14[2,i]], main = pairwise14[1:2,i])
}
dev.off()

pdf('/home/davide/università/tesi magistrale/tesi magistrale/figure/cap3/Kplots_cluster14_2.pdf', width = 10, height = 10) #imposta width heigth uguali
par(mfrow=c(3,3), pty='s')  #soluzione temporanea con più plot
for (i in (length(pairwise14)/4):(length(pairwise14)/2)){ 
  BiCopKPlot(clean_rain[,pairwise14[1,i]],clean_rain[,pairwise14[2,i]], main = pairwise14[1:2,i])
}
dev.off()

pdf('/home/davide/università/tesi magistrale/tesi magistrale/figure/cap3/scatter_cluster14.pdf', width = 10, height = 10) #imposta width heigth uguali
pairs2(cluster14, main="")
dev.off()






