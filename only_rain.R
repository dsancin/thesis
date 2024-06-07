library(xts)
library(copula)
library(zoo)
library(tseries) 
library(forecast)

#use stations with 20 years of data
#setwd("/home/davide/università/tesi magistrale/dati/arpa/extreme value all_20") #set the working directory
load("/home/davide/università/tesi magistrale/dati/arpa/only_rain.RData")

#rain<- read.csv("max_rain_extreme_value.csv", sep=",", header=TRUE) #import data from the csv

#rain[,"date"] <- paste(rain[,"date"], "01", sep = "-") #add dummy string for day (without this it does not work)
#rain[,"date"]<-as.Date(rain[,"date"], format="%Y-%m-%d")  #convert to date type
#rain_xts<-xts(rain[, -1], order.by = as.POSIXct(rain$date)) #convert to xts

#clean_rain_xts<-rain_xts[, colSums(is.na(rain_xts)) == 0]

#create dataframe of the coordinates of the station
coordinates<- data.frame(
  station = colnames(clean_rain_xts),
  alt = c(22, 85, 8, 127, 37, 438, 148, 0, 184, 29, 600, 5, 21, 268, 16, 794, 91, 142),
  lat=c(45.91792, 45.958094, 45.849486, 46.080442, 45.952356, 46.410416, 46.101692, 45.714768, 46.261297, 45.889791, 46.312663, 45.805720, 45.895661, 45.738004, 45.882311, 46.510775, 46.035212, 46.076529),
  long=c(12.545003, 13.512333, 13.337015, 13.420014, 13.002742, 12.862536, 13.073886, 13.458865, 13.122088, 13.481807, 13.274682, 13.052598, 12.814989, 13.742056, 13.15779, 13.551886, 13.226672, 12.768814)
)

#plot of the time series
trace(plot.zoo, 
      quote(mtext <- function(...) graphics::mtext(..., cex = 0.5, las = 1)))
plot.zoo(clean_rain_xts, oma = c(6, 5, 5, 0))
untrace(plot.zoo)


#SERIAL INDEPENDENCE TEST
par(mfrow=c(2,3))  #soluzione temporanea con più plot
for (x in colnames(clean_rain_xts)){
  acf(ts(clean_rain_xts[,x]), main=x)    #plot of the ACF
}

for (x in colnames(clean_rain_xts)){#round on p-value not working
  print(paste0("p-value of Ljung-Box test on rain in ", x , " is:", Box.test(clean_rain_xts[,x], lag = 20, type = "Ljung-Box")$p.value))
  print(paste0("p-value of Ljung-Box test on squared rain in ", x , " is:", Box.test(clean_rain_xts[, x]^2, lag = 20, type = "Ljung-Box")$p.value))
}


clean_rain<-as.matrix(clean_rain_xts)
colnames(clean_rain)<-colnames(clean_rain_xts)


colnames(clean_rain_xts)
#seasonally adjust time series
auto.arima(clean_rain_xts[,6])#enemonzo
auto.arima(clean_rain_xts[,9])#gemona 
auto.arima(clean_rain_xts[,11])#gemona 
auto.arima(clean_rain_xts[,16])#tarvisio

acf(residuals(arima(clean_rain_xts[,6],c(1,0,0))))
Box.test(residuals(arima(clean_rain_xts[,6],c(1,0,0))))

acf(residuals(arima(clean_rain_xts[,9],c(1,0,0))))
Box.test(residuals(arima(clean_rain_xts[,9],c(1,0,0))))

acf(residuals(arima(clean_rain_xts[,11],c(2,0,2))))
Box.test(residuals(arima(clean_rain_xts[,11],c(2,0,2))))

acf(residuals(arima(clean_rain_xts[,16],c(1,0,1))))
Box.test(residuals(arima(clean_rain_xts[,16],c(1,0,1))))

clean_rain<-pobs(clean_rain)
clean_rain[,6]<-pobs(residuals(arima(clean_rain_xts[,6],c(1,0,0))))
clean_rain[,9]<-pobs(residuals(arima(clean_rain_xts[,9],c(1,0,0))))
clean_rain[,11]<-pobs(residuals(arima(clean_rain_xts[,11],c(2,0,2))))
clean_rain[,16]<-pobs(residuals(arima(clean_rain_xts[,16],c(1,0,1))))


#for (i in colnames(clean_rain)){
#hist(clean_rain[,i])
# }



#pairwise scatterplots of stations for a given variable
#pdf('/home/davide/università/tesi magistrale/tesi magistrale/figure/cap3/pairwise_rain.pdf', width = 10, height = 10)
pairs2(clean_rain, cex=0.4)

#kendall_tau heatmap
library(heatmaply)
heatmaply_cor(corKendall(clean_rain), xlab = "Stations", 
              ylab = "Stations", dendrogram = "none", scale="none")

#STATISTICAL TESTS
#test of independency (basta guardare a hetmap del kendall tau, tutti maggiori di 0)
#test of extreme value dependence
pairwise<-combn(colnames(clean_rain), 2) #this create all the combinations of the stations (no repetitions)
#this allow to works only on the lower triangular matrix (all elements below the diagonal)

#pairwise<-t(expand.grid(colnames(clean_rain),colnames(clean_rain)) #this create all the combinations of the stations with repetitions
#clean_rain[,pairwise[1:2,1]]#this to call the first couple
#clean_rain[,pairwise[1:2,2]]#this to call the second couple

#this is the pairwise test (slow)
#test of exchangeability
exchangeability_test<-list()
pb = txtProgressBar(min = 0, max = length(pairwise)/2, initial = 0, style = 3) #initialize progress bar
for (i in 1:(length(pairwise)/2)){
  set.seed(42)
  setTxtProgressBar(pb,i)
  exchangeability_test[i]<-exchTest(clean_rain[,pairwise[1:2,i]])$p.value
}
close(pb)

exchangeability_test_matrix<-rbind(pairwise, exchangeability_test)
exchangeability_test_matrix <- matrix(unlist(exchangeability_test_matrix), nrow = 3, byrow = FALSE)
which(as.numeric(exchangeability_test_matrix[3,]) < 0.05)

#test of radial symmetry
radial_test<-list()
pb = txtProgressBar(min = 0, max = length(pairwise)/2, initial = 0, style = 3) #initialize progress bar
for (i in 1:(length(pairwise)/2)){
  set.seed(42)
  setTxtProgressBar(pb,i)
  radial_test[i]<-radSymTest(clean_rain[,pairwise[1:2,i]])$p.value
}
close(pb)

radial_test_matrix<-rbind(pairwise, radial_test)
radial_test_matrix <- matrix(unlist(radial_test_matrix), nrow = 3, byrow = FALSE)
which(as.numeric(radial_test_matrix[3,]) < 0.05)


#test of extreme value dependence
extreme_value_test<-list()
pb = txtProgressBar(min = 0, max = length(pairwise)/2, initial = 0, style = 3) #initialize progress bar
for (i in 1:(length(pairwise)/2)){ #took more or less 10 minutes while running other stuff, only 4 otherwise
  set.seed(42)
  extreme_value_test[i]<-evTestK(clean_rain[,pairwise[1:2,i]])$p.value
  setTxtProgressBar(pb,i)
}
close(pb)
extreme_value_test_matrix<-rbind(pairwise,as.numeric(extreme_value_test)) #all low value
extreme_value_test_matrix<-matrix(unlist(extreme_value_test_matrix), nrow = 3, byrow = FALSE)

which(as.numeric(extreme_value_test_matrix[3,]) > 0.05)


#KENDALL DISTRIBUTION FUNCTION (VIA K-plot)
library(VineCopula)

par(mfrow=c(4,4))  #soluzione temporanea con più plot
for (i in 1:(length(pairwise)/2)){ 
  BiCopKPlot(clean_rain[,pairwise[1,i]],clean_rain[,pairwise[2,i]], main = pairwise[1:2,i])
}

#TAIL DEPENDENCE COEFFICIENT ESTIMATION

lower_tdc<-fitLambda(clean_rain, lower.tail  = TRUE)
upper_tdc<-fitLambda(clean_rain, lower.tail  = FALSE)

colnames(lower_tdc)=colnames(clean_rain)
rownames(lower_tdc)=colnames(clean_rain)

colnames(upper_tdc)=colnames(lower_tdc)
rownames(upper_tdc)=colnames(lower_tdc)


heatmaply_cor(lower_tdc, xlab = "Stations", 
          ylab = "Stations", dendrogram = "none", scale="none")

heatmaply_cor(upper_tdc, xlab = "Stations", 
          ylab = "Stations", dendrogram = "none", scale="none")

#DISSIMILARITY MATRIX AND CLUSTERING
dissimilarity_diff_upper<- (1-upper_tdc)^0.5

upper_clustering_diff_average<-hclust(d=as.dist(dissimilarity_diff_upper), method="average")
upper_clustering_diff_complete <- hclust(d=as.dist(dissimilarity_diff_upper), method = "complete")

plot(upper_clustering_diff_average, main="Average")
plot(upper_clustering_diff_complete, main="Complete")
abline(h=0.91, col="blue")

#AVERAGE CLUSTER log
cluster11<-clean_rain[,c("Enemonzo", "Tarvisio")]
cluster12<-clean_rain[,c("Brugnera", "Vivaro")]
cluster13<-clean_rain[,c("Fagagna", "Gemona.del.Friuli", "Musi")]  
cluster14<-clean_rain[,c("Talmassons", "Cividale.del.Friuli", "Udine.S.O.", "Palazzolo.dello.Stella","San.Vito.al.Tgl.", "Codroipo")]
cluster15<-clean_rain[,c("Fossalon", "Sgonico.Zgonik", "Gradisca.d.Is.", "Capriva.del.Friuli", "Cervignano.del.Friuli")]

#library(png)
#panel.Kplot<-function(x,y)
#{
#  temp_file <- tempfile(fileext = ".png")
#  png(temp_file, width = 240, height = 240)
#  BiCopKPlot(x,y)
#  dev.off()
#  img <- readPNG(temp_file)
#  rasterImage(img, xleft = par("usr")[1], xright = par("usr")[2], 
#              ybottom = par("usr")[3], ytop = par("usr")[4])  
#}
#pairs(cluster1, lower.panel = panel.Kplot)

#pairs2(cluster1, main="cluster1")
#pairs2(cluster2, main="cluster2")
#pairs2(cluster3, main="cluster3")
#pairs2(cluster4, main="cluster4")

#heatmaply_cor(corKendall(cluster1), xlab = "Stations", 
#              ylab = "Stations", dendrogram = "none", scale="none", main="cluster1")
#heatmaply_cor(corKendall(cluster2), xlab = "Stations", 
#              ylab = "Stations", dendrogram = "none", scale="none", main="cluster2")
#heatmaply_cor(corKendall(cluster3), xlab = "Stations", 
#              ylab = "Stations", dendrogram = "none", scale="none", main="cluster3")
#heatmaply_cor(corKendall(cluster4), xlab = "Stations", 
#              ylab = "Stations", dendrogram = "none", scale="none", main="cluster4")

dev.off()
pairs2(cluster11, main="cluster11")
pairs2(cluster12, main="cluster12")
pairs2(cluster13, main="cluster13")
pairs2(cluster14, main="cluster14")
pairs2(cluster15, main="cluster15")


heatmaply_cor(fitLambda(cluster11, lower.tail  = FALSE), xlab = "Stations", 
              ylab = "Stations", dendrogram = "none", scale="none", main="cluster11")
heatmaply_cor(fitLambda(cluster12, lower.tail  = FALSE), xlab = "Stations", 
              ylab = "Stations", dendrogram = "none", scale="none", main="cluster12")
heatmaply_cor(fitLambda(cluster13, lower.tail  = FALSE), xlab = "Stations", 
              ylab = "Stations", dendrogram = "none", scale="none", main="cluster13")
heatmaply_cor(fitLambda(cluster14, lower.tail  = FALSE), xlab = "Stations", 
              ylab = "Stations", dendrogram = "none", scale="none", main="cluster14")
heatmaply_cor(fitLambda(cluster15, lower.tail  = FALSE), xlab = "Stations", 
              ylab = "Stations", dendrogram = "none", scale="none", main="cluster14")


heatmaply_cor(corKendall(cluster11), xlab = "Stations", 
              ylab = "Stations", dendrogram = "none", scale="none", main="cluster11")
heatmaply_cor(corKendall(cluster12), xlab = "Stations", 
              ylab = "Stations", dendrogram = "none", scale="none", main="cluster12")
heatmaply_cor(corKendall(cluster13), xlab = "Stations", 
              ylab = "Stations", dendrogram = "none", scale="none", main="cluster13")
heatmaply_cor(corKendall(cluster14), xlab = "Stations", 
              ylab = "Stations", dendrogram = "none", scale="none", main="cluster14")
heatmaply_cor(corKendall(cluster15), xlab = "Stations", 
              ylab = "Stations", dendrogram = "none", scale="none", main="cluster14")



#================================
#FIT DA AGGIORNARE....

#COPULA FITTING ON CLUSTERS
#AVERAGE CLUSTERS LOOK BETTER
#CLUSTER 11
N<-1000
radSymTest(cluster11)$p.value 
exchTest(cluster11)$p.value

t.copula_cluster11 <- tCopula(dim = dim(cluster11)[2], dispstr = "un", df.fixed = TRUE)
gofCopula(t.copula_cluster11, x=cluster11, N=N, simulation="mult")  

#CLUSTER 12
radSymTest(cluster12)$p.value #low p-value
exchTest(cluster12)$p.value 

fit_cluster12 <- fitCopula(gumbelCopula(dim = dim(cluster12)[2]),
                           data = cluster12, method = "mpl")  
gofCopula(fit_cluster12@copula, x = cluster12) 


#CLUSTER 13
radSymTest(cluster13)$p.value #lowish
exchTest(cluster13)$p.value
#both exchangeable and radially symmetric
#FIT and GOF A t-COPULA
t.copula_cluster13 <- tCopula(dim = dim(cluster13)[2], dispstr = "un", df.fixed = TRUE)
gofCopula(t.copula_cluster13, x=cluster13, N=N, simulation="mult")  
#very high p-value


#CLUSTER 14
radSymTest(cluster14)$p.value #low p-value
exchTest(cluster14)$p.value #lowish p-value
fit_cluster14 <- fitCopula(claytonCopula(dim = dim(cluster12)[2]),
                                data = cluster14, method = "mpl")  
gofCopula(fit_cluster14@copula, x=cluster14)  


#CLUSTER 15


#save(coordinates, clean_rain_xts, extreme_value_test_matrix, exchangeability_test_matrix, pairwise,
#     lower_tdc, upper_tdc, radial_test_matrix, clean_rain,
#     file = "/home/davide/università/tesi magistrale/dati/arpa/only_rain.RData")











#==============================================================================#
#STESSA ROBA FACENDO CLUSTER SU KENDALL TAU INVECE DI LAMBDA
tau_matrix<-corKendall(clean_rain)

#dissimilarity_log_ken<- -log(tau_matrix)

dissimilarity_diff_ken<- (2*(1-tau_matrix))^0.5

#test different clustering methods
#ken_clustering_log_average<-hclust(d=as.dist(dissimilarity_log_ken), method="average")
#ken_clustering_log_complete <- hclust(d=as.dist(dissimilarity_log_ken), method = "complete")

ken_clustering_diff_average<-hclust(d=as.dist(dissimilarity_diff_ken), method="average")
ken_clustering_diff_complete <- hclust(d=as.dist(dissimilarity_diff_ken), method = "complete")

#dendrogram plots
#par(mfrow = c(2, 2))
#brutti cluster
#plot(ken_clustering_log_average, main="Average")
#abline(h=0.8, col="blue")

#plot(ken_clustering_log_complete, main="Complete")
#abline(h=1.0, col="blue")

#average meh, complete imilar to two above(seems geographically worse)
#brutti cluster
plot(ken_clustering_diff_average, main="Average")

plot(ken_clustering_diff_complete, main="Complete")
abline(h=1.1, col="blue")

cluster1<-clean_rain[,c("Enemonzo", "Tarvisio", "Gemona.del.Friuli", "Musi")]
cluster2<-clean_rain[,c("Talmassons", "Cividale.del.Friuli", "Udine.S.O.", "Palazzolo.dello.Stella","San.Vito.al.Tgl.", "Codroipo", "Cervignano.del.Friuli", "Brugnera", "Vivaro", "Fagagna")]
cluster3<-clean_rain[,c("Fossalon", "Sgonico.Zgonik", "Gradisca.d.Is.", "Capriva.del.Friuli")]

pairs2(cluster1)
heatmaply_cor(fitLambda(cluster1, lower.tail  = FALSE), xlab = "Stations", 
              ylab = "Stations", dendrogram = "none", scale="none", main="lambda cluster1")
heatmaply_cor(corKendall(cluster1), xlab = "Stations",  
              ylab = "Stations", dendrogram = "none", scale="none", main="kendall tau cluster1")

pairs2(cluster2)
heatmaply_cor(fitLambda(cluster2, lower.tail  = FALSE), xlab = "Stations", 
              ylab = "Stations", dendrogram = "none", scale="none", main="lambda cluster2")
heatmaply_cor(corKendall(cluster2), xlab = "Stations", 
                         ylab = "Stations", dendrogram = "none", scale="none", main="kendall tau cluster2")

pairs2(cluster3)
heatmaply_cor(fitLambda(cluster3, lower.tail  = FALSE), xlab = "Stations", 
              ylab = "Stations", dendrogram = "none", scale="none", main="lambda cluster3")
heatmaply_cor(corKendall(cluster3), xlab = "Stations", 
                         ylab = "Stations", dendrogram = "none", scale="none", main="kendall tau cluster3")




