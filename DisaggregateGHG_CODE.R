# DIS-AGGREGATED GHG CODE #
#############################################################################
req_packages<-c("devtools","glmnet","igraph","readxl","grid","gridExtra","Rmisc","ggplot2","sandwich","bootUR","mFilter","psych","parallel","scales","MASS")
install.packages(req_packages)
lapply(req_packages, require, character.only=TRUE)
devtools::install_github("Marga8/HDGCvar")
library(HDGCvar)

#load dataset and plot it
# AggregateGHG <- readRDS("own_path/AggregateGHG.rds")
DATA_R_clima3<-DisaggregateGHG[,-1] #-1 to remove the dates column

Dates<-c(1871:2014)
dataset_x_ggplot<-cbind(Dates, DATA_R_clima3)

#Plot GHG
plot_GHG<-list()
plot_GHG[[1]]<-ggplot()+ geom_line(data = dataset_x_ggplot, aes(x = Dates, y = CO2, color = "#00AFBB")) +
  geom_line(data = dataset_x_ggplot, aes(x = Dates, y = CH4, color = "coral1")) +
  geom_line(data = dataset_x_ggplot, aes(x = Dates, y = N2O, color = "deeppink4"))+
  xlab('Years') +
  ylab('W/m^2')+ scale_color_discrete(name = "GHG", labels = c("CO2", "CH4", "N2O"))
grid.arrange(grobs=plot_GHG, ncol=1)


#Plot all other series (except aggregate GHG)
plots_series<-list()

plots_series[[1]]<-ggplot(data = dataset_x_ggplot, aes(x = Dates, y = T))+ geom_line(color = "#00AFBB", size = 0.5) +
  xlab("Years")  + ylab("Temperature Anomalies")
plots_series[[2]]<-ggplot(data = dataset_x_ggplot, aes(x = Dates, y = S))+ geom_line(color = "#00AFBB", size = 0.5) +
  xlab("Years")  + ylab("Solar Activity")
plots_series[[3]]<-ggplot(data = dataset_x_ggplot, aes(x = Dates, y = V))+ geom_line(color = "#00AFBB", size = 0.5) +
  xlab("Years")  + ylab("Stratospheric Aerosols")
plots_series[[4]]<-ggplot(data = dataset_x_ggplot, aes(x = Dates, y = Y))+ geom_line(color = "#00AFBB", size = 0.5) +
  xlab("Years")  + ylab("GDP (log 2010 USD)")
plots_series[[5]]<-ggplot(data = dataset_x_ggplot, aes(x = Dates, y = A))+ geom_line(color = "#00AFBB", size = 0.5) +
  xlab("Years")  + ylab("Tropospheric Aerosols")
plots_series[[6]]<-ggplot(data = dataset_x_ggplot, aes(x = Dates, y = N))+ geom_line(color = "#00AFBB", size = 0.5) +
  xlab("Years")  + ylab("El Nino SOI")
plots_series[[7]]<-ggplot(data = dataset_x_ggplot, aes(x = Dates, y = O))+ geom_line(color = "#00AFBB", size = 0.5) +
  xlab("Years")  + ylab("Ocean Heat Content")
grid.arrange(grobs=plots_series, ncol=2)

#lag-length selection (VAR order)
selected_lag<-HDGCvar::lags_upbound_BIC(DATA_R_clima3,p_max=10)
print(selected_lag)

#compute network
network<-HDGCvar::HDGC_VAR_all(DATA_R_clima3, p = selected_lag, d = 2, bound = 0.5 * nrow(DATA_R_clima3),
                               parallel = TRUE, n_cores = NULL)

#plot network
HDGCvar::Plot_GC_all(network, Stat_type="FS_cor",alpha=0.10, multip_corr=list(T,"none"),
                     directed=T, layout=layout.circle, main="Climate Network",edge.arrow.size=.2,
                     vertex.size=16, vertex.color=c("lightblue"), vertex.frame.color="blue",
                     vertex.label.size=4,vertex.label.color="black",vertex.label.cex=0.8,
                     vertex.label.dist=c(0,0,0,0,0,0,0,0,0,0), edge.curved=0,cluster=list(F,5,"black",0.8,1,0))

## HEAT MAP table for p-values ##
##########################################################################
pvals<-as.matrix(network[["tests"]][,,2,2])
pv_max <- 0.15
trunc_pvals <- pmin(pvals, pv_max)
series <-colnames(DATA_R_clima3)
GC_df <- data.frame(pvalue = c(trunc_pvals),
                    series_to = factor(rep(series, ncol(DATA_R_clima3)), levels = series),
                    series_from = factor(rep(series, each = ncol(DATA_R_clima3)), levels = series))

ggplot(GC_df, aes(x = series_from, y = series_to, fill = pvalue)) +
  geom_tile(colour = "grey50") +
  scale_y_discrete(limits = rev(levels(GC_df$series_from))) +
  labs(x = "Granger causality from", y = "Granger causality to", fill = "p-value") +
  scale_fill_gradient2(low = "darkblue", mid = "blue", high = "white", midpoint = pv_max / 2) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 0, size = 9),
        axis.text.y = element_text(size = 9),
        legend.position = "none")


#highlight paths
##########################################################################
input<-as.matrix(network[["tests"]][,,2,2])#p_values
input<-t(input)
input[input < 0.1] <- 1 #put =1 values < alpha
input[is.na(input)] <- 0 #put =0 the diagonal
input[input != 1] <- 0 #put =0 values > alpha
my_graph=graph_from_adjacency_matrix(input, mode='directed',diag=F,add.rownames = TRUE )
V(my_graph)$label = rownames(input)
# Preferences: use V(my_graph) to see vertices numbers
FROM_V <- 8 #put =9 for CH4, =10 for N2O
TO_V <- 5 # temperature

# Calculate all simple paths from FROM_V to TO_V as list of vertecy sequences
graph_outlet <- all_simple_paths(my_graph,from=FROM_V,to=TO_V)
graph_outlet

V(my_graph)$color <- "white"
V(my_graph)$color[unique(unlist(graph_outlet))] <- "lightblue"
V(my_graph)$color[c(FROM_V,TO_V)] <- "yellow"

# Colour each of the paths
E(my_graph)$color <- "gray"
E(my_graph)$width<- 1

lapply(graph_outlet, function(x) E(my_graph, path=x)$color <<- "red")
lapply(graph_outlet, function(x) E(my_graph, path=x)$width <<- 3)

# Plot all paths and mark the group of vertecies through which paths flow
plot(my_graph,layout=layout.circle,main="Climate Network",edge.arrow.size=.5,
     vertex.size=c(16,16,16,16,25,16,16,25,16,16),vertex.label.size=4)



################ BLOCK Granger Causality ################
#########################################################
#define lag-length
selected_lag=3 #=15

HDGC_VAR(GCpair=list("GCto"="S", "GCfrom"=c("CO2","N2O","CH4")), DATA_R_clima3, p=selected_lag, d=2)$test
HDGC_VAR(GCpair=list("GCto"="V", "GCfrom"=c("CO2","N2O","CH4")), DATA_R_clima3, p=selected_lag, d=2)$test
HDGC_VAR(GCpair=list("GCto"="Y", "GCfrom"=c("CO2","N2O","CH4")), DATA_R_clima3, p=selected_lag, d=2)$test
HDGC_VAR(GCpair=list("GCto"="A", "GCfrom"=c("CO2","N2O","CH4")), DATA_R_clima3, p=selected_lag, d=2)$test
HDGC_VAR(GCpair=list("GCto"="T", "GCfrom"=c("CO2","N2O","CH4")), DATA_R_clima3, p=selected_lag, d=2)$test
HDGC_VAR(GCpair=list("GCto"="N", "GCfrom"=c("CO2","N2O","CH4")), DATA_R_clima3, p=selected_lag, d=2)$test
HDGC_VAR(GCpair=list("GCto"="O", "GCfrom"=c("CO2","N2O","CH4")), DATA_R_clima3, p=selected_lag, d=2)$test



Adjamat<-matrix(0,ncol=8,nrow=8)
rownames(Adjamat)<-c("S","V","Y","A","T","N","O","GHGs")
colnames(Adjamat)<-c("S","V","Y","A","T","N","O","GHGs")
Adjamat[8,c(2,4,5,6)]<-1
Grafico<-graph_from_adjacency_matrix(Adjamat, mode='directed',diag=F,add.rownames = TRUE )
V(Grafico)$label = rownames(Adjamat)
plot(Grafico,directed=T, main="Climate Network",edge.arrow.size=.5,
     vertex.size=16, vertex.color=c("lightblue"), vertex.frame.color="blue",
     vertex.label.size=4,vertex.label.color="black",vertex.label.cex=0.8,
     vertex.label.dist=c(0,0,0,0,0,0,0,0), edge.curved=0)
