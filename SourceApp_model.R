################################################################################
# 1) Takes input data from a file excel with the following columns: 
#    filtro,	tc,	ec,	oc,	pmc_tc,	errpmc_tc,	pmc_ec,	errpmc_ec,	pmc_oc (empty)
# 2) Assumptions: 
#     - Fossil fuel fraction of modern carbon fm_ff = 0;
#     - Biomass burning fraction of modern carbon fm_bb = 1.03 
#       calculated applying the Chapman model on the historical atmospheric 
#       pmc series (Jungfraujoch);
#     - Biogenic fraction of modern carbon fm_bio = 0.997 from the pmc 
#       Jungfraujoch time series 2020-2021;
#     - Ratio OC/EC_bb = 5.4, from a PMF source apportionment available on the 
#       same filters
################################################################################
library(lhs)
library(ggplot2)
library(viridis)

#Load data
df <- read.csv("/Users/fra/Desktop/Tesi/Modellistica/dati.csv", header = T, sep = ";", 
               dec = ".", na.strings = "NaN")
df[,c(5:9)] <- df[,c(5:9)]/100

# Calculate the pmc of OC having TC, EC, pmc_EC, pmc_TC with propagation errors
df$fm_oc <- (df$fm_tc*df$tc-df$fm_ec*df$ec)*1/df$oc
df$errfm_oc <- sqrt((df$fm_tc*0.1*df$tc/df$oc)^2 + 
                      (df$fm_ec*0.1*df$ec/df$oc)^2 +
                      (df$tc*df$errfm_tc/df$oc)^2 + 
                      (df$ec*df$errfm_ec/df$oc)^2 + 
                      (df$fm_ec*df$ec/df$oc^2 - df$fm_tc*df$tc/df$oc^2)^2*(0.1*df$oc)^2)
df$errfm_oc_perc <- df$errfm_oc/df$fm_oc*100

# Reference values for the fm and for OC/EC from biomass
fm_bb <- 1.03
errfm_bb <- 0.02
fm_bio <- 0.997
OCsuEC_bb <- 5.4 
sigma <- 1
#-------------------------------------------------------------------------------
# Calculate EC from biomass burning using fm_bb with error
df$ec_bb <- df$ec*df$fm_ec/fm_bb
df$errec_bb <- sqrt((df$fm_ec*0.1*df$ec/fm_bb)^2 + 
                      (df$ec*df$errfm_ec/fm_bb)^2 +
                      (df$ec*df$fm_ec*errfm_bb/(fm_bb)^2)^2)

# Calculate EC from fossil fuel by difference
df$ec_ff <- df$ec-df$ec_bb
df$errec_ff <- sqrt((df$ec*0.1)^2 + (df$errec_bb)^2)

# Calculate biogenic OC
df$oc_bio <- (df$oc*df$fm_oc-df$oc_bb*fm_bb)/fm_bio

# Calculate fossil OC by difference
df$oc_ff <- df$oc - df$oc_bb - df$oc_bio

# Calculate the relative results in percentage
ECff <- 100*df$ec_ff/df$ec
ECbb <- 100*df$ec_bb/df$ec
OCff <- 100*df$oc_ff/df$oc
OCbb <- 100*df$oc_bb/df$oc
OCbio <- 100*df$oc_bio/df$oc

results <- data.frame(ECff, ECbb, OCff, OCbb, OCbio)
results_fm <- data.frame(cbind(df[,1], df[,c(5:10)]))

file_name <- "Results_fm.csv" 
write.csv(results_fm, file = file_name, row.names = FALSE)


################################################################################
# LHS
################################################################################
# sampling with latin hypercube sampling
n_samples <- 100
n_parameter <- 1
sampled_values01 <- randomLHS(n_samples, n_parameter)

# transform the sampled values in the desired range
min <- OCsuEC_bb - sigma
max <- OCsuEC_bb + sigma

# calculate the range width
range <- max - min

# transform the sampled values
sampled_values <- OCsuEC_bb + range * (sampled_values01 - 0.5)

#-------------------------------------------------------------------------------
sampled_values <- as.matrix(read.csv("/Users/fra/Desktop/Tesi/Modellistica/SampledOCsuEC.csv", header = T, sep = ";", 
                                     dec = ".", na.strings = "NaN"))

# Sensitivity of the results to the variation of OCsuEC_bb
path <- "/Users/fra/Desktop/Tesi/Modellistica/LHS/"
matrix <- matrix(nrow = n_samples, ncol = 5)
matrix <- cbind(sampled_values, matrix)
colnames(matrix) <- c("lhs", "ec_bb", "ec_ff", "oc_bb", "oc_bio", "oc_ff")
# loop on filters
for(j in 1:nrow(df)){
  # loop on the extracted values for the OC/EC ratio
  for(i in 1:n_samples){
    val <- sampled_values[i]
    matrix[i, 2] <- df$ec_bb[j]
    matrix[i, 3] <- df$ec_ff[j]
    matrix[i, 4] <- val*df$ec_bb[j] # OCbb
    matrix[i, 5] <- (df$oc[j]*df$fm_oc[j]-matrix[i, 4]*fm_bb)/fm_bio # OCbio
    matrix[i, 6] <- df$oc[j] - matrix[i, 4] - matrix[i, 5] # OCff
  }
  
  # Calculate relative apportionments in percentage
  matrix[,2] <- matrix[,2]*100/df$ec[j]
  matrix[,3] <- matrix[,3]*100/df$ec[j]
  matrix[,4] <- matrix[,4]*100/df$oc[j]
  matrix[,5] <- matrix[,5]*100/df$oc[j]
  matrix[,6] <- matrix[,6]*100/df$oc[j]
  
  file_name <- paste(df$filtro[j], ".csv", sep = "")
  write.csv(matrix, paste(path,file_name, sep = ""), row.names = FALSE)
  
  #Scatterplot
  grafico<-data.frame(matrice)
  scatter <- ggplot(grafico, aes(x=lhs))+
      theme_bw()+
      geom_point(aes(y = oc_bio, color = "OCbio")) +
      geom_point(aes(y = oc_ff, color = "OCff")) +
      geom_point(aes(y = oc_bb, color = "OCbb")) +
      labs(color = "Series") +
      scale_color_manual(values = c("red", "blue", "green"))+
      xlab("(OC/EC)_bb")+
      ylab("Apportionments in %")+
      ggtitle(df$filtro[j])
  ggsave(scatter,file=paste(df$filtro[j],".png",sep=""),
         path=path,dpi=500,
         limitsize = FALSE,width=25,height=20,units="cm")  
  
  png(paste(df$filtro[j],"box.png",sep=""))
  box <- boxplot(grafico[4:6], main=df$filtro[j])
  dev.off()
  print(box$stats)
}

################################################################################
# Comparisons with optical
################################################################################
# Calculate percentage contribution on tc of ff and bb

err_bbrel <- 100*df$errec_bb/df$ec
err_ffrel <- 100*df$errec_ff/df$ec
confronto_ottico <- cbind(ECbb, df$errec_bb, err_bbrel, ECff, df$errec_ff, err_ffrel)




#latex
xtable::xtable(risultati)
xtable::xtable(df[,c(1:8)])



readr