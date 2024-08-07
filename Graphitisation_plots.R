################################################################################
# What does it do? 
# On each iteration, it plots and saves in .png format the pressure time-serie 
# observed in the graphitisation chamber during the Bosch reaction, for the 
# sample in question.

# Usage:
# 1) Select input path from which to read the data.
# 2) Select output path where to save the plots (remember to create the folder, 
#    if necessary).
# 3) Choose a save name for the plot (Note: if not updated, any old plots 
#    already present in the same folder will be overwritten).
# 4) Execute.
################################################################################
library(ggplot2)

# Paths dir_list[i] and output to update
output <- "/Users/fra/Desktop/Tesi/Preparazione_campioni_14C_21nov22/Plot"

main_dir <- "/Users/fra/Desktop/Tesi/Preparazione_campioni_14C_21nov22"
dir_list <- list.dirs(main_dir, recursive = FALSE)

check <- paste(main_dir, "/RH", sep = "")
for (i in 1:length(dir_list)) {
  if (substring(dir_list[i], 1, nchar(check)) == check) {
    
    # Read pressure and Peltier temperature files
    sample_name <- substring(dir_list[i], nchar(main_dir)+1, nchar(dir_list[i]))
    print(i)
    print(sample_name)
    pressure <- read.table(paste(dir_list[i],"/Trasduttore_pressione.txt",sep=""), 
                           header = T, sep = "", dec = ",", na.strings = "NaN")
    peltier <- read.table(paste(dir_list[i],"/Temperatura_Peltier.txt",sep=""), 
                          header = F, sep = "", dec = ",", na.strings = "NaN")
    
    # Remove date column and fix header
    pressure <- as.data.frame(pressure[,-1])
    colnames(pressure) <- c("time", "p")
    if (nchar(pressure$time[1]) > 5) {
      pressure$time <- as.POSIXct(strptime(pressure$time, format = "%H:%M:%S"))
      start <- as.POSIXct(strptime(peltier[2,2], format = "%H:%M:%S"))
    } else {
      pressure$time <- as.POSIXct(strptime(pressure$time, format = "%H:%M"))
      start <- as.POSIXct(strptime(peltier[2,2], format = "%H:%M"))
    }
    
    # Select the pressures measured after Peltier ignition
    plot_data <- pressure[which(pressure$time > start), ]
    neg <- which(sign(pressure$p) == -1)
    if (length(neg) > 0) {
      plot_data <- plot_data[neg[length(neg)]:length(plot_data$p), ] 
    }
    max_index <- which.max(plot_data$p)
    if (length(max_index:length(plot_data$p)) != 0) {
      
      plot_data <- plot_data[c(max_index:length(plot_data$p)), ]
      # Plot structure
      plot <- ggplot()+
        geom_point(data = plot_data, aes(x = plot_data[,"time"], y = plot_data[,"p"]))+
        theme_bw()+
        xlab("Time [hh:mm]")+
        ylab("P [mbar]")+
        ggtitle(sample_name)
      ggsave(plot, file = paste(sample_name, ".png", sep = ""), path = output, dpi = 500, limitsize = FALSE)
      
    }
    
  }
}