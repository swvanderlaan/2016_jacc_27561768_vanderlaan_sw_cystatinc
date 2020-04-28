#####               CALCULATE & PLOT POWER OF MR CYSTATIN C FOR CVD             #####
#
# Version: POWER_MR.v3.20160419
#
# Last update: 2016-04-19
# Written by: Sander W. van der Laan (s.w.vanderlaan-2@umcutrecht.nl)
#                                                    
# Description: Script to calculate power for MR (for a binary trait) in 
#              various situation and plot it.
#
#
### -------------------------------------------------------------------------------------
### CLEAR THE BOARD
rm(list = ls())

### -------------------------------------------------------------------------------------
### GENERAL R SETUP 

### FUNCTION TO INSTALL PACKAGES, VERSION A -- This is a function found by 
### Sander W. van der Laan online from @Samir: 
### http://stackoverflow.com/questions/4090169/elegant-way-to-check-for-missing-packages-and-install-them
### Compared to VERSION 1 the advantage is that it will automatically check in both CRAN and Bioconductor
install.packages.auto <- function(x) { 
     x <- as.character(substitute(x)) 
     if (isTRUE(x %in% .packages(all.available = TRUE))) { 
          eval(parse(text = sprintf("require(\"%s\")", x)))
     } else { 
          # Update installed packages - this may mean a full upgrade of R, which in turn
          # may not be warrented. 
          #update.packages(ask = FALSE) 
          eval(parse(text = sprintf("install.packages(\"%s\", dependencies = TRUE)", x)))
     }
     if (isTRUE(x %in% .packages(all.available = TRUE))) { 
          eval(parse(text = sprintf("require(\"%s\")", x)))
     } else {
          source("http://bioconductor.org/biocLite.R")
          # Update installed packages - this may mean a full upgrade of R, which in turn
          # may not be warrented.
          #biocLite(character(), ask = FALSE) 
          eval(parse(text = sprintf("biocLite(\"%s\")", x)))
          eval(parse(text = sprintf("require(\"%s\")", x)))
     }
}
install.packages.auto("gtx")


# Create datestamp
Today = format(as.Date(as.POSIXlt(Sys.time())), "%Y%m%d")

### -------------------------------------------------------------------------------------
### SETUP ANALYSIS
# Assess where we are
getwd()
# Set locations
ROOT_loc = "/Users/swvanderlaan/PLINK/analyses/cystatinc/"
OUT_loc = "/Users/swvanderlaan/PLINK/analyses/cystatinc/power"
setwd(OUT_loc)

### -------------------------------------------------------------------------------------
### CALCULATE THE POWERRRR
# Reference: http://ije.oxfordjournals.org/content/43/3/922.abstract, 
# Appendix: http://ije.oxfordjournals.org/content/suppl/2014/03/06/dyu005.DC1/binpower050713app.pdf

# Take the inverse of the logit(Z) from some logistic regression
expit <- function(x) { 
     return(exp(x)/(1 + exp(x))) 
}

### Calculations of one specific scenario.

#
## Setting the things
#rsq   = 0.0275     # squared correlation of SNP with biomarker
##b1    = 0.2      # causal effect (log odds ratio per SD
#b1    = log(1.79) #                or log of OR per SD)
#sig = 0.05       # significance level (alpha)
#pow = 0.8        # power level       (1-beta)
#ratio = 1000/30000       # ratio of cases:controls = 1:ratio
#
## When you don't know your total sample size, but you know the power
#cat("Sample size required for", pow*100, "% power: \n")
#samplesize <- (qnorm(1 - sig/2) + qnorm(pow)) ^ 2/b1 ^ 2/rsq/(ratio/(1 + ratio))/(1/(1 + ratio))
#samplesize
#
## When you do have a sample size, but you don't know the power.
#n = 40000        #  Sample size
#
#cat("Power of analysis with", n, "participants: \n")
#power <- pnorm(sqrt(n*rsq*(ratio/(1 + ratio))*(1/(1 + ratio)))*b1 - qnorm(1 - sig/2))
#power*100


### Calculation of a range of scenarios
# SNP vs. CystC levels: 
beta = -0.091227

# Let's make a dataframe for multiple scenarios
# Sample sizes      Ratio               Cases     Controls
# HF = 36788        9.694186047         3,440     33,348
# MI = 41171        6.234405201         5,691     35,480
# CHD = 149804      2.478313365         43,068    106,736
# Stroke = 119036   6.092230696         16,784    102,252
# CVD = 252216      2.984958605         63,292    188,924

# Fully adjusted - Observational analysis Logistic & Cox
#Trait	analysis_type	                    OR	     low CI    high CI	beta	     low CI	high CI	cases	controls	n	     k	          ratio	     R2	     p
#CHD	     logistic; incident + prevalent	1.53	     1.29	     1.82	     0.427	0.255	0.599	43068	106736	149804	0.403500225	2.478313365	0.0275	0.0100
#Stroke	logistic; incident + prevalent	1.52	     1.22	     1.89	     0.418	0.199	0.637	16784	102252	119036	0.164143489	6.092230696	0.0275	0.0100
#HF	     logistic; incident + prevalent	4.77	     3.53	     6.44	     1.562	1.262	1.862	3440	     33348	36788	0.103154612	9.694186047	0.0275	0.0100
#MI	     logistic; incident + prevalent	1.21	     0.95	     1.53	     0.188	-0.052	0.428	5691	     35480	41171	0.160400225	6.234405201	0.0275	0.0100
#CVD	     logistic; incident + prevalent	1.82	     1.56	     2.13	     0.601	0.447	0.755	63292	188924	252216	0.335013021	2.984958605	0.0275	0.0100
#CHD	     cox; incident only	               1.93	     1.54	     2.43	     0.659	0.431	0.887	43068	106736	149804	0.403500225	2.478313365	0.0275	0.0100
#Stroke	cox; incident only	               1.50	     1.09	     2.05	     0.403	0.090	0.716	16784	102252	119036	0.164143489	6.092230696	0.0275	0.0100
#HF	     cox; incident only	               4.51	     3.26	     6.22	     1.505	1.182	1.828	3440	     33348	36788	0.103154612	9.694186047	0.0275	0.0100
#MI	     cox; incident only	               1.34	     0.93	     1.95	     0.296	-0.077	0.669	5691	     35480	41171	0.160400225	6.234405201	0.0275	0.0100
#CVD	     cox; incident only	               2.19	     1.82	     2.63	     0.782	0.597	0.966	63292	188924	252216	0.335013021	2.984958605	0.0275	0.0100

pval = 0.05/5
r2 = 0.0275
# Logistic
power_CHD = pnorm(sqrt(149804*r2*(2.478313365/(1 + 2.478313365))*(1/(1 + 2.478313365)))*log(1.29) - qnorm(1 - pval/2))
power_Stroke = pnorm(sqrt(119036*r2*(6.092230696/(1 + 6.092230696))*(1/(1 + 6.092230696)))*log(1.22) - qnorm(1 - pval/2))
power_HF = pnorm(sqrt(36788*r2*(9.694186047/(1 + 9.694186047))*(1/(1 + 9.694186047)))*log(3.53) - qnorm(1 - pval/2))
power_MI = pnorm(sqrt(41171*r2*(6.234405201/(1 + 6.234405201))*(1/(1 + 6.234405201)))*log(0.95) - qnorm(1 - pval/2))
power_CVD = pnorm(sqrt(252216*r2*(2.984958605/(1 + 2.984958605))*(1/(1 + 2.984958605)))*log(1.56) - qnorm(1 - pval/2))
# Cox
power_CHD_cox = pnorm(sqrt(149804*r2*(2.478313365/(1 + 2.478313365))*(1/(1 + 2.478313365)))*log(1.54) - qnorm(1 - pval/2))
power_Stroke_cox = pnorm(sqrt(119036*r2*(6.092230696/(1 + 6.092230696))*(1/(1 + 6.092230696)))*log(1.09) - qnorm(1 - pval/2))
power_HF_cox = pnorm(sqrt(36788*r2*(9.694186047/(1 + 9.694186047))*(1/(1 + 9.694186047)))*log(3.26) - qnorm(1 - pval/2))
power_MI_cox = pnorm(sqrt(41171*r2*(6.234405201/(1 + 6.234405201))*(1/(1 + 6.234405201)))*log(0.93) - qnorm(1 - pval/2))
power_CVD_cox = pnorm(sqrt(252216*r2*(2.984958605/(1 + 2.984958605))*(1/(1 + 2.984958605)))*log(1.82) - qnorm(1 - pval/2))

# CALCULATION FOR TABLE AND FIGURE
# Some presettings
causaleffect = seq(log(0.9), log(6.5), by = 0.001)
significance = 0.05/5 # 5 outcomes
rsquared = 0.0275 # CystC phenotypic variance explained by SNP rs911119
diseases = c("HF", "MI", "CHD", "IS", "CVD")
samplesizes = c(36788, 41171, 149804, 119036, 252216)
ratio = c(9.694186047, 6.234405201, 2.478313365, 6.092230696, 2.984958605)
#ratio = c(0.103154612, 0.160400225, 0.403500225, 0.164143489, 0.335013021)

power.mr.bin <- function(N, r2, ccratio, effect, pval, disease){
     cat("Calculating power...\n")
     power = pnorm(sqrt(N*r2*(ccratio/(1 + ccratio))*(1/(1 + ccratio)))*effect - qnorm(1 - pval/2))
     effectsize = effect
     ccratio = ccratio
     rsqrd = r2
     n = N
     P = pval
     disease = disease
     output = c(disease, n, ccratio, rsqrd, P, effectsize, power)
     cat("Disease........................................:", disease, "\n")
     cat("Total sample size..............................:", n, "\n")
     cat("Case/control ratio.............................:", ccratio, "\n")
     cat("The r^2 of the variant with the biomarker......:", round(rsqrd, 3), "\n")
     cat("Significance threshold.........................:", signif(P, 3), "\n")
     cat("Effect size....................................:", round(effectsize, 3), "\n")
     cat("Power..........................................:", round(power, 2), "\n\n")
     return(output)
     print(output)
}

# How to start an empty dataframe: http://stackoverflow.com/questions/10689055/create-an-empty-data-frame
power.df <- data.frame(matrix(NA, ncol = 7, nrow = 0))

# Power for diseases
for (d in 1:length(diseases)) {
     for (e in 1:length(causaleffect)) {
          power.df.TEMP <- data.frame(matrix(NA, ncol = 7, nrow = 0))
          power.df.TEMP[1,] = power.mr.bin(samplesizes[d], rsquared, ratio[d], causaleffect[e], significance, diseases[d])
          power.df = rbind(power.df, power.df.TEMP)
     }
}

# Edit the column names
colnames(power.df) = c("Disease", "N", "Case/Control Ratio", "r^2", "P threshold", "EffectSize", "Power")

# Correct the variable types
power.df$Disease <- as.character(power.df$Disease)
power.df$N <- as.numeric(power.df$N)
power.df$`Case/Control Ratio` <- as.numeric(power.df$`Case/Control Ratio`)
power.df$`r^2` <- as.numeric(power.df$`r^2`)
power.df$`P threshold` <- as.numeric(power.df$`P threshold`)
power.df$EffectSize <- as.numeric(power.df$EffectSize)
power.df$Power <- as.numeric(power.df$Power)

power.df.80 <- power.df[ which(power.df$Power >= 0.80 & power.df$Power < 0.81), ]

# Save the data
write.table(power.df, 
            paste0(OUT_loc,"/",Today,"_MR_CystC_Powerrr.txt"), 
            quote = FALSE , row.names = FALSE, col.names = TRUE, sep = "\t")
write.table(power.df.80, 
            paste0(OUT_loc,"/",Today,"_MR_CystC_Powerrr_at80.txt"), 
            quote = FALSE , row.names = FALSE, col.names = TRUE, sep = "\t")


### Plots for all data
#?axis
#?par
#?legend
#?plot
pdf(paste0(Today,"_MR_CystC_Power_Binary_ratio.pdf"), width = 7, height = 7, onefile = TRUE)
# get the number of colors to plot
powercol <- c("#F59D10", "#DB003F", "#9A3480", "#2F8BC9", "#9FC228")
power.df.HF <- power.df[ which(power.df$Disease == "HF"), ]
power.df.MI <- power.df[ which(power.df$Disease == "MI"), ]
power.df.CHD <- power.df[ which(power.df$Disease == "CHD"), ]
power.df.Stroke <- power.df[ which(power.df$Disease == "IS"), ]
power.df.CVD <- power.df[ which(power.df$Disease == "CVD"), ]

# start the actual plot for the first category
plot(exp(power.df.HF$EffectSize), power.df.HF$Power, type = "l", pch = 20, lty = 1, lwd = 2, col = "#F59D10", 
     #main = "Power to detect a causal effect", 
     xlab = "odds ratio\nper log2 increase of cystatin C",
     ylab = "power", 
     xlim = c(1, 1.80),
     ylim = c(0,1),
     axes = FALSE,
     cex.axis = 1,
     cex.lab = 1,
     cex.main = 1,
     bty = "n")

# add in the additional lines
lines(exp(power.df.MI$EffectSize), power.df.MI$Power, type = "l", pch = 20, lty = 1, lwd = 2, col = "#DB003F")
lines(exp(power.df.CHD$EffectSize), power.df.CHD$Power, type = "l", pch = 20, lty = 1, lwd = 2, col = "#9A3480")
lines(exp(power.df.Stroke$EffectSize), power.df.Stroke$Power, type = "l", pch = 20, lty = 1, lwd = 2, col = "#2F8BC9")
lines(exp(power.df.CVD$EffectSize), power.df.CVD$Power, type = "l", pch = 20, lty = 1, lwd = 2, col = "#9FC228")

# add in the axes
#axis(1, at = c(1.0, 1.10, 1.20, 1.30, 1.40, 1.50, 1.60), pos = 0, cex.axis = 0.75)
axis(1, at = c(1.0, 1.10, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70, 1.80), pos = 0, cex.axis = 1)
axis(2, at = c(0, 0.20, 0.40, 0.60, 0.80, 1.0), pos = 1, cex.axis = 1)


diseases = c("HF", "MI", "CHD", "IS", "CVD")
samplesizes = c(36788, 41171, 149804, 119036, 252216)
ratio = c(9.694186047, 6.234405201, 2.478313365, 6.092230696, 2.984958605)
#ratio = c(0.103154612, 0.160400225, 0.403500225, 0.164143489, 0.335013021)
# add in a legend
# Ratio
legend(1.50, 0.25, 
       legend = c("HF (36,788 - 1:9.7)", "MI (41,171 - 1:6.2)", "CHD (149,804 - 1:2.5)", "IS (119,036 - 1:6.1)", "CVD (252,216 - 1:3.0)"), 
       title = expression(paste("Outcome (",italic(N)," - ratio)")),
       cex = 0.75,
       lty = 1,
       lwd = 2,
       #pch = 20,
       col = c(powercol), 
       bty = "n")
# Kappa
#legend(1.50, 0.25, 
#       legend = c("HF (36,788 - 0.10)", "MI (41,171 - 0.16)", "CHD (149,804 - 0.40)", "Stroke (119,036 - 0.16)", "CVD (252,216 - 0.34)"), 
#       title = expression(paste("Outcome (",italic(N)," - ",kappa,")")),
#       cex = 0.75,
#       lty = 1,
#       lwd = 2,
#       #pch = 20,
#       col = c(powercol), 
#       bty = "n")

abline(h = 0.8, lwd = 2, lty = 2, col = "#595A5C")
dev.off()

### -------------------------------------------------------------------------------------
### SAVE THE DATA

save.image("~/PLINK/analyses/cystatinc/power/MR_Power_Binary.RData")


### -------------------------------------------------------------------------------------
### DEFINE THE COLORS

###  UtrechtSciencePark Colours Scheme
###
### Website to convert HEX to RGB: http://hex.colorrrs.com.
### For some functions you should divide these numbers by 255.
###
###  Color                         HEX       RGB                 CHR       MAF/INFO
### ---------------------------------------------------------------------------------------------------
###	yellow                        #FBB820 (251,184,32)	=>	1 		or 1.0 > INFO
###	gold                          #F59D10 (245,157,16)	=>	2		
###	salmon                        #E55738 (229,87,56) 	=>	3 		or 0.05 < MAF < 0.2 or 0.4 < INFO < 0.6
###	darkpink                      #DB003F ((219,0,63)		=>	4		
###	lightpink                     #E35493 (227,84,147)	=>	5 		or 0.8 < INFO < 1.0
###	pink                          #D5267B (213,38,123)	=>	6		
###	hardpink                      #CC0071 (204,0,113)		=>	7		
###	lightpurple                   #A8448A (168,68,138)	=>	8		
###	purple                        #9A3480 (154,52,128)	=>	9		
###	lavendel                      #8D5B9A (141,91,154)	=>	10		
###	bluepurple                    #705296 (112,82,150)	=>	11		
###	purpleblue                    #686AA9 (104,106,169)	=>	12		
###	lightpurpleblue               #6173AD (97,115,173)	=>	13		
###	seablue                       #4C81BF (76,129,191)	=>	14		
###	skyblue                       #2F8BC9 (47,139,201)	=>	15		
###	azurblue                      #1290D9 (18,144,217)	=>	16		 or 0.01 < MAF < 0.05 or 0.2 < INFO < 0.4
###	lightazurblue                 #1396D8 (19,150,216)	=>	17		
###	greenblue                     #15A6C1 (21,166,193)	=>	18		
###	seaweedgreen                  #5EB17F (94,177,127)	=>	19		
###	yellowgreen                   #86B833 (134,184,51)	=>	20		
###	lightmossgreen                #C5D220 (197,210,32)	=>	21		
###	mossgreen                     #9FC228 (159,194,40)	=>	22		or MAF > 0.20 or 0.6 < INFO < 0.8
###	lightgreen                    #78B113 (120,177,19)	=>	23/X
###	green                         #49A01D (73,160,29)		=>	24/Y
###	grey                          #595A5C (89,90,92)		=>	25/XY	or MAF < 0.01 or 0.0 < INFO < 0.2
###	lightgrey                     #A2A3A4	(162,163,164)	=> 	26/MT
### ---------------------------------------------------------------------------------------------------

