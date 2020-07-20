# set the working directory to where your input file is:
setwd("C:\\Users\\salvador\\Desktop\\Nam\\TMT")

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("MSstatsTMT")

library(MSstatsTMT)
# version 1.4.6 and 1.6.0 work for this tutorial
?MSstatsTMT

# we also use another two pacakges for data manipulation
# install.packages(c("tidyr", "dplyr"))
library(tidyr)
library(dplyr)

input <-  read.table("census-out-21063-NM_msstatsTMT.txt", sep = "\t", header = TRUE)

## 3.4 Preliminary check
# number of proteins:
length(unique(input$ProteinName)) 


# 4. Normalizing and summarizing data with proteinSummarization
# **! Always pay attention to the default options **  
# 
# After reading the datasets, `MSstatsTMT` performs 
# 
# * 1) logarithm transformation of `Intensity` column
# 
# * 2) global median normalization between channels
# 
# * 3) channel-level protein summarization
# 
# * 4) Local protein-level normalization with reference channel
# 
# To get started with this function, visit the help section of `proteinSummarization ` first: 
# 
?proteinSummarization 
 
## 4.1 Default summarization and normalization options

# `proteinSummarization` perform first (1) global_norm = FALSE: a global equal median normalization between channels will be performed to account for differences in labeling efficiency and other techinical artifacts.
# 
# Then, (2) method = "msstats": missing value imputation and protein summarization will be performed, which is robust parameter estimation by TMP (Tukey's median polish).
# 
# Last, (3) reference_norm = TRUE: normalization between MS runs. It needs at least one normalization channel in each MS run, annotated by 'Norm' in Condition column. If there are multiple normalization channels, all the normalization channels are averaged for the normalization. FALSE will not perform normalization step.
# 
# Below show the default for all options in proteinSummarization.
# 
 # use MSstats for protein summarization
quant <- proteinSummarization(input,
                                 method="msstats",
                                 global_norm=FALSE,
                                 reference_norm=FALSE,
                                 remove_norm_channel = FALSE,
                                 remove_empty_channel = TRUE)

# This table includes normalized channel-level log2 intensities. (column : Abundance)
# Now we have one summarized log2 intensities per Protein, Run and Channel.
# to see the top part of the table:
head(quant)
# if you already run proteinSummarization before and saved the quant object, dont save it and run load
save(quant, file='quant.rda')
# it will save you a lot of time
# load('quant.rda')


  # 5. Visualization of protein summaries
  
  # Profile plot is good visualization to check individual measurements. Each dot means one intensity per Run per Channel. Each panel represents one MS run and each dot within one panel is one channel within one Run. The dots are linked with line per feature. If line is disconnected, that means there is no value (missing value). Color means different peptides and charge stages. 

# if you have many MS runs, adjust width of plot (make wider)
# Profile plot for the normalized data
# this will take some time:
dataProcessPlotsTMT(data.peptide = input, # PSM-level data
                    data.summarization = quant, # protein-level data
                    type = 'ProfilePlot', # choice of visualization
                    width = 21,
                    height = 7,
                    address="pd_norm_") 

# `pd_norm_ProfilePlot.pdf` and `pd_norm_ProfilePlot_wSummarization.pdf` are generated in the current directory.

# Then, Let's go though profile plots to see overall quality of data. There are two pdfs for each protein, first is profile plot and second plot is profile plot with summarized and normalized data. `pd_norm_ProfilePlot.pdf` shows each peptide ions across runs and channels, grouped per condition. Each peptide has a different colour/type layout. `pd_norm_ProfilePlot_wSummarization.pdf` shows the same peptide ions in grey, with the values as summarized by the model overlayed in red.

# Instead of making all profile plots for all proteins, we can make plot for individual protein. Here is the example of background protein, `O75844`
 dataProcessPlotsTMT(data.peptide = input, # PSM-level data
                    data.summarization = quant, # protein-level data
                    type='ProfilePlot', # choice of visualization
                    width = 21,
                    height = 7,
                    which.Protein = 'A0A096MIV5',
                    address="pd_norm_A0A096MIV5_") 
 




# 7. Finding differentially abundant proteins across conditions

# After we summarized each protein's behavior across conditions and normalized the data between runs in `proteinSummarization` step, we are all set to compare protein changes between groups of conditions. Within MSstatsTMT we can do this with the `groupComparisonTMT` function, which takes as input the output of the `proteinSummarization` function. 

?groupComparisonTMT

## 7.1.Pairwise comparison

# If you want to make all the pairwise comparison,`MSstatsTMT` has an easy option for it. Setting `contrast.matrix = pairwise` compares all the possible pairs between two conditions.

# quant <- quant %>%  filter(Condition != "Norm")
# quant <- quant %>%  filter(Condition != "Control")
# this will not work for your data Name, because you don;t have more than 2 conditions
# test.pd.pairwise <- groupComparisonTMT(data = quant,
#                                        contrast.matrix = "pairwise",
#                                        moderated = TRUE, # do moderated t test
#                                        adj.method = "BH") # multiple comparison adjustment

# show the comparisons
# unique(test.pd.pairwise$Label)

# Let's check the output.
# colnames(test.pd.pairwise)

# Show test result
# Label : which comparison is used
# log2FC : estimated log2 fold change between two conditions (the contrast)
# adj.pvalue : adjusted p value
# head(test.pd.pairwise)

# Let's save the testing result as .csv file.
# save(test.pd.pairwise, file='pd.result.rda')
# write.csv(test.pd.pairwise, file='testResult_pd.csv')

## 7.2 Volcano plot

# Let's inspect the results to see what proteins are changing significantly between two concentrations.
# we need to use MSstats library now to make volcano plots and so on
library(MSstats)

# groupComparisonPlots(data=test.pd.pairwise, 
#                      type="VolcanoPlot", 
#                      logBase.pvalue=2, 
#                      ProteinName=TRUE, # only for small protein number
#                      address="pd_pairwise_")


## 7.3 Assign contrast matrix

# If you would like to compare some specific combination of conditions, you need to tell `groupComparisonTMT` the contrast of the conditions to compare. You can make your `contrast.matrix` in R in a text editor. We define our contrast matrix by adding a column for every condition. We add a row for every comparison we would like to make between groups of conditions.  

# **0** is for conditions we would like to ignore.
# **1** is for conditions we would like to put in the numerator of the ratio or fold-change.
# **-1** is for conditions we would like to put in the denumerator of the ratio or fold-change.

# If you have multiple groups, you can assign any group comparisons you are interested in.

# check unique conditions and check order of condition information
# In this case, 2 conditions, Control and NMcLPT
unique(quant$Condition)

# 'Norm' will be removed during tesing and should be not considered in the contrast
# comparison1<-matrix(c(-1,0,0,1),nrow=1) # 0.5-0.125
comparison<-matrix(c(-1,1),nrow=1) # NMcLTP-Control
# comparison2<-matrix(c(0,-1,1,0),nrow=1) # 0.667-0.5
# comparison<-rbind(comparison1, comparison2)
# Set the column names
# colnames(comparison)<- c("0.125", "0.5", "0.667", "1")
colnames(comparison)<- c("NMcLTP", "Control")
# Set the names of each row
row.names(comparison)<-c("NMcLTP-Control")

comparison

test.pd <- groupComparisonTMT(data = quant, 
                               contrast.matrix = comparison,
                               moderated = TRUE, # do moderated t test
                               adj.method = "BH") # multiple comparison adjustment
groupComparisonPlots(data=test.pd, 
                     type="VolcanoPlot", 
                     logBase.pvalue=2, 
                     ProteinName=TRUE, # only for small protein number
                     address="pd_")

# only valid for at least having 2 comparisons (so at least 3 conditions)
# groupComparisonPlots(data=test.pd, 
#                      type="Heatmap", 
#                      logBase.pvalue=2, 
#                      ProteinName=TRUE, # only for small protein number
#                      address="pd_")

# this works but it is only useful when having more than one comparison (so at least 3 conditions)
groupComparisonPlots(data=test.pd, 
                     type="ComparisonPlot", 
                     logBase.pvalue=2, 
                     ProteinName=TRUE, # only for small protein number
                     which.Protein = c("Q6PCU8","Q4QQT4"),
                     address="pd_")

# 8. msstatstmt.log and sessionInfo.txt

# These two files are important to keep the records of package versions and options in functions.
