# assumes files with pattern "*T_.bg" in the current directory.
# assumption of chromosome names : excluding chromosomes containing "[_YM]" in their name for loess smoothing

args = commandArgs(TRUE)
path_to_files = args[1]
merged_bedgraph_file = args[2] # "merge_RT.txt"
bg_files = args[3] # c("my_sample_1_T_.bg","my_sample_2_T_.bg")


setwd(path_to_files)
library(preprocessCore) 

#  Import the merged bedgraph files: 
merge<-read.table(merged_bedgraph_file, header=FALSE)
colnames(merge)<-c(c("chr","start","end"), list.files(path=".",pattern="*T_.bg"))
merge_values<-as.matrix(merge[,4:ncol(merge)]) 

# Set the datasets to use for quantile normalization (bold names have to be adapted): 
# normalization on all datasets: 
ad<-stack(merge[,4:ncol(merge)])$values 

if (length(bg_files)==1){
  # normalisation on one datasets:
  ad<-merge[, bg_files[1]]
} else {
  # normalisation on multiple datasets (You can add as many datasets as you want):
  ad<-stack(merge[, bg_files])$values 
}

# Normalise the data: 
norm_data<-normalize.quantiles.use.target(merge_values,ad)
merge_norm<-data.frame(merge[,1:3],norm_data)
colnames(merge_norm)<-colnames(merge) 

# Register the quantile normalized data into bedgraph files: 
for(i in 4:ncol(merge_norm)){
  write.table( merge_norm[complete.cases(merge_norm[,i]), c(1,2,3,i)], gsub(".bg", "qnorm.bedGraph", colnames(merge_norm)[i]), sep="\t",row.names=FALSE, quote=FALSE, col.names=FALSE)
} 

# Select the chromosome for Loess smoothing (You can modify the pattern option to select different chromosomes. Here, the pattern selects all chromosomes except chromosomes containing "_" in their name (which can be problematics for Loess smoothing), and the chromosomes Y and M (mitochondrial)):
chrs=grep(levels(merge_norm$chr),pattern="[_YM]",invert=TRUE,value=TRUE) 

# Initialise an R-list to stock your datasets: 
AllLoess=list() 

# Perform Loess smoothing (This smoothing is similar to the Loess smoothing used by Ryba et al. 18 for repli-chip analysis). The window size used for span value (in bold) can be adapted. We are using here 300kb windows but you can increase this value to increase the smoothing. 
for(i in 1:(ncol(merge_norm)-3)){
  AllLoess[[i]]=data.frame();
  cat("Current dataset:", colnames(merge_norm)[i+3], "\n");
  for(Chr in chrs){
    RTb=subset(merge_norm, merge_norm$chr==Chr);
    lspan=300000/(max(RTb$start)-min(RTb$start));
    cat("Current chrom:", Chr, "\n");
    RTla=loess(RTb[,i+3] ~ RTb$start, span=lspan);
    RTl=data.frame(c(rep(Chr,times=RTla$n)), RTla$x, merge_norm[which(merge_norm$chr==Chr & merge_norm$start %in% RTla$x),3],RTla$fitted); colnames(RTl)=c("chr","start","end",colnames(RTb)[i+3]);
    if(length(AllLoess[[i]])!=0){
      AllLoess[[i]]=rbind(AllLoess[[i]],RTl)
    };
    if(length(AllLoess[[i]])==0){
      AllLoess[[i]]=RTl
    }
  }
} 

# Register the Loess smoothed data into bedgraph files:
for(i in 1:length(AllLoess)){
  write.table(AllLoess[[i]][complete.cases(AllLoess[[i]]),], gsub(".bg","Loess.bedGraph", colnames(AllLoess[[i]]))[4], sep="\t", row.names=FALSE, quote=FALSE, col.names=FALSE)
}
