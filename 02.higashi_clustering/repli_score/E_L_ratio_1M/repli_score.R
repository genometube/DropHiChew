library(preprocessCore)
library(ggplot2)

#args = commandArgs(trailingOnly=TRUE)
#in_file=args[1]
#in_file="E_L_293_replicate_1_20000.bedgraph"
in_file="E_L_293_replicate_1_1000000.bedgraph"

merge<-read.table('merge_RT.txt', header=TRUE)
merge_values<-as.matrix(merge[,4:ncol(merge)])

ad<-merge[,in_file]
# options(mc.cores = 1)

norm_data<-normalize.quantiles.use.target(merge_values,ad)
merge_norm<-data.frame(merge[,1:3],norm_data)
colnames(merge_norm)<-colnames(merge)


for(i in 4:ncol(merge_norm)){write.table(merge_norm[complete.cases(merge_norm[,i]),
c(1,2,3,i)], gsub(".bedgraph" , "qnorm.bedGraph", colnames(merge_norm)[i]), sep= "\t" ,row.names=FALSE, quote=FALSE, col.names = FALSE)}

table(merge_norm$chr)
merge_norm$chr <- factor(merge_norm$chr)
chrs=grep(levels(merge_norm$chr),pattern= "chr",invert=FALSE,value=TRUE)

AllLoess=list()

for(i in 1:(ncol(merge_norm)-3)){
  AllLoess[[i]]=data.frame();
  cat("Current dataset:", colnames(merge_norm)[i+3], "\n");
  
  for(Chr in chrs){
    RTb=subset(merge_norm, merge_norm$chr==Chr);
    lspan=300000/(max(RTb$start)-min(RTb$start));
    cat("Current chrom:" , Chr, "\n");
    RTla=loess(RTb[,i+3] ~ RTb$start, span=lspan);
    RTl=data.frame(c(rep(Chr,times=RTla$n)), RTla$x, merge_norm[which( merge_norm$chr==Chr & merge_norm$start %in% RTla$x),3],RTla$fitted);
    colnames(RTl)=c("chr" , "start" , "end" ,colnames(RTb)[i+3]);
    if(length(AllLoess[[i]])!=0){
      AllLoess[[i]]=rbind(AllLoess[[i]],RTl)};
    if(length(AllLoess[[i]])==0){
      AllLoess[[i]] = RTl}}}

for(i in 1:length(AllLoess)){write.table(AllLoess[[i]][complete.cases(AllLoess[[i]]),], 
                                         gsub(".bedgraph" , "Loess.bedGraph" , colnames(AllLoess[[i]]))[4], sep= "\t" , row.names=FALSE, quote=FALSE, col.names = FALSE)}



normed<-merge_norm
normed$locus<-as.numeric(rownames(normed))
ggplot(normed,aes(x=locus,y=normed$E_L_293_replicate_1_20000.bedgraph))+geom_line()+xlim(5000,10000)+ylim(-5,5)

loess_corrected<-AllLoess[[1]]
loess_corrected$locus<-as.numeric(rownames(loess_corrected))
ggplot(loess_corrected,aes(x=locus,y=loess_corrected$E_L_293_replicate_1_20000.bedgraph))+geom_line()+xlim(5000,10000)+ylim(-5,5)
    
