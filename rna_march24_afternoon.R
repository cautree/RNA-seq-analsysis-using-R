library(ggplot2)
getwd()
library(GGally)
setwd("/home/yanyan/Dropbox")
list.files()
rm(list=ls())
df=read.delim("seq.tsv",header = TRUE)
dim(df)
head(df,1)[1:5]
#dev.off()
colnames(df)


colnames(df)
rownames(df)
rownames(df) = df$gene
df$gene =NULL


# df_test=head(df)[1:5,1:5]
# df_test.t=as.data.frame(t(df_test))
# df_test.t



df.transposed=as.data.frame(t(df))
df.transposed$P
dim(df.transposed)





df.transposed$pd1_flag = ifelse(!is.na(df.transposed$PDCD1), 1,0)


# df.transposed$cd4flag = ifelse(!is.na(df.transposed$CD4),1,0)
# df.transposed$pd1flag = ifelse(!is.na(df.transposed$CXCL13),1,0)
# 
# par(mfrow = c(2, 2))
# boxplot(df.transposed$CXCL13, df.transposed$cd3eflag, main="cxcl13 on cells with or without cd3e ")
# boxplot(df.transposed$CXCL13, df.transposed$cd4flag, main="cxcl13 on cells with or without cd4 ")
# boxplot(df.transposed$CD4, df.transposed$cxcl13flag, main="cd4 on cells with or without cxcl13")
# boxplot(df.transposed$CD3E, df.transposed$cxcl13flag, main="cd3e on cells with or without cxcl13")
# 
# dev.off()

df.transposed$cd4andcd3e_flag = ifelse(!is.na(df.transposed$CD3E) |!is.na(df.transposed$CD4), 1,0)

cells_with_cd4andcd3e = subset(df.transposed, cd4andcd3e_flag==1)
cells_with_cd4andcd3e$pd1_flag = ifelse(!is.na(cells_with_cd4andcd3e$PDCD1), 1,0)
dim(cells_with_cd4andcd3e)
cells_without_cd4andcd3e = subset(df.transposed, cd4andcd3e_flag==0)
cells_without_cd4andcd3e$pd1_flag = ifelse(!is.na(cells_without_cd4andcd3e$PDCD1), 1,0)
dim(cells_without_cd4andcd3e)




# cells_with_cd4andcd3e_withpd1_ = subset(cells_with_cd4andcd3e, pd1_flag==1)
# cells_with_cd4andcd3e_withoutpd1_ = subset(cells_with_cd4andcd3e, pd1_flag==0)
# cells_without_cd4andcd3e_withpd1_ = subset(cells_without_cd4andcd3e, pd1_flag==1)
# cells_without_cd4andcd3e_withoutpd1_ = subset(cells_without_cd4andcd3e, pd1_flag==0)
# 
# dim(cells_with_cd4andcd3e_withpd1_)
# dim(cells_with_cd4andcd3e_withoutpd1_)
# dim(cells_without_cd4andcd3e_withpd1_)
# dim(cells_without_cd4andcd3e_withoutpd1_)

summary(cells_with_cd4andcd3e$PDCD1)

PD1_33 = quantile(cells_with_cd4andcd3e$PDCD1,0.33, na.rm = TRUE)
PD1_66= quantile(cells_with_cd4andcd3e$PDCD1,0.66, na.rm = TRUE)

cells_with_cd4andcd3e$PD1_bucket = cut(cells_with_cd4andcd3e$PDCD1, c(0,1,2,23))
cells_with_cd4andcd3e$PD1_bucket =as.factor(cells_with_cd4andcd3e$PD1_bucket)

cells_with_cd4andcd3e$PD1_bucket

#first PDCD1 NA
df_a = subset(cells_with_cd4andcd3e,is.na(PDCD1) )
dim(df_a)
df_aa = subset(cells_with_cd4andcd3e,!is.na(PDCD1) )
dim(df_aa )
#first PDCD1 <=1
df_b = subset(df_aa,PDCD1<=1 )
dim(df_b)

#first PDCD1 between 1 and 2
df_c = subset(df_aa, PDCD1>1 & PDCD1 <=2 )
dim(df_c)

#first PDCD1 bigger than 2
df_d = subset(df_aa, PDCD1>2)
dim(df_d)

cytokine = read.csv("Cytokine.csv", header = TRUE)
track=cytokine$tracking_id
track
track2 = c(as.character(track),"CD4","CXCL13","CD3E")
length(track2)
length(unique(track2))
dim(cytokine)

##cells_with_cd4andcd3e_withcxcl13
cytokine_with_a =df_a[,colnames(df_a) %in% track2]
dim(cytokine_with_a)

b=as.data.frame(t(cytokine_with_a))
dim(b)

d=as.data.frame(t(b))



getColMedian = function(x){
  x=x[!is.na(x)]
  return (median(x))
}

colMedian = as.data.frame(sapply(d,getColMedian))

dim(colMedian)
colMedian$genes = rownames(colMedian)
names(colMedian) = c("median","gene_name")
head(colMedian)


colMedian_ordered = colMedian[order(colMedian$median,decreasing = TRUE),]
colMedian_ordered= transform(colMedian_ordered, gene_name=reorder(gene_name, median) )

colMedian_ordered_top60_with = colMedian_ordered[1:60,]



ggplot(data=colMedian_ordered_top60_with, aes(x=gene_name, y=median, fill=gene_name)) +
  geom_bar(stat="identity", show.legend = FALSE)+
  coord_flip()

col_mean = as.data.frame(colMeans(d,na.rm = TRUE))
col_mean$genes = rownames(col_mean)
names(col_mean) = c("mean","gene_name")


colMean_ordered = col_mean[order(col_mean$mean,decreasing = TRUE),]
colMean_ordered = transform(colMean_ordered, gene_name=reorder(gene_name, mean) )
colMean_ordered_top60_with = colMean_ordered[1:60,]
colMean_ordered_top60_with
colMean_ordered_a = colMean_ordered

ggplot(data=colMean_ordered_top60_with, aes(x=gene_name, y=mean, fill=gene_name)) +
  geom_bar(stat="identity", show.legend = FALSE)+
  coord_flip()


#cells_with_cd4andcd3e_withoutcxcl13
cytokine_with_b =df_b[,colnames(df_b) %in% track2]

b=as.data.frame(t(cytokine_with_b))
dim(b)

d=as.data.frame(t(b))



getColMedian = function(x){
  x=x[!is.na(x)]
  return (median(x))
}

colMedian = as.data.frame(sapply(d,getColMedian))

dim(colMedian)
colMedian$genes = rownames(colMedian)
names(colMedian) = c("median","gene_name")
head(colMedian)


colMedian_ordered = colMedian[order(colMedian$median,decreasing = TRUE),]
colMedian_ordered= transform(colMedian_ordered, gene_name=reorder(gene_name, median) )
colMedian_ordered_top60_with = colMedian_ordered[1:60,]



ggplot(data=colMedian_ordered_top60_with, aes(x=gene_name, y=median, fill=gene_name)) +
  geom_bar(stat="identity", show.legend = FALSE)+
  coord_flip()

col_mean = as.data.frame(colMeans(d,na.rm = TRUE))
col_mean$genes = rownames(col_mean)
names(col_mean) = c("mean","gene_name")


colMean_ordered = col_mean[order(col_mean$mean,decreasing = TRUE),]
colMean_ordered = transform(colMean_ordered, gene_name=reorder(gene_name, mean) )
colMean_ordered_top60_with = colMean_ordered[1:60,]


colMean_ordered_b = colMean_ordered 


ggplot(data=colMean_ordered_top60_with, aes(x=gene_name, y=mean, fill=gene_name)) +
  geom_bar(stat="identity", show.legend = FALSE)+
  coord_flip()




#cells_without_cd4andcd3e_withcxcl13
cytokine_with_c =df_c[,colnames(df_c) %in% track2]

b=as.data.frame(t(cytokine_with_c))
dim(b)

d=as.data.frame(t(b))



getColMedian = function(x){
  x=x[!is.na(x)]
  return (median(x))
}

colMedian = as.data.frame(sapply(d,getColMedian))

dim(colMedian)
colMedian$genes = rownames(colMedian)
names(colMedian) = c("median","gene_name")
head(colMedian)


colMedian_ordered = colMedian[order(colMedian$median,decreasing = TRUE),]
colMedian_ordered= transform(colMedian_ordered, gene_name=reorder(gene_name, median) )
colMedian_ordered_top60_with = colMedian_ordered[1:60,]



ggplot(data=colMedian_ordered_top60_with, aes(x=gene_name, y=median, fill=gene_name)) +
  geom_bar(stat="identity", show.legend = FALSE)+
  coord_flip()

col_mean = as.data.frame(colMeans(d,na.rm = TRUE))
col_mean$genes = rownames(col_mean)
names(col_mean) = c("mean","gene_name")


colMean_ordered = col_mean[order(col_mean$mean,decreasing = TRUE),]
colMean_ordered = transform(colMean_ordered, gene_name=reorder(gene_name, mean) )
colMean_ordered_top60_with = colMean_ordered[1:60,]
colMean_ordered_top60_with

colMean_ordered_c = colMean_ordered

ggplot(data=colMean_ordered_top60_with, aes(x=gene_name, y=mean, fill=gene_name)) +
  geom_bar(stat="identity", show.legend = FALSE)+
  coord_flip()


#cells_without_cd4andcd3e_withoutcxcl13
cytokine_with_d =df_d[,colnames(df_d) %in% track2]

b=as.data.frame(t(cytokine_with_d))
dim(b)

d=as.data.frame(t(b))



getColMedian = function(x){
  x=x[!is.na(x)]
  return (median(x))
}

colMedian = as.data.frame(sapply(d,getColMedian))

dim(colMedian)
colMedian$genes = rownames(colMedian)
names(colMedian) = c("median","gene_name")
head(colMedian)


colMedian_ordered = colMedian[order(colMedian$median,decreasing = TRUE),]
colMedian_ordered= transform(colMedian_ordered, gene_name=reorder(gene_name, median) )
colMedian_ordered_top60_with = colMedian_ordered[1:60,]



ggplot(data=colMedian_ordered_top60_with, aes(x=gene_name, y=median, fill=gene_name)) +
  geom_bar(stat="identity", show.legend = FALSE)+
  coord_flip()

col_mean = as.data.frame(colMeans(d,na.rm = TRUE))
col_mean$genes = rownames(col_mean)
names(col_mean) = c("mean","gene_name")


colMean_ordered = col_mean[order(col_mean$mean,decreasing = TRUE),]
colMean_ordered = transform(colMean_ordered, gene_name=reorder(gene_name, mean) )
colMean_ordered_top60_with = colMean_ordered[1:60,]
colMean_ordered_top60_with

colMean_ordered_d = colMean_ordered

ggplot(data=colMean_ordered_top60_with, aes(x=gene_name, y=mean, fill=gene_name)) +
  geom_bar(stat="identity", show.legend = FALSE)+
  coord_flip()


merged1 = merge(colMean_ordered_a,colMean_ordered_b, by=0, all=TRUE)

merged2 = merge(colMean_ordered_c,colMean_ordered_d, by=0, all=TRUE)
merged_all = merge(merged1,merged2, by=0, all=TRUE)
dim(merged_all)

head(merged_all)
names(merged_all)
merged_all$Row.names=NULL
merged_all$gene_name.x.x=NULL
merged_all$gene_name.y.x=NULL
merged_all$gene_name.x.y=NULL
merged_all$gene_name.y.y=NULL
merged_all$Row.names.y=NULL

head(merged_all)
names(merged_all) = c("genename","A","B","C","D")

merged_all 

write.csv(merged_all, "heat_map.csv ")

head(merged_all)
dim(merged_all)

#names(merged_all) = c("genename","CD4/CD3E w CXCL13","CD4/CD3E w/o CXCL13","w/o CD4/CD3E w CXCL13","w/o CD4/CD3E w/o CXCL13")

rownames(merged_all) = merged_all$genename
merged_all$genename=NULL
rownames(merged_all)
merged_matrix = as.matrix(merged_all, dimnames=list(rownames(merged_all), colnames(merged_all)))
#merged_matrix[41,]

merged_matrix = merged_matrix[rowSums(!is.na(merged_matrix))!=0, colSums(!is.na(merged_matrix))!=0]

dim(merged_matrix)
merged_matrix1 = merged_matrix[1:40,]

merged_matrix1 

#get rid of 41
merged_matrix2 = merged_matrix[41:80,]

merged_matrix3 = merged_matrix[81:115,]


dev.new(width=5, height=10)
heatmap(merged_matrix1)

dev.new(width=5, height=10)
heatmap(merged_matrix2)

dev.new(width=5, height=10)
heatmap(merged_matrix3)

