library('DESeq2')

count1 = read.csv('simAR_sizeFactor_data1_1.txt', header=F);
count2 = read.csv('simAR_sizeFactor_data2_1.txt', header=F);

t1=cbind(count1[,1:4], count2[, 1:4])
t2=cbind(count1[,5:8], count2[, 5:8])
t3=cbind(count1[,9:12], count2[, 9:12])
t4=cbind(count1[,13:16], count2[, 13:16])
t5=cbind(count1[,17:20], count2[, 17:20])
coln = c( 'c1_rep1', 'c1_rep2', 'c1_rep3', 'c1_rep4', 'c2_rep1', 'c2_rep2', 'c2_rep3', 'c2_rep4')
rown = paste('g', 1:1000, sep='')
colnames(t1)=coln
colnames(t2)=coln
colnames(t3)=coln
colnames(t4)=coln
colnames(t5)=coln
rownames(t1)=rown
rownames(t2)=rown
rownames(t3)=rown
rownames(t4)=rown
rownames(t5)=rown

coldata = read.table('coldata.txt', header=T, row.names=1)

dds1 = DESeqDataSetFromMatrix(countData = t1 
                             , colData = coldata,
                             , design = ~ condition
                             ) 
dds1 = DESeq(dds1)
res1 = results(dds1)

dds2 = DESeqDataSetFromMatrix(countData = t2 
                             , colData = coldata,
                             , design = ~ condition
                             ) 
dds2 = DESeq(dds2)
res2 = results(dds2)
dds3 = DESeqDataSetFromMatrix(countData = t3 
                             , colData = coldata,
                             , design = ~ condition
                             ) 
dds3 = DESeq(dds3)
res3 = results(dds3)
dds4 = DESeqDataSetFromMatrix(countData = t4 
                             , colData = coldata,
                             , design = ~ condition
                             ) 
dds4 = DESeq(dds4)
res4 = results(dds4)
dds5 = DESeqDataSetFromMatrix(countData = t5 
                             , colData = coldata,
                             , design = ~ condition
                             ) 
dds5 = DESeq(dds5)
res5 = results(dds5)


pmat = cbind(res1$pvalue, res2$pvalue, res3$pvalue, res4$pvalue, res5$pvalue)

minp = apply(pmat, 1, min)
write.table(file='min_AR.txt', minp, row.names=F, quote=F, col.names=F)

sump = apply(pmat, 1, mean)
write.table(file='sum_AR.txt', sump, row.names=F, quote=F, col.names=F)
