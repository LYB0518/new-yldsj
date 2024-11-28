# 设置工作目录并加载GenABEL包
setwd("E:\\医疗大数据")
library(GenABEL)

# 加载GWAS数据
genodata <- load.gwaa.data(
  phenofile = "data/phenotype.dat",
  genofile = "data/genotype.raw",
  makemap = F,
  sort = F)

# 获取个体数量和标记数量
nids(genodata)
nsnps(genodata)

# 获取前5个标记的摘要信息
summary(genodata@gtdata)[1:5, ]

# 绘制第4号染色体上前500个标记的次要等位基因频率（MAF）
chr4.idx <- which(genodata@gtdata@chromosome == 4)
maf <- summary(genodata@gtdata)[chr4.idx[1:500],"Q.2"]
plot(maf, type = 'l', main="MAF on chr 4", xlab = "marker", ylab = "MAF", col = "slateblue")

# 汇总表型信息
summary(genodata@phdata)

# Fisher确切检验，用于测试性别与病例对照状态之间的关联
attach(phdata(genodata))
tab <- table(pheno, sex)
rownames(tab) <- c("control","case")
colnames(tab) <- c("male","female")
tab
fisher.test(tab)
detach(phdata(genodata))

# 质量控制过程
qc0 <- check.marker(genodata, call = 0.95, perid.call = 0.95,
                    maf = 1e-08, p.lev = 1e-05,
                    hweidsubset = genodata@phdata$pheno == 0)

# 移除质量控制中发现的问题标记和个体
genodata.qc0 <- genodata[qc0$idok, qc0$snpok]

# 计算基因组亲缘关系矩阵并进行多维尺度分析
autosomalMarkerNames <- autosomal(genodata.qc0)
genodata.qc0.gkin <- ibs(genodata.qc0[, autosomalMarkerNames], weight = "freq")
genodata.qc0.dist <- as.dist(0.5 - genodata.qc0.gkin)
genodata.qc0.pcs  <- cmdscale(genodata.qc0.dist, k=10)
genodata.qc0.pcs <- cbind(genodata.qc0.pcs, genodata.qc0@phdata$pheno)
plot(x=genodata.qc0.pcs[, 1], y=genodata.qc0.pcs[, 2], xlab = "PC 1", ylab = "PC 2",
     col=c("grey","slateblue")[as.factor(genodata.qc0.pcs[, 11])], main = "Genomic kinship")
legend(x=0.06,y=0.09,legend=c("control","case"), col=c("grey","slateblue"),pch=19)

# 进行GWAS分析
an <- qtscore(pheno, genodata.qc0, trait="binomial")

# 绘制曼哈顿图
plot(an, col=c("olivedrab","slateblue"), cex=.5, main="Manhattan plot")

# 计算λ值并绘制Q-Q图
estlambda(an[, "P1df"], plot=TRUE)

# 校正多重检验
an.pca <- qtscore(pheno~genodata.qc0.pcs[, 1:10], genodata.qc0, trait="binomial")
plot(an.pca, col=c("olivedrab","slateblue"), cex=.5, main="Manhattan plot")
estlambda(an.pca[, "P1df"], plot=TRUE)

# 总结顶部SNPs并添加显著性线到曼哈顿图
plot(an.pca, col=c("olivedrab","slateblue"), cex=.5, main="Manhattan plot")
abline(h = -log10(5e-8), col = "red")
summary(an.pca)

# 计算Bonferroni校正阈值
bonferroni <- -log10(0.05/nsnps(genodata.qc0))
plot(an.pca, col=c("olivedrab","slateblue"), cex=.5, main="Manhattan plot")
abline(h = bonferroni, col = "red")

# 执行置换测试
an.pca.per <- qtscore(pheno~genodata.qc0.pcs[, 1:10], genodata.qc0, trait="binomial", times = 10000)