library(vegan)
defaultdir <- getwd()
##读取数据
#读入物种数据，以细菌 OTU 水平丰度表为例
otu <- read.delim(paste(defaultdir,"/","otu_table.txt",sep=""), row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
otu <- data.frame(t(otu))

#读取环境数据
env <- read.delim(paste(defaultdir,"/","env_table.txt",sep=""), row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)

#############################
##CCA，通过 vegan 包的 CCA 函数 cca() 执行，详情 ?cca

#调用格式 1，cca(Y, X, W)：Y，响应变量矩阵；X，解释变量矩阵；W，协变量矩阵（偏 CCA 时使用）
#这里无协变量矩阵，所以直接输入响应变量矩阵和解释变量矩阵
otu_cca <- cca(otu, env)

#或者格式 2，cca(Y~var1+var2+var3+factorA+var2*var3+Condition(var4))
#var1、var2 等，数值型解释变量；factorA，因子型解释变量；var2*var3，考虑变量间的交互作用；Condition(var4)，变量 4 作为协变量
#Y~. 是 Y~var1+var2+...+varn 的简写，不涉及交互作用及协变量
otu_cca <- cca(otu~., env)

#############################
##CCA 结果解读
#查看统计结果信息，以 I 型标尺为例
otu_cca.scaling1 <- summary(otu_cca, scaling = 1)
otu_cca.scaling1

##作图查看排序结果，详情 ?plot.cca
#三序图，包含 I 型标尺和 II 型标尺，样方坐标展示为使用物种加权计算的样方得分
par(mfrow = c(1, 2))
plot(otu_cca, scaling = 1, main = 'I 型标尺', display = c('wa', 'sp', 'cn'))
plot(otu_cca, scaling = 2, main = 'II 型标尺', display = c('wa', 'sp', 'cn'))

#I 型标尺，样方得分是物种得分加权平均，故排序图中物种通常分布在样方范围之外
#II 型标尺，物种得分是样方得分的平均值，并按出现该物种的所有样方中该物种丰度加权，故排序图中样方通常分布在物种范围之外

#比较分别使用物种加权计算的样方坐标以及拟合的样方坐标的差异
#隐藏物种，以 I 型标尺为例展示双序图
par(mfrow = c(1, 2))
plot(otu_cca, scaling = 1, main = 'I 型标尺，加权', display = c('wa', 'cn'))
plot(otu_cca, scaling = 1, main = 'I 型标尺，拟合', display = c('lc', 'cn'))

##CCA 结果提取
#scores() 提取排序得分（坐标），以 I 型标尺为例，前四轴为例
#使用物种加权和计算的样方得分
otu_cca_site.scaling1 <- scores(otu_cca, choices = 1:4, scaling = 1, display = 'wa')	
#物种变量（响应变量）得分
otu_cca_sp.scaling1 <- scores(otu_cca, choices = 1:4, scaling = 1, display = 'sp')
#环境变量（解释变量）得分
otu_cca_env.scaling1 <- scores(otu_cca, choices = 1:4, scaling = 1, display = 'bp')

#或者在 summary() 后提取，以 I 型标尺为例，前四轴为例
otu_cca.scaling1 <- summary(otu_cca, scaling = 1)
#使用物种加权和计算的样方得分
otu_cca_site.scaling1 <- otu_cca.scaling1$site[ ,1:4]
#物种
otu_cca_sp.scaling1 <- otu_cca.scaling1$species[ ,1:4]
#环境
otu_cca_env.scaling1 <- otu_cca.scaling1$biplot[ ,1:4]

#若需要输出在本地
#样方
write.table(data.frame(otu_cca_site.scaling1), 'otu_cca_site.scaling1.txt', sep = '\t', col.names = NA, quote = FALSE)
#物种
write.table(data.frame(otu_cca_sp.scaling1), 'otu_cca_sp.scaling1.txt', sep = '\t', col.names = NA, quote = FALSE)
#环境
write.table(data.frame(otu_cca_env.scaling1), 'otu_cca_env.scaling1.txt', sep = '\t', col.names = NA, quote = FALSE)

#不建议直接在原始数据集中提取，因为这样提取的坐标数值未经标尺缩放处理，不利于反映生物学问题
#otu_cca$CCA$u[ ,1:4]
#otu_cca$CCA$v[ ,1:4]
#otu_cca$CCA$biplot[ ,1:4]

##coef() 提取 CCA 典范系数
cca_coef <- coef(otu_cca)

#############################
##R2 校正
#RsquareAdj() 提取 R2，详情 ?RsquareAdj() 
r2 <- RsquareAdj(otu_cca)
otu_cca_noadj <- r2$r.squared	#原始 R2
otu_cca_adj <- r2$adj.r.squared	#校正后的 R2

#关于约束轴承载的特征值或解释率，应当在 R2 校正后重新计算
otu_cca_exp_adj <- otu_cca_adj * otu_cca$CCA$eig/sum(otu_cca$CCA$eig)
otu_cca_eig_adj <- otu_cca_exp_adj * otu_cca$tot.chi

##置换检验
#所有约束轴的置换检验，即全局检验，基于 999 次置换，详情 ?anova.cca
otu_cca_test <- anova.cca(otu_cca, permutations = 999)

#各约束轴逐一检验，基于 999 次置换
otu_cca_test_axis <- anova.cca(otu_cca, by = 'axis', permutations = 999)

#p 值校正（Bonferroni 为例）
otu_cca_test_axis$`Pr(>F)` <- p.adjust(otu_cca_test_axis$`Pr(>F)`, method = 'bonferroni')

#############################
##变量选择
#计算方差膨胀因子，详情 ?vif.cca
vif.cca(otu_cca)

#前向选择，这里以 ordiR2step() 的方法为例，基于 999 次置换检验，详情 ?ordiR2step
otu_cca_forward_pr <- ordiR2step(cca(otu~1, env), scope = formula(otu_cca), R2scope = TRUE, direction = 'forward', permutations = 999)

#以 otu_cca 和 otu_cca_forward_pr 为例，简要绘制双序图比较变量选择前后结果
par(mfrow = c(1, 2))
plot(otu_cca, scaling = 1, main = '原始模型，I 型标尺', display = c('wa', 'cn'))
plot(otu_cca_forward_pr, scaling = 1, main = '前向选择后，I 型标尺', display = c('wa', 'cn'))

#细节部分查看
summary(otu_cca_forward_pr, scaling = 1)

#比较选择前后校正后 R2 的差异，详情 ?RsquareAdj
#可以看到变量选择后，尽管去除了很多环境变量，但总 R2 并未损失很多
RsquareAdj(otu_cca)$adj.r.squared
RsquareAdj(otu_cca_forward_pr)$adj.r.squared

#所有约束轴的全局检验，999 次置换，详情 ?anova.cca
otu_cca_forward_pr_test <- anova.cca(otu_cca_forward_pr, permutations = 999)

#各约束轴逐一检验，999 次置换
otu_cca_forward_pr_test_axis <- anova.cca(otu_cca_forward_pr, by = 'axis', permutations = 999)

#p 值校正（Bonferroni 为例）
otu_cca_forward_pr_test_axis$`Pr(>F)` <- p.adjust(otu_cca_forward_pr_test_axis$`Pr(>F)`, method = 'bonferroni')

##提取或输出变量选择后的排序坐标
#提取方式和上文一致，这里通过 scores() 提取，以 I 型标尺为例，前两轴为例

#使用物种加权和计算的样方得分
otu_cca_forward_pr_site.scaling1 <- scores(otu_cca_forward_pr, choices = 1:2, scaling = 1, display = 'wa')
write.table(data.frame(otu_cca_forward_pr_site.scaling1), 'otu_cca_forward_pr_site.scaling1.txt', sep = '\t', col.names = NA, quote = FALSE)
#物种变量（响应变量）得分
otu_cca_forward_pr_sp.scaling1 <- scores(otu_cca_forward_pr, choices = 1:2, scaling = 1, display = 'sp')
write.table(data.frame(otu_cca_forward_pr_sp.scaling1), 'otu_cca_forward_pr_sp.scaling1.txt', sep = '\t', col.names = NA, quote = FALSE)
#环境变量（解释变量）得分
otu_cca_forward_pr_env.scaling1 <- scores(otu_cca_forward_pr, choices = 1:2, scaling = 1, display = 'bp')
write.table(data.frame(otu_cca_forward_pr_env.scaling1), 'otu_cca_forward_pr_env.scaling1.txt', sep = '\t', col.names = NA, quote = FALSE)

#############################
##变差分解 varpart()
#以两组环境变量（DOC、AP+AK）为例，运行变差分解，详情 ?varpart
#参数 chisquare = TRUE，执行 CCA 的变差分解；chisquare = FALSE，执行 RDA 的变差分解
#注：如果 varpart() 不支持 CCA，请更新 R 版本（如 R3.6 的 vegan）
otu_cca_vp <- varpart(otu, env['DOC'], env[c('AP', 'AK')], chisquare = TRUE)
otu_cca_vp

plot(otu_cca_vp, digits = 2, Xnames = c('DOC', 'AP+AK'), bg = c('blue', 'red'))

#解释变差的置换检验，以 DOC 所能解释的全部变差为例；999 次置换
anova.cca(cca(otu~DOC, env), permutations = 999)

#若考虑 DOC 单独解释的变差部分，需将其它变量作为协变量；999 次置换
anova.cca(cca(otu~DOC+Condition(AP+AK), env), permutations = 999)

#############################
##以上述前向选择后的简约模型 otu_cca_forward_pr 为例作图展示前两轴

#计算校正 R2 后的约束轴解释率
exp_adj <- RsquareAdj(otu_cca_forward_pr)$adj.r.squared * otu_cca_forward_pr$CCA$eig/sum(otu_cca_forward_pr$CCA$eig)
cca1_exp <- paste('CCA1:', round(exp_adj[1]*100, 2), '%')
cca2_exp <- paste('CCA2:', round(exp_adj[2]*100, 2), '%')

#plot() 作图，详情 ?plot.cca
#样方展示为点，物种暂且展示为“+”，环境变量为向量
par(mfrow = c(1, 2))

plot(otu_cca_forward_pr, type = 'n', display = c('wa', 'cn'), choices = 1:2, scaling = 1, main = 'I型标尺，双序图', xlab = cca1_exp, ylab = cca2_exp)
points(otu_cca_forward_pr, choices = 1:2, scaling = 1, display = 'wa', pch = 19, col = c(rep('red', 9), rep('orange', 9), rep('green3', 9)), cex = 1)
text(otu_cca_forward_pr, choices = 1:2, scaling = 1, display = 'cn', col = 'blue', cex = 0.8)

plot(otu_cca_forward_pr, type = 'n', display = c('wa', 'sp', 'cn'), choices = 1:2, scaling = 1, main = 'I型标尺，三序图', xlab = cca1_exp, ylab = cca2_exp)
points(otu_cca_forward_pr, choices = 1:2, scaling = 1, display = 'sp', pch = 3, col = 'gray', cex = 1)
points(otu_cca_forward_pr, choices = 1:2, scaling = 1, display = 'wa', pch = 19, col = c(rep('red', 9), rep('orange', 9), rep('green3', 9)), cex = 1)
text(otu_cca_forward_pr, choices = 1:2, scaling = 1, display = 'cn', col = 'blue', cex = 0.8)

#下面是 ggplot2 方法
#提取样方和环境因子排序坐标，前两轴，I 型标尺
otu_cca_forward_pr.scaling1 <- summary(otu_cca_forward_pr, scaling = 1)
otu_cca_forward_pr.site <- data.frame(otu_cca_forward_pr.scaling1$sites)[1:2]
otu_cca_forward_pr.env <- data.frame(otu_cca_forward_pr.scaling1$biplot)[1:2]

#手动添加分组
otu_cca_forward_pr.env$name <- rownames(otu_cca_forward_pr.env)
otu_cca_forward_pr.site$name <- rownames(otu_cca_forward_pr.site)
otu_cca_forward_pr.site$group <- c(rep('A', 9), rep('B', 9), rep('C', 9))

#ggplot2 作图
library(ggplot2)
library(ggrepel)	#用于 geom_label_repel() 添加标签

p <- ggplot(otu_cca_forward_pr.site, aes(CCA1, CCA2)) +
geom_point(aes(color = group)) +
stat_ellipse(aes(color = group), level = 0.95, show.legend = FALSE, linetype = 2) +
scale_color_manual(values = c('red', 'orange', 'green3')) +
theme(panel.grid.major = element_line(color = 'gray', size = 0.1), panel.background = element_rect(color = 'black', fill = 'transparent'), 
	legend.title = element_blank(), legend.key = element_rect(fill = 'transparent'), plot.title = element_text(hjust = 0.5)) + 
labs(x = cca1_exp, y = cca2_exp, title = 'I型标尺，双序图') +
geom_vline(xintercept = 0, color = 'gray', size = 0.5) + 
geom_hline(yintercept = 0, color = 'gray', size = 0.5) +
geom_segment(data = otu_cca_forward_pr.env, aes(x = 0, y = 0, xend = CCA1, yend = CCA2), arrow = arrow(length = unit(0.2, 'cm')), size = 0.3, color = 'blue') +
geom_text(data = otu_cca_forward_pr.env, aes(CCA1 * 1.2, CCA2 * 1.2, label = name), color = 'blue', size = 3) +
geom_label_repel(aes(label = name, color = group), size = 3, box.padding = unit(0.5, 'lines'), show.legend = FALSE)

p

