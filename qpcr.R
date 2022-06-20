
# This file contains the code used to analyze qRT-PCR data for Berger & Tarrant 2022.


library(ggplot2)
library(plotrix)
library(dplyr)
library(tidyr)
library(lomb)
library(reshape2)
library(permute)

setwd("F:/qpcr/")

#normalization

df = read.csv("Norm_constant.csv")
#df = read.csv("Norm_cycle.csv")
names(df)[1] = "tissue"

df = melt(df)

df$time = 0
df[startsWith(as.character(df$tissue), "T1B"),]$time = 1
df[startsWith(as.character(df$tissue), "T2B"),]$time = 2
df[startsWith(as.character(df$tissue), "T3B"),]$time = 3
df[startsWith(as.character(df$tissue), "T4B"),]$time = 4
df[startsWith(as.character(df$tissue), "T5B"),]$time = 5
df[startsWith(as.character(df$tissue), "T6B"),]$time = 6
df[startsWith(as.character(df$tissue), "T7B"),]$time = 7
df[startsWith(as.character(df$tissue), "T8B"),]$time = 8
df[startsWith(as.character(df$tissue), "T9B"),]$time = 9
df[startsWith(as.character(df$tissue), "T10B"),]$time = 10
df[startsWith(as.character(df$tissue), "T11B"),]$time = 11
df[startsWith(as.character(df$tissue), "T12B"),]$time = 12
df[startsWith(as.character(df$tissue), "T13B"),]$time = 13

df = df[!(is.na(df$value)),]
#names(df)[3]= "N0"
store_df <- df

df=df[df$time!=0,]
names(df)[2] = "gene"

q <- ggplot(df, aes(x = time, y = value, color = gene)) + geom_point(position = position_dodge(width=0.5)) + theme_bw()
q <- q + scale_x_continuous(breaks = seq(1,13), labels=c(12, 16, 20, 0, 4, 8, 12, 16, 20, 0, 4, 8, 12))
q

test = data.frame(time=rep(seq(1,13), each =5))

temp = group_by(df, time, gene) %>% summarise(means = mean(value)) 
test$means = temp$means
test$gene = temp$gene

temp = group_by(df, time, gene) %>% summarise(sem = std.error(value))
test$sem = temp$sem

tmp = test %>% filter(gene == "Cry2")

p <- ggplot(tmp, aes(x=time, y=means*100000)) + geom_point() + geom_line() + theme_bw() + geom_errorbar(aes(ymin=(means-sem)*100000, ymax=(means+sem)*100000)) + scale_x_continuous(breaks = seq(1,13), labels=c(12, 16, 20, 0, 4, 8, 12, 16, 20, 0, 4, 8, 12))

p <- p + theme(text = element_text(size=15), axis.text.x=element_text(vjust=0.5)) + ylab("Relative Expression") + xlab("Circadian Time (ZT)") + ggtitle("Hes Expression (cycle, norm)")

p

df$circ_time = ((df$time + 2)  %% 6) * 4

b <- ggplot(df, aes(y=value, x=circ_time, color = gene)) + geom_point() + theme_bw() + scale_x_continuous(breaks=sort(unique(df$circ_time)), labels=c(0,4,8,12,16,20)) + geom_smooth()
b


lomb = lsp(tmp$means, type='period', from=5, to=7, ofac=30)
lomb

subset = df[df$gene == "Clk",]

perm_lomb_fun <- function(df,nperm = 2000, start=5, end=7) {
  lsp_fun <- function(df){
    temp = group_by(df, time) %>% summarise(means = mean(value))
    lomb = lsp(temp$means, type='period', from=start, to=end, ofac=30, plot=FALSE)
    return(lomb$peak)
  }
  
  obs_power = lsp_fun(df)
  
  mean_shuffle_fun = function(perm){
    shuff = df
    shuff$value = df[perm,]$value
    return(lsp_fun(shuff))
  }
  
  ## generate the required set of permutations
  pset <- shuffleSet(nrow(df), nset = nperm)
  
  ## iterate over the set of permutations applying meanDif
  
  D <- apply(pset, 1, mean_shuffle_fun)
  
  D <- c(obs_power, D)
  
  ## compute & return the p-value
  Ds <- sum(D >= D[1]) # how many >= to the observed diff?
  return(Ds / (nperm + 1)) # what proportion of perms is this (the pval)
}

a = perm_lomb_fun(subset)
a

Clk = a
test$gene <- as.character(test$gene)
test[test$gene=="Clk",]$gene <- "Clock"
test[test$gene=="Hes",]$gene <- "Helt"

ID = "Cry2"
tmp = test %>% filter(gene == ID)

p <- ggplot(tmp, aes(x=time, y=means*100000)) + geom_point() + geom_line() + theme_bw() + geom_errorbar(aes(ymin=(means-sem)*100000, ymax=(means+sem)*100000)) + scale_x_continuous(breaks = seq(1,13), labels=c(12, 16, 20, 0, 4, 8, 12, 16, 20, 0, 4, 8, 12))

p <- p + theme(text = element_text(size=15), axis.text.x=element_text(vjust=0.5)) + ylab("Relative Expression") + xlab("Circadian Time (ZT)") + 
  ggtitle(paste0(ID, " expression (cycle)"))

png(res=600, width=6, height=4, units="in", file=file.path("F:/circadian/qpcr/", paste0(ID, "_cycle.png") ))
p + theme(text=element_text(size=20, color="black"))
dev.off()


