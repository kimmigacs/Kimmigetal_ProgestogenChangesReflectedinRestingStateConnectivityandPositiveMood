# plots_hormonesRSAmanuscript.R

## info
# - E2, EE, P4 (progesteron), P (progestin), T

##################
##### set up #####
##################

rm(list=ls())

### set wd
#setwd("path_to_data_if_not_same_file_as_code")

### libraries
library(ggplot2)
library(ggpubr)
#library(ggbreak)
#library(ggprism)
library(glue)
#library(rstatix)
library(exactRankTests)

### helper functions
pstar <- function(pvals){
    psym <- ifelse(pvals < .05, "*", "ns")
    psym[pvals < .01] <- "**"
    psym[pvals < .001] <- "***"
    psym[is.na(pvals)] <- ""
    return(psym)
}

pw.wilcox.exact <- function(x, grp, correction="bonferroni"){
    grp <- factor(grp)
    combs <- t(combn(1:nlevels(grp), 2))
    ps <- numeric(nrow(combs))
    for(i in 1:nrow(combs)){
        ps[i] <- wilcox.exact(x[grp == combs[i, 1]],
                              x[grp == combs[i, 2]],
                              exact=TRUE)$p.value*nrow(combs)
    }
    ps[ps > 1] <- 1
    ret <- as.data.frame(cbind(combs, p.adj=ps))
    names(ret)[1:2] <- c("group1", "group2")
    ret$p.adj.signif <- pstar(ret$p.adj)
    return(ret)
}

single.means <- function(x, grp, funs){
    ret <- data.frame(group1=1:length(funs),
                      group2=1:length(funs),
                      p.adj=NA)
    grp <- factor(grp)
    for(i in 1:length(funs)){
        if(funs[i] == "") next
        ret[i, "p.adj"] <- do.call(funs[i], 
                                   list(x[grp == levels(grp)[i]]))$p.value*4
    }
    ret$p.adj[ret$p.adj > 1] <- 1
    ret$p.adj.signif <- pstar(ret$p.adj)
    return(ret)
}

########################
##### prepare data #####
########################

# read in data
mydata_delta <- read.csv("Barplot_SourceData_anonymized.csv", header=TRUE, stringsAsFactors=TRUE)

# setup variables for graphs
Group_delta <-c(mydata_delta$group_noONC, mydata_delta$group_noONC, 
                mydata_delta$group_noONC, mydata_delta$group_noONC, 
                mydata_delta$group_noONC)
Group_delta <- factor(Group_delta, levels = c(1, 2, 3, 4), 
                      labels = c("fNC", "OC", "stopOC", "startOC"))

### adjust units of hormones
# piko to nano
mydata_delta$Testosterone_delta <- 
    mydata_delta$Testosterone_delta * 3.467 * 0.001 
mydata_delta$Progesterone_delta <- 
    mydata_delta$Progesterone_delta * 3.18 * 0.001 
mydata_delta$Gestagen_delta <- 
    mydata_delta$Gestagen_delta * 3.18 * 0.001
# stay with piko
mydata_delta$Estradiol_delta <- mydata_delta$Estradiol_delta * 3.671 
mydata_delta$EE_delta <- mydata_delta$EE_delta * 3.671

# get hormone data together for graph
hormonesdelta_data <-c(mydata_delta$Testosterone_delta, 
                       mydata_delta$Estradiol_delta, 
                       mydata_delta$Progesterone_delta, 
                       mydata_delta$EE_delta, 
                       mydata_delta$Gestagen_delta)

hormones_delta <- rep(c("Testosterone", "Estradiol", "Progesterone", 
                        "Ethinylestradiol", "Progestin") , each = 88)
hormones_delta <- factor(hormones_delta, 
                         levels = c("Testosterone", "Estradiol", 
                                    "Progesterone", "Ethinylestradiol", 
                                    "Progestin"),
                   labels = c("Testosterone (T1-T0, nmol/l)", 
                              "Estradiol (T1-T0, pmol/l)", 
                              "Progesterone (T1-T0, nmol/l)", 
                              "Ethinylestradiol (T1-T0, pmol/l)",
                              "Progestin (T1-T0, nmol/l)"))

mydatadelta_plots <- data.frame(Group_delta, hormonesdelta_data, hormones_delta)

##################################
##### bar charts per hormone #####
##################################
theme_set(ggpubr::theme_pubr())

## prepare plot parameters
# bar width
barw <- 0.9
# y-labels
labs <- c("fNC", "OC", "dOC", "sOC")
# pointsize axis labs + strip.text
pts.x <- 12
pts.y <- 10
pts.st <- 14
# jitter extent
jitt <- 0.7
# significance label size, vertical adjustment
ssize <- 5
sva <- 0.6
# error bar size, width
ebs <- 0.55
ebw <- 0.4
# plot margins in cm (top, right, bottom, left); 0.03 to avoid cut offs when
# zooming out
margs <- c(0, 0.03, 0, 0)

# -------------- #
# --- plot 1 --- #
# -------------- #
aa <- subset(mydatadelta_plots, hormones_delta=="Testosterone (T1-T0, nmol/l)")
aa$Group_deltaN <- as.numeric(aa$Group_delta)

## ANOVA + TUKEY #
#aov.res <- aov(hormonesdelta_data ~ as.factor(Group_deltaN), aa)
#tuk <- aov.res %>% tukey_hsd() # 95% familywise alpha
#tuk$p.adj <- sprintf("%.3f", round(tuk$p.adj, 3)) 
#tuk$p.adj.signif[tuk$p.adj.signif == "****"] <- "***"

## PAIRWISE WILCOXON #
#res <- pairwise.wilcox.test(aa$hormonesdelta_data, aa$Group_deltaN, 
#                            p.adjust.method="bonferroni")
#d <- as.data.frame(cbind(t(combn("1":"4", 2)), 
#                         p.adj=as.vector(na.omit(c(res$p.value)))))
#names(d)[1:2] <- c("group1", "group2")
#d$p.adj.signif <- pstar(d$p.adj)

# PAIRWISE EXACT WILCOXON #
res <- pw.wilcox.exact(aa$hormonesdelta_data, aa$Group_deltaN)

# SINGLE MEAN #
# t-test - fNC, OC-, OC+
# wilcox - OC
# none - 
# group order: fNC, OC, OC-, OC+
funs <- c("t.test", "wilcox.exact", "t.test", "t.test")
res.single <- single.means(aa$hormonesdelta_data, aa$Group_deltaN, funs)

# PLOT #
testosteronedelta_graph <- 

  ggbarplot(aa,
            x = "Group_deltaN", y = "hormonesdelta_data", width = barw,
            color = "Group_deltaN", fill = "Group_deltaN",
            position = position_dodge(0.25),
            add = c("mean_se"), 
            add.params = list(width = ebw, size = ebs)) +
  geom_jitter(aes(Group_deltaN, hormonesdelta_data, fill = Group_deltaN), 
              shape = 21, size = 2, alpha = 0.5,
              position = position_jitterdodge(jitter.height = 0, 
                                              jitter.width = jitt)) +
  stat_pvalue_manual(res[res$p.adj.signif != "ns", ], label = "p.adj.signif",
                     y.position = 1, tip.length = 0,
                     bracket.shorten = 0.1, size = ssize, vjust = sva) +
  add_pvalue(res.single[res.single$p.adj.signif != "ns", ], 
             label = "p.adj.signif", label.size = ssize,
             y.position = -1.40, vjust = 2.5 * sva,
             tip.length = 0, remove.bracket = TRUE) +
  facet_grid(. ~ glue('"T ("*Delta*", nmol/l)"'), 
             labeller = label_parsed) +
  scale_color_manual(values = c("black", "black", "black", "black")) +
  scale_x_discrete(labels = labs) +
  ylim(-1.5, 1) +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(size = pts.x),
        axis.title.y = element_blank(), 
        axis.text.y = element_text(size = pts.y),
        axis.line = element_line(linewidth = 0),
        strip.background = element_rect(linewidth = 0),
        legend.position = 'none',
        strip.text.x = element_text(size = pts.st),
        plot.margin = unit(margs, "cm")) 

# -------------- #
# --- plot 2 --- #
# -------------- #
aa <- subset(mydatadelta_plots, hormones_delta=="Estradiol (T1-T0, pmol/l)")
aa$Group_deltaN <- as.numeric(aa$Group_delta)

# PAIRWISE EXACT WILCOXON #
res <- pw.wilcox.exact(aa$hormonesdelta_data, aa$Group_deltaN)

# SINGLE MEAN #
# t-test - 
# wilcox - fNC, OC, OC-, OC+
# none - 
# group order: fNC, OC, OC-, OC+
funs <- rep("wilcox.exact", 4)
res.single <- single.means(aa$hormonesdelta_data, aa$Group_deltaN, funs)

# PLOT #
estradioldelta_graph <- 

    ggbarplot(aa,
              x = "Group_deltaN", y = "hormonesdelta_data", width=barw,
              color = "Group_deltaN", fill="Group_deltaN",
              position = position_dodge(0.25),
              add = c("mean_se"), 
              add.params = list(width=ebw, size = ebs)) +
    geom_jitter(aes(Group_deltaN, hormonesdelta_data, fill = Group_deltaN), 
                shape = 21, size=2, alpha=0.5,
                position = position_jitterdodge(jitter.height = 0, 
                                                jitter.width = jitt)) +
    stat_pvalue_manual(res[res$p.adj.signif != "ns", ], label="p.adj.signif",
                       y.position=c(415, 505, 370, 460, 370), vjust=sva,
                       tip.length=0, bracket.shorten=0.1, size=ssize) +
    add_pvalue(res.single[res.single$p.adj.signif != "ns", ], 
               label="p.adj.signif", label.size=ssize,
               y.position=rep(-340, 2), vjust=2.5*sva,
               tip.length=0, remove.bracket=TRUE) +
    facet_grid(. ~ glue('"E2 ("*Delta*", pmol/l)"'), 
               labeller=label_parsed) +
    scale_color_manual(values = c("black", "black", "black", "black")) +
    scale_x_discrete(labels=labs) +
    ylim(-370, 510) +
    theme(axis.title.x = element_blank(), 
          axis.text.x = element_text(#face="bold", 
                                     size = pts.x),
          axis.title.y = element_blank(), 
          axis.text.y = element_text(size = pts.y),
          axis.line = element_line(linewidth=0),
          strip.background = element_rect(linewidth=0),
          legend.position = 'none',
          strip.text.x = element_text(size = pts.st),
          plot.margin = unit(margs, "cm")) 

# -------------- #
# --- plot 3 --- #
# -------------- #
aa <- subset(mydatadelta_plots, hormones_delta=="Progesterone (T1-T0, nmol/l)")
aa$Group_deltaN <- as.numeric(aa$Group_delta)

# PAIRWISE EXACT WILCOXON #
res <- pw.wilcox.exact(aa$hormonesdelta_data, aa$Group_deltaN)

# SINGLE MEAN #
# t-test - fNC
# wilcox - OC, OC-, OC+
# none - fNC (no variance)
# group order: fNC, OC, OC-, OC+
funs <- c("t.test", "wilcox.exact", "wilcox.exact", "wilcox.exact")
res.single <- single.means(aa$hormonesdelta_data, aa$Group_deltaN, funs)

# PLOT #
progesteronedelta_graph <- 

    ggbarplot(aa,
              x = "Group_deltaN", y = "hormonesdelta_data", width=barw,
              color = "Group_deltaN", fill="Group_deltaN",
              position = position_dodge(0.25),
              add = c("mean_se"), 
              add.params = list(width=ebw, size = ebs)) +
    geom_jitter(aes(Group_deltaN, hormonesdelta_data, fill = Group_deltaN), 
                shape = 21, size=2, alpha=0.5,
                position = position_jitterdodge(jitter.height = 0, 
                                                jitter.width = jitt)) +
    stat_pvalue_manual(res[res$p.adj.signif != "ns", ], label="p.adj.signif",
                       y.position=c(15.25, 13.5, 17, 13.5), tip.length=0, 
                       bracket.short=0.1, size=ssize, vjust=sva) +
    add_pvalue(res.single[res.single$p.adj.signif != "ns", ], 
               label="p.adj.signif", label.size=ssize,
               y.position=rep(-18, 2), vjust=2.5*sva,
               tip.length=0, remove.bracket=TRUE) +
    facet_grid(. ~ glue('"P4 ("*Delta*", nmol/l)"'), 
               labeller=label_parsed) +
    scale_color_manual(values = c("black", "black", "black", "black")) +
    scale_x_discrete(labels=labs) +
    ylim(-19.1, 17.2) +
    theme(axis.title.x = element_blank(), 
          axis.text.x = element_text(#face="bold", 
                                     size = pts.x),
          axis.title.y = element_blank(), 
          axis.text.y = element_text(size = pts.y),
          axis.line = element_line(linewidth=0),
          strip.background = element_rect(linewidth=0),
          legend.position = 'none',
          strip.text.x = element_text(size = pts.st),
          plot.margin = unit(margs, "cm")) 

# -------------- #
# --- plot 4 --- #
# -------------- #
aa <- subset(mydatadelta_plots,
             hormones_delta=="Ethinylestradiol (T1-T0, pmol/l)")
aa$Group_deltaN <- as.numeric(aa$Group_delta)

# PAIRWISE COMPARISONS #
# wilcox - all except with group fNC (no variance)
#   => set to "ns" s.t. it wont be printed in plot
res <- pw.wilcox.exact(aa$hormonesdelta_data, aa$Group_deltaN)
res[1:3, "p.adj.signif"] <- "ns"

# SINGLE MEAN #
# t-test - 
# wilcox - OC, OC-, OC+
# none - fNC (no variance)
# group order: fNC, OC, OC-, OC+
funs <- c("", "wilcox.exact", "wilcox.exact", "wilcox.exact")
res.single <- single.means(aa$hormonesdelta_data, aa$Group_deltaN, funs)

# PLOT #
EEdelta_graph <- 

    ggbarplot(aa,
              x = "Group_deltaN", y = "hormonesdelta_data", width=barw,
              color = "Group_deltaN", fill="Group_deltaN",
              position = position_dodge(0.25),
              add = c("mean_se"), 
              add.params = list(width=ebw, size = ebs)) +
    geom_jitter(aes(Group_deltaN, hormonesdelta_data, fill = Group_deltaN), 
                shape = 21, size=2, alpha=0.5,
                position = position_jitterdodge(jitter.height = 0, 
                                                jitter.width = jitt)) +
    stat_pvalue_manual(res[res$p.adj.signif != "ns", ], label="p.adj.signif",
                       y.position=c(400, 445, 400) - 5, vjust=sva,
                       tip.length=0, bracket.shorten=0.1, size=ssize) +
    add_pvalue(res.single[res.single$p.adj.signif != "ns", ], 
               label="p.adj.signif", label.size=ssize,
               y.position=rep(-340, 3), vjust=2.5*sva,
               tip.length=0, remove.bracket=TRUE) +
    facet_grid(. ~ glue('"EE ("*Delta*", pmol/l)"'), 
               labeller=label_parsed) +
    scale_color_manual(values = c("black", "black", "black", "black")) +
    scale_x_discrete(labels=labs) +
    ylim(-370, 447) +
    theme(axis.title.x = element_blank(), 
          axis.text.x = element_text(#face="bold", 
                                     size = pts.x),
          axis.title.y = element_blank(), 
          axis.text.y = element_text(size = pts.y),
          axis.line = element_line(linewidth=0),
          strip.background = element_rect(linewidth=0),
          legend.position = 'none',
          strip.text.x = element_text(size = pts.st),
          plot.margin = unit(margs, "cm")) 

# -------------- # 
# --- plot 5 --- # 
# -------------- #
aa <- subset(mydatadelta_plots, hormones_delta=="Progestin (T1-T0, nmol/l)")
aa$Group_deltaN <- as.numeric(aa$Group_delta)

# PAIRWISE COMPARISONS #
# wilcox - all except with group fNC (no variance)
#   => set to "ns" s.t. it wont be printed in plot
res <- pw.wilcox.exact(aa$hormonesdelta_data, aa$Group_deltaN)
res[1:3, "p.adj.signif"] <- "ns"

# SINGLE MEAN #
# t-test - 
# wilcox - OC, OC-, OC+
# group order: fNC, OC, OC-, OC+
funs <- c("", "wilcox.exact", "wilcox.exact", "wilcox.exact")
res.single <- single.means(aa$hormonesdelta_data, aa$Group_deltaN, funs)

# PLOT #
Gestagendelta_graph <- 

    ggbarplot(aa,
              x = "Group_deltaN", y = "hormonesdelta_data", width=barw,
              color = "Group_deltaN", fill="Group_deltaN",
              position = position_dodge(0.25),
              add = c("mean_se"), 
              add.params = list(width=ebw, size = ebs)) +
    geom_jitter(aes(Group_deltaN, hormonesdelta_data, fill = Group_deltaN), 
                shape = 21, size=2, alpha=0.5,
                position = position_jitterdodge(jitter.height = 0, 
                                                jitter.width = jitt)) +
    stat_pvalue_manual(res[res$p.adj.signif != "ns", ], label="p.adj.signif",
                       y.position=c(130, 145, 130) + 10, vjust=sva,
                       tip.length=0, bracket.shorten=0.1, size=ssize) +
    add_pvalue(res.single[res.single$p.adj.signif != "ns", ], 
               label="p.adj.signif", label.size=ssize,
               y.position=rep(-102, 3), vjust=2.5*sva,
               tip.length=0, remove.bracket=TRUE) +
    facet_grid(. ~ glue('"P ("*Delta*", nmol/l)"'), 
               labeller=label_parsed) +
    scale_color_manual(values = c("black", "black", "black", "black")) +
    scale_x_discrete(labels=labs) +
    ylim(-112, 157.5) +
    theme(axis.title.x = element_blank(), 
          axis.text.x = element_text(#face="bold", 
                                     size = pts.x),
          axis.title.y = element_blank(), 
          axis.text.y = element_text(size = pts.y),
          axis.line = element_line(linewidth=0),
          strip.background = element_rect(linewidth=0),
          legend.position = 'none',
          strip.text.x = element_text(size = pts.st),
          plot.margin = unit(margs, "cm")) 

# --- total plot -- #
figure_hormonesdelta <- ggarrange(estradioldelta_graph, NULL, 
                                  EEdelta_graph, NULL, 
                                  progesteronedelta_graph, NULL,
                                  Gestagendelta_graph, NULL, 
                                  testosteronedelta_graph, 
                                  nrow = 1, 
                                  widths = c(1, 0, 1, 0, 1, 0, 1, 0, 1))

## save plot
hi <- 3
png("HormonesGroupDiff.png", units="in", width=10, height=3, res=300)
#tiff("HormonesGroupDiff.tiff", units="in", width=10, height=3, res=300)
#pdf("5plots_Hormones.pdf", height=hi, width=10)

print(figure_hormonesdelta)

dev.off()

#source("plots_hormonesRSAmanuscript_PW.R")

