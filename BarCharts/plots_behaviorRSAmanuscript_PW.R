# plots_behaviorRSAmanuscript.R

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

pw.t.test <- function(x, grp, var.eq, correction="bonferroni"){
    grp <- factor(grp)
    combs <- t(combn(1:nlevels(grp), 2))
    ps <- numeric(nrow(combs))
    for(i in 1:nrow(combs)){
        ps[i] <- t.test(x[grp == combs[i, 1]],
                        x[grp == combs[i, 2]],
                        var.equal=var.eq)$p.value*nrow(combs)
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

# transform ER RT from ms to s
mydata_delta$EmoRec_TotRT_delta <- mydata_delta$EmoRec_TotRT_delta/1000

# setup variables for graphs
Group_delta <-c(mydata_delta$group_noONC, mydata_delta$group_noONC, 
                mydata_delta$group_noONC, mydata_delta$group_noONC, 
                mydata_delta$group_noONC)
Group_delta <- factor(Group_delta, levels = c(1, 2, 3, 4), 
                      labels = c("fNC", "OC", "stopOC", "startOC"))

# get behavior data together for graph
behaviordelta_data <-c(mydata_delta$EmoRec_TotAccuracy_delta, 
                       mydata_delta$EmoRec_TotRT_delta,
                       mydata_delta$PANAS_pos_delta, 
                       mydata_delta$PANAS_neg_delta,
                       mydata_delta$BDI_delta)

behavior_delta <- rep(c("ER_tAcc", "ER_tRT", "PANAS_p", 
                        "PANAS_n", "BDI") , each = 88)
behavior_delta <- factor(behavior_delta, 
                         levels = c("ER_tAcc", "ER_tRT", "PANAS_p", 
                        "PANAS_n", "BDI"))

mydatadelta_plots <- data.frame(Group_delta, behaviordelta_data, behavior_delta)

##################################
##### bar charts per behavior #####
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
aa <- subset(mydatadelta_plots, behavior_delta=="ER_tAcc")
aa$Group_deltaN <- as.numeric(aa$Group_delta)

## ANOVA + TUKEY #
#aov.res <- aov(behaviordelta_data ~ as.factor(Group_deltaN), aa)
#tuk <- aov.res %>% tukey_hsd() # 95% familywise alpha
#tuk$p.adj <- sprintf("%.3f", round(tuk$p.adj, 3)) 
#tuk$p.adj.signif[tuk$p.adj.signif == "****"] <- "***"

## PAIRWISE WILCOXON #
#res <- pairwise.wilcox.test(aa$behaviordelta_data, aa$Group_deltaN, 
#                            p.adjust.method="bonferroni")
#d <- as.data.frame(cbind(t(combn("1":"4", 2)), 
#                         p.adj=as.vector(na.omit(c(res$p.value)))))
#names(d)[1:2] <- c("group1", "group2")
#d$p.adj.signif <- pstar(d$p.adj)

# PAIRWISE COMPARISON #
# t-test - all
res <- pw.t.test(aa$behaviordelta_data, aa$Group_deltaN, var.eq=TRUE)

# SINGLE MEAN #
# t-test - fNC, OC, OC-
# wilcox - OC+
# none - 
# group order: fNC, OC, OC-, OC+
funs <- c("t.test", "t.test", "t.test", "wilcox.exact")
res.single <- single.means(aa$behaviordelta_data, aa$Group_deltaN, funs)

# PLOT #
ERtA_graph <- 

    ggbarplot(aa,
              x = "Group_deltaN", y = "behaviordelta_data", width=barw,
              color = "Group_deltaN", fill="Group_deltaN",
              position = position_dodge(0.25),
              add = c("mean_se"), 
              add.params = list(width=ebw, size = ebs)) +
    geom_jitter(aes(Group_deltaN, behaviordelta_data, fill = Group_deltaN), 
                shape = 21, size=2, alpha=0.5,
                position = position_jitterdodge(jitter.height = 0, 
                                                jitter.width = jitt)) +
#    stat_pvalue_manual(res[res$p.adj.signif != "ns", ], label="p.adj.signif",
#                       #y.position=c(0.95, 0.8, 0.8) + 0.2, tip.length=0, TUKEY
#                       y.position=1, tip.length=0,
#                       bracket.shorten=0.1, size=ssize, vjust=sva) +
#    add_pvalue(res.single[res.single$p.adj.signif != "ns", ], 
#               label="p.adj.signif", label.size=ssize,
#               y.position=-1.40, vjust=2.5*sva,
#               tip.length=0, remove.bracket=TRUE) +
    facet_grid(. ~ glue('"Total ER-A ("*Delta*", %)"'), 
               labeller=label_parsed) +
    scale_color_manual(values = c("black", "black", "black", "black")) +
    scale_x_discrete(labels=labs) +
#    ylim(, ) +
    theme(axis.title.x = element_blank(), 
          axis.text.x = element_text(#face="bold", 
                                     size = pts.x),
          axis.title.y = element_blank(), 
          axis.text.y = element_text(size = pts.y),
          legend.position = 'none',
          strip.text.x = element_text(size = pts.st),
          axis.line = element_line(linewidth=0),
          strip.background = element_rect(linewidth=0),
          plot.margin = unit(margs, "cm"))

# -------------- #
# --- plot 2 --- #
# -------------- #
aa <- subset(mydatadelta_plots, behavior_delta=="ER_tRT")
aa$Group_deltaN <- as.numeric(aa$Group_delta)

# PAIRWISE COMPARISON #
# wilcox - all
res <- pw.wilcox.exact(aa$behaviordelta_data, aa$Group_deltaN)

# SINGLE MEAN #
# t-test - 
# wilcox - fNC, OC, OC-, OC+
# none - 
# group order: fNC, OC, OC-, OC+
funs <- rep("wilcox.exact", 4)
res.single <- single.means(aa$behaviordelta_data, aa$Group_deltaN, funs)

# PLOT #
ERtRT_graph  <- 

    ggbarplot(aa,
              x = "Group_deltaN", y = "behaviordelta_data", width=barw,
              color = "Group_deltaN", fill="Group_deltaN",
              position = position_dodge(0.25),
              add = c("mean_se"), 
              add.params = list(width=ebw, size = ebs)) +
    geom_jitter(aes(Group_deltaN, behaviordelta_data, fill = Group_deltaN), 
                shape = 21, size=2, alpha=0.5,
                position = position_jitterdodge(jitter.height = 0, 
                                                jitter.width = jitt)) +
#    stat_pvalue_manual(res[res$p.adj.signif != "ns", ], label="p.adj.signif",
#                       y.position=c(415, 505, 370, 460, 370), vjust=sva,
#                       tip.length=0, bracket.shorten=0.1, size=ssize) +
    add_pvalue(res.single[res.single$p.adj.signif != "ns", ], 
               label="p.adj.signif", label.size=ssize,
               y.position=rep(-2.5, 3), vjust=2.5*sva,
               tip.length=0, remove.bracket=TRUE) +
    facet_grid(. ~ glue('"Total ER-T ("*Delta*", s)"'), 
               labeller=label_parsed) +
    scale_color_manual(values = c("black", "black", "black", "black")) +
    scale_x_discrete(labels=labs) +
    ylim(-2.6, 0.6) +
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
aa <- subset(mydatadelta_plots, behavior_delta=="PANAS_p")
aa$Group_deltaN <- as.numeric(aa$Group_delta)

# PAIRWISE COMPARISON #
# wilcox - all
res <- pw.wilcox.exact(aa$behaviordelta_data, aa$Group_deltaN)

# SINGLE MEAN #
# t-test - OC
# wilcox - fNC, OC-, OC+
# none - 
# group order: fNC, OC, OC-, OC+
funs <- c("wilcox.exact", "t.test", "wilcox.exact", "wilcox.exact")
res.single <- single.means(aa$behaviordelta_data, aa$Group_deltaN, funs)

# PLOT #
PANASp_graph <- 

    ggbarplot(aa,
              x = "Group_deltaN", y = "behaviordelta_data", width=barw,
              color = "Group_deltaN", fill="Group_deltaN",
              position = position_dodge(0.25),
              add = c("mean_se"), 
              add.params = list(width=ebw, size = ebs)) +
    geom_jitter(aes(Group_deltaN, behaviordelta_data, fill = Group_deltaN), 
                shape = 21, size=2, alpha=0.5,
                position = position_jitterdodge(jitter.height = 0, 
                                                jitter.width = jitt)) +
#    stat_pvalue_manual(res[res$p.adj.signif != "ns", ], label="p.adj.signif",
#                       y.position=c(15.25, 13.5, 17, 13.5), tip.length=0, 
#                       bracket.short=0.1, size=ssize, vjust=sva) +
#    add_pvalue(res.single[res.single$p.adj.signif != "ns", ], 
#               label="p.adj.signif", label.size=ssize,
#               y.position=rep(-18, 2), vjust=2.5*sva,
#               tip.length=0, remove.bracket=TRUE) +
    facet_grid(. ~ glue('"Mood+ ("*Delta*")"'), 
               labeller=label_parsed) +
    scale_color_manual(values = c("black", "black", "black", "black")) +
    scale_x_discrete(labels=labs) +
#    ylim(, ) +
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
             behavior_delta=="PANAS_n")
aa$Group_deltaN <- as.numeric(aa$Group_delta)

# PAIRWISE COMPARISONS #
# wilcox - all
res <- pw.wilcox.exact(aa$behaviordelta_data, aa$Group_deltaN)

# SINGLE MEAN #
# t-test - fNC
# wilcox - OC, OC-, OC+
# none - 
# group order: fNC, OC, OC-, OC+
funs <- c("t.test", "wilcox.exact", "wilcox.exact", "wilcox.exact")
res.single <- single.means(aa$behaviordelta_data, aa$Group_deltaN, funs)

# PLOT #
PANASn_graph <- 

    ggbarplot(aa,
              x = "Group_deltaN", y = "behaviordelta_data", width=barw,
              color = "Group_deltaN", fill="Group_deltaN",
              position = position_dodge(0.25),
              add = c("mean_se"), 
              add.params = list(width=ebw, size = ebs)) +
    geom_jitter(aes(Group_deltaN, behaviordelta_data, fill = Group_deltaN), 
                shape = 21, size=2, alpha=0.5,
                position = position_jitterdodge(jitter.height = 0, 
                                                jitter.width = jitt)) +
#    stat_pvalue_manual(res[res$p.adj.signif != "ns", ], label="p.adj.signif",
#                       y.position=c(400, 445, 400) - 5, vjust=sva,
#                       tip.length=0, bracket.shorten=0.1, size=ssize) +
#    add_pvalue(res.single[res.single$p.adj.signif != "ns", ], 
#               label="p.adj.signif", label.size=ssize,
#               y.position=rep(-340, 3), vjust=2.5*sva,
#               tip.length=0, remove.bracket=TRUE) +
    facet_grid(. ~ glue('"Mood- ("*Delta*")"'), 
               labeller=label_parsed) +
    scale_color_manual(values = c("black", "black", "black", "black")) +
    scale_x_discrete(labels=labs) +
#    ylim(, ) +
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
aa <- subset(mydatadelta_plots, behavior_delta=="BDI")
aa$Group_deltaN <- as.numeric(aa$Group_delta)

# PAIRWISE COMPARISONS #
# wilcox - all
res <- pw.wilcox.exact(aa$behaviordelta_data, aa$Group_deltaN)

# SINGLE MEAN #
# t-test - fNC, OC
# wilcox - OC-, OC+
# group order: fNC, OC, OC-, OC+
funs <- c("t.test", "t.test", "wilcox.exact", "wilcox.exact")
res.single <- single.means(aa$behaviordelta_data, aa$Group_deltaN, funs)

# PLOT #
BDI_graph <- 

    ggbarplot(aa,
              x = "Group_deltaN", y = "behaviordelta_data", width=barw,
              color = "Group_deltaN", fill="Group_deltaN",
              position = position_dodge(0.25),
              add = c("mean_se"), 
              add.params = list(width=ebw), size = ebs) +
    geom_jitter(aes(Group_deltaN, behaviordelta_data, fill = Group_deltaN), 
                shape = 21, size=2, alpha=0.5,
                position = position_jitterdodge(jitter.height = 0, 
                                                jitter.width = jitt)) +
    stat_pvalue_manual(res[res$p.adj.signif != "ns", ], label="p.adj.signif",
                       y.position=c(28, 25, 25), vjust=sva,
                       tip.length=0, bracket.shorten=0.1, size=ssize) +
    add_pvalue(res.single[res.single$p.adj.signif != "ns", ], 
               label="p.adj.signif", label.size=ssize,
               y.position=-26.5, vjust=2.5*sva,
               tip.length=0, remove.bracket=TRUE) +
    facet_grid(. ~ glue('"Depression ("*Delta*")"'), 
               labeller=label_parsed) +
    scale_color_manual(values = c("black", "black", "black", "black")) +
    scale_x_discrete(labels=labs) +
    ylim(-28.5, 28.5) +
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
figure_behaviordelta <- ggarrange(PANASp_graph, NULL,
                                  PANASn_graph, NULL, 
                                  BDI_graph, NULL,
                                  ERtA_graph, NULL, 
                                  ERtRT_graph, 
                                  nrow = 1, 
                                  widths = c(1, 0, 1, 0, 1, 0, 1, 0, 1))

## save plot
hi <- 3
png("BehaviourGroupDiff.png", units="in", width=10, height=3, res=300)
#tiff("BehaviourGroupDiff.tiff", units="in", width=10, height=3, res=300)
#pdf("5plots_behav.pdf", height=hi, width=10)
#ggsave(filename = "C:\Users\UKPP\Documents\HormonesRestingStateRSA\Manuskript_preprintRSA\BehaviourGroupDiff.jpg",
#       plot=figure_behaviordelta,
#       width = 10, height = hi, dpi=300, units = "cm")
print(figure_behaviordelta)

dev.off()

#source("plots_behaviorRSAmanuscript_PW.R")

