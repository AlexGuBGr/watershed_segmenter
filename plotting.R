library(ggplot2)   #"v3.4.3"
library(ggsignif)  #"v0.6.4"


p_n_adjp <- function(lst, dat, labs, test, alt) {
  
  ps <- rep(NA, length(lst))
  
  for (i in 1:length(lst)) {
    tmpt <- test(x = dat[ labs == lst[[i]][[1]] ], y = dat[ labs == lst[[i]][[2]] ], alternative=alt)
    print(tmpt)
    ps[[i]] <- tmpt$p.value
  }
  ps
}

permtest <- function(st32vis, conds) {
  set.seed(1234)
  st32vis <- st32vis[st32vis$lables %in% conds,]
  st32vis$lables <- factor(st32vis$lables, levels = conds, ordered = TRUE)
  
  perm_test <- independence_test(Cell_size ~ lables, data = st32vis, distribution = approximate(nresample = 10000))
  pvalue(perm_test)
}

#pp <- permtest(st32vis, c("CT0", "60"))


boxplotter <- function(strr, st32vis, logg, fact, test, alt, colr, ylab, name) {
  
  writexl::write_xlsx(rbind( data.frame( table(st32vis$lables)),  data.frame( Var1 = "total", Freq = dim(st32vis)[[1]] ) ),
                      path = paste0(name, "_count.xlsx") )
  
  
  rgb_colors <- c(
    rgb(128	,128,	128, maxColorValue = 255),    
    rgb(144,	191,	249, maxColorValue = 255),    
    rgb(5,	190,	120, maxColorValue = 255),
    rgb(255,	192,	128, maxColorValue = 255),
    rgb(249,	100,	149, maxColorValue = 255)
  )
  
  
  if (strr == "morning") {
    lvl <- c("CT0", "36", "60","84", "108")
    comps <- list(c("CT0", "36"),c("CT0", "60"),c("CT0", "84"),c("CT0", "108"))
  } else if (strr == "evening") {
    lvl <- c("CT12", "48", "72", "96", "120" )
    comps <- list( c("CT12", "48"), c("CT12","72"), c("CT12","96"), c("CT12","120") )
  }
  
  if (logg == T ) {
    st32vis$Cell_size <- log10(st32vis$Cell_size)
  }
  
  
  ns <- unique(st32vis$lables)
  df <- data.frame(t(sapply(ns, FUN = function(x) { tmp <- st32vis$Cell_size[st32vis$lables == x];  c( mean( tmp), median( tmp) )  } ) ))
  print(df)
  colnames(df) <- c("mean (log10)", "median (log10")
  rownames(df) <- ns
  df[["group"]] <- ns
  print(df)
  writexl::write_xlsx(df, path = paste0(name, "_mean_median.xlsx"))
  
  
  if (class(test) != "function") {
    
    ps <- rep(NA, length(comps))
    
    for (i in 1:length(comps)) {
      ps[[i]] <- permtest(st32vis, comps[[i]])
    }
    
  } else {
    ps <- p_n_adjp(comps, 
                   st32vis$Cell_size,
                   st32vis$lables,
                   test, alt
    )
  }
  
  pad <- p.adjust(ps, method = "bonferroni")
  
  writexl::write_xlsx(data.frame(comps = sapply(comps, FUN = function(x){ paste(x, collapse = "_vs_") } ),
                                   p.val = ps, 
                                   padj = pad), path = paste0(name, ".xlsx") )
  
  anno <- sapply(pad, FUN = function(x) { if (x <= 0.0001) {return("***")
  } else if (x <= 0.001) {return("**")
  } else if (x <= 0.01) {return("*")
  } else {return("ns")}
  })
  
  siglvl <- c("***"=0.0001, "**"=0.001, "*"=0.01)
  
  
  ypox <- ( ( (1:length(lvl)) - 0.5 ) * fact) + max(st32vis$Cell_size)
  
  st32vis$lables <- factor(st32vis$lables,
                           levels = lvl, ordered = TRUE)
  
  if (colr == "fill"){ 
    
    pl <- ggplot(st32vis[grepl(strr,st32vis$time),], aes(x=lables, y=Cell_size, fill = lables )) + 
      geom_boxplot() +
      geom_signif(comparisons = comps, annotations = anno, color = "black",
                  #map_signif_level=siglvl, 
                  y_position= ypox, test = wilcox.test ) #+ stat_summary(fun=mean, geom="point", shape="_", size=9, color="#4C4C4C", fill="#4C4C4C") + stat_summary(fun=mean, geom="point", shape="|", size=4, color="#4C4C4C", fill="#4C4C4C")
    pl <- pl + scale_fill_manual(values = rgb_colors) + theme_classic() + labs(y = ylab, x = NULL) + ylim(0.7, 4)
    
  } else {
    
    pl <- ggplot(st32vis[grepl(strr,st32vis$time),], aes(x=lables, y=Cell_size, color = lables )) + 
      geom_boxplot() +
      geom_signif(comparisons = comps, annotations = anno, color = "black",
                  #map_signif_level=siglvl, 
                  y_position= ypox, test = wilcox.test ) #+ stat_summary(fun=mean, geom="point", shape="_", size=9, color="#4C4C4C", fill="#4C4C4C") + stat_summary(fun=mean, geom="point", shape="|", size=4, color="#4C4C4C", fill="#4C4C4C")
    pl <- pl + scale_color_manual(values = rgb_colors) + theme_classic() + labs(y = ylab, x = NULL) + ylim(0.7, 4)
  }
  
  pl
  
  ggsave(name, units = "cm", width = 10, height = 10, dpi=300)
  
}



boxplotter44 <- function(strr, st32vis, logg, fact, test, alt, colr, ylab, name) {
  
  
  writexl::write_xlsx(rbind( data.frame( table(st32vis$lables)),  data.frame( Var1 = "total", Freq = dim(st32vis)[[1]] ) ),
                                 path = paste0(name, "_count.xlsx") )
  
  
  
  rgb_colors <- c(
    "#000000",
    rgb(212	,212,	212, maxColorValue = 255),    
    rgb(255,	177,	100, maxColorValue = 255),    
    
    rgb(192,	96,	0, maxColorValue = 255),
    
    rgb(255,	141,	28, maxColorValue = 255),
    
    rgb(123,	62,	0, maxColorValue = 255),
    
    rgb(96,	96	,96, maxColorValue = 255),
    rgb(236,	118,	0, maxColorValue = 255),
    rgb(79,	39,	0, maxColorValue = 255)
    
  )
  
  
  lvl <- c("Group before running","sedentary mice", "Low-runners", "High-runners", 
           "3-day Inactive Low-runners", "3-day Inactive High-runners",
           "sedentary mice*", "3-week Low-runners", "3-week High-runners"
  )
  
  subcatvec <- c("-", "sed\n    ", "low\nChronic activity\nperiod", "high\n     ", 
                 
                 "low\n           3-day inactivity\n            period", "high", 
                 
                 "sed", "low\n3-week inactivity\nperiod", "high")
  
  cats <- list(
    `Chronic activity period` = c( "sedentary mice", "Low-runners", "High-runners"),
    `3-day inactivity period` = c("3-day Inactive Low-runners", "3-day Inactive High-runners"),
    `3-week inactivity period` = c("sedentary mice*", "3-week Low-runners", "3-week High-runners"),
    `before run` = c("Group before running")
  )
  
  subcats <- list(
    `sed` = c( "sedentary mice", "sedentary mice*"),
    `low` = c("Low-runners", "3-day Inactive Low-runners", "3-week Low-runners"),
    `high` = c("High-runners", "3-day Inactive High-runners", "3-week High-runners"),
    `-` = c("Group before running")
  )
  
  st32vis[["cat"]] <- ""
  for (i in names(cats)) {
    st32vis[["cat"]][ st32vis$lables %in% cats[[i]] ] <- i
  }
  
  st32vis[["subcat"]] <- ""
  for (i in names(subcats)) {
    st32vis[["subcat"]][ st32vis$lables %in% subcats[[i]] ] <- i
  }
  
  comps <- list( c("sedentary mice", "Low-runners"), 
                 c("sedentary mice", "High-runners"),
                 c("sedentary mice", "3-day Inactive Low-runners"),
                 c("sedentary mice",  "3-day Inactive High-runners"),
                 
                 c("sedentary mice*", "3-week Low-runners"),
                 c("sedentary mice*", "3-week High-runners"),
                 
                 c("Group before running", "sedentary mice"),
                 c("Group before running", "High-runners")
  )
  
  ordr <- c(1,2,3,4, 1,2, 5,6)
  
  if (logg == T ) {
    st32vis$Cell_size <- log10(st32vis$Cell_size)
  }
  
  
  ns <- unique(st32vis$lables)
  df <- data.frame(t(sapply(ns, FUN = function(x) { tmp <- st32vis$Cell_size[st32vis$lables == x];  c( mean( tmp), median( tmp) )  } ) ))
  print(df)
  colnames(df) <- c("mean (log10)", "median (log10")
  rownames(df) <- ns
  df[["group"]] <- ns
  print(df)
  writexl::write_xlsx(df, path = paste0(name, "_mean_median.xlsx"))
  
  
  ps <- p_n_adjp(comps, 
                 st32vis$Cell_size,
                 st32vis$lables,
                 test, alt
  )
  pad <- p.adjust(ps, method = "bonferroni")
  
  writexl::write_xlsx(data.frame(comps = sapply(comps, FUN = function(x){ paste(x, collapse = "_vs_") } ),
                                 p.val = ps, 
                                 padj = pad), path = paste0(name, ".xlsx") )
  
  anno <- sapply(pad, FUN = function(x) { if (x <= 0.0001) {return("***")
  } else if (x <= 0.001) {return("**")
  } else if (x <= 0.01) {return("*")
  } else {return("ns")}
  })
  
  siglvl <- c("***"=0.0001, "**"=0.001, "*"=0.01)
  
  ypos <- ( ( (ordr) - 0.5 ) * fact) + max(st32vis$Cell_size)
  
  st32vis$lables <- factor(st32vis$lables, levels = lvl, ordered = TRUE)
  
  st32vis$cat <- factor(st32vis$cat, levels = c("before run", "Chronic activity period", "3-day inactivity period", "3-week inactivity period"), ordered = TRUE) 
  
  
  pl <- ggplot(st32vis, aes(x=lables, y=Cell_size, fill = lables )) + 
    geom_boxplot() +
    geom_signif(comparisons = comps, annotations = anno, color = "black",
                y_position= ypos, test = wilcox.test ) #+ stat_summary(fun=mean, geom="point", shape="_", size=7, color="#4C4C4C", fill="#4C4C4C") + stat_summary(fun=mean, geom="point", shape="|", size=2.5, color="#4C4C4C", fill="#4C4C4C")
  
  pl <- pl + scale_fill_manual(values = rgb_colors) + theme_classic() + labs(y = ylab, x = NULL) +
    scale_x_discrete(labels= subcatvec) + ylim(0.7, 4.3)
  
  pl
  ggsave(name, units = "cm", width = 16, height = 10, dpi=300)
  
}


st32vis <- readxl::read_excel("st32vis_10x_cleaned.xlsx")
boxplotter("evening",st32vis, logg = T, fact = 0.2, test = wilcox.test, "two.sided", colr = "fill", ylab = "Cells size (log10)", name = "st32vis_evening.png") # 50 for non-log
boxplotter("morning",st32vis, logg = T, fact = 0.2, test = wilcox.test, "two.sided", colr = "", ylab = "Cells size (log10)", name = "st32vis_morning.png") # 50 for non-log
#boxplotter("evening",st32vis, logg = F, fact = 50, test = wilcox.test, "two.sided", colr = "fill", ylab = "Cells size", name = "st32vis_nl_evening.png") # 50 for non-log
#boxplotter("morning",st32vis, logg = F, fact = 50, test = wilcox.test, "two.sided", colr = "", ylab = "Cells size", name = "st32vis_nl_morning.png") # 50 for non-log
rm(st32vis)

st32sub <- readxl::read_excel("st32sub_10x_cleaned.xlsx")
boxplotter("evening",st32sub, logg = T, fact = 0.2, test = wilcox.test, "two.sided", colr = "fill", ylab = "Cells size (log10)", name = "st32sub_evening.png") # 50 for non-log
boxplotter("morning",st32sub, logg = T, fact = 0.2, test = wilcox.test, "two.sided", colr = "", ylab = "Cells size (log10)", name = "st32sub_morning.png") # 50 for non-log
#boxplotter("evening",st32sub, logg = F, fact = 50, test = wilcox.test, "two.sided", colr = "fill", ylab = "Cells size", name = "st32sub_nl_evening.png") # 50 for non-log
#boxplotter("morning",st32sub, logg = F, fact = 50, test = wilcox.test, "two.sided", colr = "", ylab = "Cells size", name = "st32sub_nl_morning.png") # 50 for non-log
rm(st32sub)


st44vis <- readxl::read_excel("st44vis_10x_cleaned.xlsx")
boxplotter44("", st44vis, logg = T, fact = 0.2, test = wilcox.test, "two.sided", colr = "fill", ylab = "Cells size (log10)", name = "st44vis.png")
rm(st44vis)

st44sub <- readxl::read_excel("st44sub_10x_cleaned.xlsx")
boxplotter44("", st44sub, logg = T, fact = 0.2, test = wilcox.test, "two.sided", colr = "fill", ylab = "Cells size (log10)", name = "st44sub.png")
rm(st44sub)
