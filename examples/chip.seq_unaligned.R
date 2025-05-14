library(animint2)
data(chip.seq)
ann.colors <- c("noDifference"='#f6f4bf',"difference"="#ff7d7d")
fp.fn.map <- c(false.positive="noDifference",false.negative="difference")
fp.fn.colors <- ann.colors[fp.fn.map]
names(fp.fn.colors) <- names(fp.fn.map)
total.vars <- 168
set.linetypes <- c(validation="dotted", train="solid")
## nonzero pair errors, for dots on the select sample plot.
pair.errors <-
rbind(data.frame(chip.seq$fp.pairs, error.type="false.positive"),
  data.frame(chip.seq$fn.pairs, error.type="false.negative")
)
nonzero.pairs <- subset(pair.errors, errors > 0)
nonzero.pairs$error.type <- factor(nonzero.pairs$error.type, 
                                  c("false.positive", "false.negative"))
## false positive/negative rates for the error/model complexity plot.
fp.fn <- rbind(
  data.frame(chip.seq$fp.set, error.type="false.positive"),
  data.frame(chip.seq$fn.set, error.type="false.negative")
)
fp.fn$error.type <- factor(fp.fn$error.type, c("false.positive", "false.negative"))
fp.fn.nonzero <- subset(fp.fn, errors > 0)

# Unaligned visualization
viz <- animint(
  title = "ChIP-seq Unaligned",
  source = "https://github.com/suhaani-agarwal/animint-examples/blob/master/examples/chip.seq_unaligned.R",
  chroms = ggplot() +
  theme_bw()+
    theme_animint(width=300, height=270, margin(80)) +
    geom_segment(aes(0, chr.int, xend=bases/1e6, yend=chr.int),
                 data=chip.seq$chroms, color="grey") +
    geom_text(aes(0, chr.int, label=paste0("chr", chr)),
              data=chip.seq$set.info, hjust=1) +
    ggtitle("select annotated set") +
    guides(color="none") +
    xlim(-25, 250) +
    xlab("position on chromosome (mega base pairs)") +
    theme(axis.line.y=element_blank(),
      axis.text.y=element_blank(),
      axis.title.y=element_blank(),
      axis.title.x = element_text(size = 1, margin = 1),
      axis.ticks.y=element_blank()
      )+
    geom_point(
      aes((chromEnd+chromStart)/2/1e6, chr.int),
      data=chip.seq$set.info, 
      size=5,
      clickSelects="set.name"
    ),
  
  samples = ggplot() +
  theme_bw()+
    scale_x_continuous("false positive/negative rate (%)",
                      breaks=c(0, 50, 100), limits=c(-250, 100)) +
    theme_animint(width=350, height=300) +
    scale_fill_manual(values=fp.fn.colors) +
    ggtitle("select samples") +
    ylab("sample") +
    theme(axis.line.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank()),
  
  error = ggplot() +
    theme_animint(height=320) +
    theme(axis.line.y=element_blank(),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 10))+
    ggtitle("select model complexity and set") +
    xlab("model complexity -log(lambda)") +
    ylab("incorrect/total annotations (percent)") +
    geom_vline(
      aes(xintercept=complexity),
      data=chip.seq$models, 
      size=15, 
      alpha=1/2,
      clickSelects="complexity.i"
    ) +
    scale_color_manual(values=fp.fn.colors) +
    geom_line(
      aes(complexity, percent.error, linetype=set.name, group=set.name),
      data=chip.seq$error, 
      size=5, 
      alpha=3/4,
      clickSelects="set.name"
    ) +
    scale_linetype_manual(values=set.linetypes) +
    geom_line(
      aes(complexity, percent/2, group=error.type),
      data=fp.fn, 
      size=4.5, 
      color="grey",
      showSelected="set.name"
    ) +
    geom_line(
      aes(complexity, percent/2, color=error.type, group=error.type),
      data=fp.fn,
      size=1.5,
      showSelected="set.name"
    ) +
    scale_size_manual(values=c(false.negative=1.5, false.positive=3.5)) +
    geom_text(
      aes(6, 40, label=sprintf("%d/%d=%.1f%% false positive bases",
                               errors, bases, percent)),
      data=chip.seq$fp.set,
      showSelected=c("complexity.i", "set.name")
    ) +
    geom_text(
      aes(6, 33, label=sprintf("%d/%d=%.1f%% false negative bases",
                               errors, bases, percent)),
      data=chip.seq$fn.set,
      showSelected=c("complexity.i", "set.name")
    ) +
    geom_text(
      aes(6, 26, label=sprintf("%d/%d nonzero coefficients",
                               variables, total.vars)),
      data=chip.seq$nonzero,
      showSelected="complexity.i"
    ),
  
  roc = ggplot() +
    theme_animint(width=220, height=240) +
    scale_linetype_manual(values=set.linetypes) +
    guides(size="none", linetype="none") +
    geom_path(
      aes(FPR, TPR, group=set.name, linetype=set.name),
      data=chip.seq$roc.curves,
      size=3,
      showSelected="complexity.i",
      clickSelects="set.name"
    ) +
    geom_point(
      aes(FPR, TPR),
      data=chip.seq$roc.points,
      color="violet",
      size=5,
      showSelected="complexity.i",
      clickSelects="set.name"
    ) +
    scale_size_manual(values=c(train=3, validation=5)) +
    ggtitle("ROC curves") +
    xlab("False positive rate") +
    ylab("True positive rate"),
  
  probability = ggplot() +
    theme_animint(width=1300, height=300) +
    ggtitle("learned difference function") +
    theme(axis.line.x=element_blank(),
          axis.text.x=element_blank(),
          axis.title.x=element_blank(),
          axis.ticks.x=element_blank()) +
    geom_rect(
      aes(xmin=min.norm, xmax=max.norm, ymin=0, ymax=1, fill=annotation),
      data=chip.seq$regions,
      showSelected="set.name"
    ) +
    geom_text(
      aes(0.5, -0.08, label=paste0(width.bp, " bases on chr", chr)),size = 15,
      data=chip.seq$bases,
      showSelected="set.name"
    ) +
    scale_fill_manual(values=ann.colors) +
    geom_text(
      aes(mid.norm, 1.05, label=sprintf("%d/%d=%.1f%% false positive bases",
                                         errors, bases, percent)),size = 15,
      data=chip.seq$fp.pairs,
      showSelected=c("complexity.i", "sample1", "sample2", "set.name")
    ) +
    geom_text(
      aes(mid.norm, 1.05, label=sprintf("%d/%d=%.1f%% false negative bases",
                                         errors, bases, percent)),size = 15,
      data=chip.seq$fn.pairs,
      showSelected=c("complexity.i", "sample1", "sample2", "set.name")
    ) +
    geom_hline(yintercept=1/2, color="grey") +
    geom_ribbon(
      aes(mid.norm, ymin=min.prob, ymax=max.prob),
      data=chip.seq$probability,
      color="blue",
      showSelected=c("sample1", "sample2", "complexity.i", "set.name"),
      chunk_vars=c("sample1","sample2","complexity.i","set.name")
    ) +
    scale_y_continuous("probability of difference", breaks=c(0, 1/2, 1)),
  
  signal = ggplot() +
    scale_y_continuous("<---- one signal ----- another signal ->",
                       breaks=c(-1, 0, 1),
                       labels=c("max", "min", "max")) +
    theme_animint(width=1300, height=300) +
    ggtitle("ChIP-seq signal pair") +
    theme(axis.line.x=element_blank(),
          axis.text.x=element_blank(),
          axis.title.x=element_blank(),
          axis.ticks.x=element_blank()) +
    ylab("<---- one signal ----- another signal ->") +
    scale_fill_manual(values=ann.colors) +
    geom_rect(
      aes(xmin=min.norm, xmax=max.norm, ymin=-1, ymax=1, fill=annotation),
      data=chip.seq$regions,
      showSelected="set.name"
    ) +
    geom_text(
      aes(-0.07, 1.02, label=sprintf("%s %s max=%.1f", cell.type, sample1, max)),size = 15,
      data=chip.seq$signal.max$sample1,
      hjust=0,
      showSelected=c("set.name", "sample1")
    ) +
    geom_rect(
      aes(xmin=normStart, xmax=normEnd, ymin=0, ymax=signal.norm),
      data=chip.seq$signal.segments$sample1,
      size=0,
      showSelected=c("set.name", "sample1")
    ) +
    geom_text(
      aes(-0.07, -1, label=sprintf("%s %s max=%.1f", cell.type, sample2, max)),size = 15,
      data=chip.seq$signal.max$sample2,
      hjust=0,
      showSelected=c("set.name", "sample2")
    ) +
    geom_rect(
      aes(xmin=normStart, xmax=normEnd, ymin=-signal.norm, ymax=0),
      data=chip.seq$signal.segments$sample2,
      size=0,
      showSelected=c("set.name", "sample2")
    ) +
    geom_hline(yintercept=0, color="white"),
  
  selector.types = list(
    set.name = "single",
    sample1 = "single",
    sample2 = "single",
    complexity.i = "single"
  ),
  duration = list(complexity.i = 2000)
)

for(selector.name in names(chip.seq$samples)){
  sample.df <- chip.seq$samples[[selector.name]]
  sample.df$x <- -5
  y.fact <- if(selector.name=="sample1") -1 else 1
  other.name <- if(selector.name=="sample1") "sample2" else "sample1"
  sample.df$y <- sample.df$y * y.fact
  sample.id <- sample.df[[selector.name]]
  sample.df$label <- paste(sample.df$cell.type, sample.id)
  rownames(sample.df) <- sample.id
  nonzero.names <- as.character(nonzero.pairs[[selector.name]])
  nonzero.pairs$y <- sample.df[nonzero.names, "y"]
  sample.df$xmin <- 0
  sample.df$xmax <- 100
  sample.df$ymin <- sample.df$y-1/2
  sample.df$ymax <- sample.df$y+1/2
  text.df <- sample.df
  text.df$y <- text.df$y-1/2
  nonzero.pairs$key <- with(nonzero.pairs, {
    paste(sample1, sample2, error.type)
  })
  
  viz$samples <- viz$samples +
    geom_text(
      aes(x, y, label=label),
      data=text.df,
      hjust=1,
      clickSelects=selector.name,
      showSelected="set.name"
    ) +
    geom_rect(
      aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
      data=sample.df,
      alpha=1/2,
      showSelected="set.name",
      clickSelects=selector.name
    ) +
    geom_point(
      aes(percent, y, fill=error.type),
      data=nonzero.pairs,
      color="black",
      size=3,
      alpha=0.55,
      key="key",
      clickSelects=other.name,
      showSelected=c("set.name", "complexity.i")
    )
}

animint2dir(viz, "chip-seq-unaligned")