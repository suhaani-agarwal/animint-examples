works_with_R(
  "3.2.2",
  ggplot2="2.1.0",
  PeakSegJoint="2016.3.2",
  "tdhock/animint@03735869af84629d269556442345b2ea506ab42a")

load("../data/PSJ.RData")

res.error <- PSJ$error.total.chunk

ann.colors <-
  c(noPeaks="#f6f4bf",
    peakStart="#ffafaf",
    peakEnd="#ff4c4c",
    peaks="#a445ee")

## prob.regions are the black segments that show which regions are
## mapped to which segmentation problems.
all.regions <- do.call(rbind, PSJ$regions.by.problem)
prob.regions.names <-
  c("bases.per.problem", "problem.i", "problem.name",
    "chromStart", "chromEnd")
prob.regions <- unique(data.frame(all.regions)[, prob.regions.names])
prob.regions$sample.id <- "problems"

all.modelSelection <- do.call(rbind, PSJ$modelSelection.by.problem)
modelSelection.errors <-
  all.modelSelection[!is.na(all.modelSelection$errors), ]
penalty.range <-
  with(all.modelSelection, c(min(max.log.lambda), max(min.log.lambda)))
penalty.mid <- mean(penalty.range)

coverage.counts <- table(PSJ$coverage$sample.id)
facet.rows <- length(coverage.counts)+1
dvec <- diff(log(res.error$bases.per.problem))
dval <- exp(mean(dvec))
dval2 <- (dval-1)/2 + 1
res.error$min.bases.per.problem <- res.error$bases.per.problem/dval2
res.error$max.bases.per.problem <- res.error$bases.per.problem*dval2

modelSelection.labels <- unique(with(all.modelSelection, {
  data.frame(problem.name=problem.name,
             bases.per.problem=bases.per.problem,
             problemStart=problemStart,
             problemEnd=problemEnd,
             min.log.lambda=penalty.mid,
             peaks=max(peaks)+0.5)
}))

cat("constructing data viz with for loops\n")
print(system.time({
  viz.for <-
    list(coverage=ggplot()+
           geom_segment(aes(chromStart/1e3, problem.i,
                            xend=chromEnd/1e3, yend=problem.i,
                            showSelected=bases.per.problem,
                            clickSelects=problem.name),
                        data=prob.regions)+
           ggtitle("select problem")+
           geom_text(aes(chromStart/1e3, problem.i,
                         showSelected=bases.per.problem,
                         label=sprintf("%d problems mean size %.1f kb",
                           problems, mean.bases/1e3)),
                     data=PSJ$problem.labels,
                     hjust=0)+
           geom_segment(aes(problemStart/1e3, problem.i,
                            showSelected=bases.per.problem,
                            clickSelects=problem.name,
                            xend=problemEnd/1e3, yend=problem.i),
                        size=5,
                        data=PSJ$problems)+
           scale_y_continuous("aligned read coverage",
                              breaks=function(limits){
                                floor(limits[2])
                              })+
           scale_linetype_manual("error type",
                                 limits=c("correct", 
                                   "false negative",
                                   "false positive"
                                          ),
                                 values=c(correct=0,
                                   "false negative"=3,
                                   "false positive"=1))+
           scale_x_continuous(paste("position on chr11",
                                    "(kilo bases = kb)"))+
           coord_cartesian(xlim=c(118167.406, 118238.833))+
           geom_tallrect(aes(xmin=chromStart/1e3, xmax=chromEnd/1e3,
                             fill=annotation),
                         alpha=0.5,
                         color="grey",
                         data=PSJ$filled.regions)+
           scale_fill_manual(values=ann.colors)+
           theme_bw()+
           theme_animint(width=1500, height=facet.rows*100)+
           theme(panel.margin=grid::unit(0, "cm"))+
           facet_grid(sample.id ~ ., labeller=function(df){
             df$sample.id <- sub("McGill0", "", sub(" ", "\n", df$sample.id))
             df
           }, scales="free")+
           geom_line(aes(base/1e3, count),
                     data=PSJ$coverage,
                     color="grey50"),

         resError=ggplot()+
           ggtitle("select problem size")+
           ylab("minimum percent incorrect regions")+
           geom_tallrect(aes(xmin=min.bases.per.problem,
                             xmax=max.bases.per.problem,
                             clickSelects=bases.per.problem),
                         alpha=0.5,
                         data=res.error)+
           scale_x_log10()+
           geom_line(aes(bases.per.problem, errors/regions*100,
                         color=chunks, size=chunks),
                     data=data.frame(res.error, chunks="this"))+
           geom_line(aes(bases.per.problem, errors/regions*100,
                         color=chunks, size=chunks),
                     data=data.frame(PSJ$error.total.all, chunks="all")),

         modelSelection=ggplot()+
           geom_segment(aes(min.log.lambda, peaks,
                            xend=max.log.lambda, yend=peaks,
                            showSelected=problem.name,
                            showSelected2=bases.per.problem),
                        data=data.frame(all.modelSelection, what="peaks"),
                        size=5)+
           geom_text(aes(min.log.lambda, peaks,
                         showSelected=problem.name,
                         showSelected2=bases.per.problem,
                         label=sprintf("%.1f kb in problem %s",
                           (problemEnd-problemStart)/1e3, problem.name)),
                     data=data.frame(modelSelection.labels, what="peaks"))+
           geom_segment(aes(min.log.lambda, as.integer(errors),
                            xend=max.log.lambda, yend=as.integer(errors),
                            showSelected=problem.name,
                            showSelected2=bases.per.problem),
                        data=data.frame(modelSelection.errors, what="errors"),
                        size=5)+
           ggtitle("select number of samples with 1 peak")+
           ylab("")+
           facet_grid(what ~ ., scales="free"),
         
         title="Animint compiler with for loops",

         first=PSJ$first)

  ## For every problem there is a selector (called problem.dot) for the
  ## number of peaks in that problem. So in this for loop we add a few
  ## layers with aes_string(clickSelects=problem.dot) or
  ## aes_string(showSelected=problem.dot) to the coverage and
  ## modelSelection plots.
  for(problem.dot in names(PSJ$modelSelection.by.problem)){
    regions.dt <- PSJ$regions.by.problem[[problem.dot]]
    regions.dt[[problem.dot]] <- regions.dt$peaks
    if(!is.null(regions.dt)){
      viz.for$coverage <- viz.for$coverage+
        geom_tallrect(aes_string(xmin="chromStart/1e3",
                                 xmax="chromEnd/1e3",
                                 linetype="status",
                                 showSelected=problem.dot,
                                 showSelected2="bases.per.problem"),
                      data=data.frame(regions.dt),
                      fill=NA,
                      color="black")
    }
    if(problem.dot %in% names(PSJ$peaks.by.problem)){
      peaks <- PSJ$peaks.by.problem[[problem.dot]]
      peaks[[problem.dot]] <- peaks$peaks
      prob.peaks.names <-
        c("bases.per.problem", "problem.i", "problem.name",
          "chromStart", "chromEnd", problem.dot)
      prob.peaks <- unique(data.frame(peaks)[, prob.peaks.names])
      prob.peaks$sample.id <- "problems"
      viz.for$coverage <- viz.for$coverage +
        geom_segment(aes_string("chromStart/1e3", "0",
                                xend="chromEnd/1e3", yend="0",
                                clickSelects="problem.name",
                                showSelected=problem.dot,
                                showSelected2="bases.per.problem"),
                     data=peaks, size=7, color="deepskyblue")+
        geom_segment(aes_string("chromStart/1e3", "problem.i",
                                xend="chromEnd/1e3", yend="problem.i",
                                clickSelects="problem.name",
                                showSelected=problem.dot,
                                showSelected2="bases.per.problem"),
                     data=prob.peaks, size=7, color="deepskyblue")
    }
    modelSelection.dt <- PSJ$modelSelection.by.problem[[problem.dot]]
    modelSelection.dt[[problem.dot]] <- modelSelection.dt$peaks
    viz.for$modelSelection <- viz.for$modelSelection+
      geom_tallrect(aes_string(xmin="min.log.lambda", 
                               xmax="max.log.lambda", 
                               clickSelects=problem.dot,
                               showSelected="problem.name",
                               showSelected2="bases.per.problem"),
                    data=modelSelection.dt, alpha=0.5)
  }
}))

cat("compiling data viz\n")
print(system.time({
  animint2dir(viz.for, out.dir="PSJ-for-new")
}))

sample.peaks <- do.call(rbind, PSJ$peaks.by.problem)
prob.peaks.names <-
  c("bases.per.problem", "problem.i", "problem.name", "peaks",
    "chromStart", "chromEnd")
problem.peaks <- unique(data.frame(sample.peaks)[, prob.peaks.names])
problem.peaks$sample.id <- "problems"

peakvar <- function(position){
  paste0(gsub("[-:]", ".", position), "peaks")
}

cat("constructing data viz with .variable .value\n")
print(system.time({
  viz <-
    list(coverage=ggplot()+
           geom_segment(aes(chromStart/1e3, problem.i,
                            xend=chromEnd/1e3, yend=problem.i,
                            showSelected=bases.per.problem,
                            clickSelects=problem.name),
                        data=prob.regions)+
           ggtitle("select problem")+
           geom_text(aes(chromStart/1e3, problem.i,
                         showSelected=bases.per.problem,
                         label=sprintf("%d problems mean size %.1f kb",
                           problems, mean.bases/1e3)),
                     data=PSJ$problem.labels,
                     hjust=0)+
           geom_segment(aes(problemStart/1e3, problem.i,
                            showSelected=bases.per.problem,
                            clickSelects=problem.name,
                            xend=problemEnd/1e3, yend=problem.i),
                        size=5,
                        data=PSJ$problems)+
           scale_y_continuous("aligned read coverage",
                              breaks=function(limits){
                                floor(limits[2])
                              })+
           scale_linetype_manual("error type",
                                 limits=c("correct", 
                                   "false negative",
                                   "false positive"
                                          ),
                                 values=c(correct=0,
                                   "false negative"=3,
                                   "false positive"=1))+
           scale_x_continuous(paste("position on chr11",
                                    "(kilo bases = kb)"))+
           coord_cartesian(xlim=c(118167.406, 118238.833))+
           geom_tallrect(aes(xmin=chromStart/1e3, xmax=chromEnd/1e3,
                             fill=annotation),
                         alpha=0.5,
                         color="grey",
                         data=PSJ$filled.regions)+
           scale_fill_manual(values=ann.colors)+
           theme_bw()+
           theme_animint(width=1500, height=facet.rows*100)+
           theme(panel.margin=grid::unit(0, "cm"))+
           facet_grid(sample.id ~ ., labeller=function(df){
             df$sample.id <- sub("McGill0", "", sub(" ", "\n", df$sample.id))
             df
           }, scales="free")+
           geom_line(aes(base/1e3, count),
                     data=PSJ$coverage,
                     color="grey50")+
           geom_tallrect(aes(xmin=chromStart/1e3,
                             xmax=chromEnd/1e3,
                             linetype=status,
                             showSelected.value=peaks,
                             showSelected.variable=peakvar(problem.name),
                             showSelected2=bases.per.problem),
                         data=all.regions,
                         fill=NA,
                         color="black")+
           geom_segment(aes(chromStart/1e3, 0,
                            xend=chromEnd/1e3, yend=0,
                            clickSelects=problem.name,
                            showSelected.variable=peakvar(problem.name),
                            showSelected.value=peaks,
                            showSelected2=bases.per.problem),
                        data=sample.peaks, size=7, color="deepskyblue")+
           geom_segment(aes(chromStart/1e3, problem.i,
                            xend=chromEnd/1e3, yend=problem.i,
                            clickSelects=problem.name,
                            showSelected.variable=peakvar(problem.name),
                            showSelected.value=peaks,
                            showSelected2=bases.per.problem),
                        data=problem.peaks, size=7, color="deepskyblue"),

         resError=ggplot()+
           ggtitle("select problem size")+
           ylab("minimum percent incorrect regions")+
           geom_tallrect(aes(xmin=min.bases.per.problem,
                             xmax=max.bases.per.problem,
                             clickSelects=bases.per.problem),
                         alpha=0.5,
                         data=res.error)+
           scale_x_log10()+
           geom_line(aes(bases.per.problem, errors/regions*100,
                         color=chunks, size=chunks),
                     data=data.frame(res.error, chunks="this"))+
           geom_line(aes(bases.per.problem, errors/regions*100,
                         color=chunks, size=chunks),
                     data=data.frame(PSJ$error.total.all, chunks="all")),

         modelSelection=ggplot()+
           geom_segment(aes(min.log.lambda, peaks,
                            xend=max.log.lambda, yend=peaks,
                            showSelected=problem.name,
                            showSelected2=bases.per.problem),
                        data=data.frame(all.modelSelection, what="peaks"),
                        size=5)+
           geom_text(aes(min.log.lambda, peaks,
                         showSelected=problem.name,
                         showSelected2=bases.per.problem,
                         label=sprintf("%.1f kb in problem %s",
                           (problemEnd-problemStart)/1e3, problem.name)),
                     data=data.frame(modelSelection.labels, what="peaks"))+
           geom_segment(aes(min.log.lambda, as.integer(errors),
                            xend=max.log.lambda, yend=as.integer(errors),
                            showSelected=problem.name,
                            showSelected2=bases.per.problem),
                        data=data.frame(modelSelection.errors, what="errors"),
                        size=5)+
           ggtitle("select number of samples with 1 peak")+
           ylab("")+
           geom_tallrect(aes(xmin=min.log.lambda, 
                             xmax=max.log.lambda, 
                             clickSelects.variable=
                               peakvar(problem.name),
                             clickSelects.value=peaks,
                             showSelected=problem.name,
                             showSelected2=bases.per.problem),
                         data=all.modelSelection, alpha=0.5)+
           facet_grid(what ~ ., scales="free"),
         
         title="Animint compiler with .variable .value aesthetics",

         first=PSJ$first)
### For every problem there is a selector (called problem.name) for
### the number of peaks in that problem. The animint2dir compiler
### creates a selection variable for every unique value of
### clickSelects.variable and showSelected.variable (and it uses
### clickSelects.value and showSelected.value to set/update the
### selected value/geoms).
}))

## Timings show that plot construction is must faster when not using
## a for loop.

## constructing data viz with for loops
##    user  system elapsed 
##  37.340   0.252  37.873

## constructing data viz with .variable .value
##    user  system elapsed 
##   0.104   0.000   0.105

## Constructing the plot using the for loop above makes 1180 geoms
## which are chunked into many small tsv files, which results in a lot
## of wasted disk space:

## ~/R/animint-examples/examples $ du -hs PSJ-for/
## 8.3M	PSJ-for/
## ~/R/animint-examples/examples $ ls PSJ-for/|wc -l
## 1767
## ~/R/animint-examples/examples $ ls PSJ-for/geom*.tsv|wc -l
## 1760

## ~/R/animint-examples/examples $ du --block-size=1 -cs PSJ-for/*.tsv|tail
## 4096	PSJ-for/geom98_segment_coverage_chunk1.tsv
## 4096	PSJ-for/geom99_segment_coverage_chunk1.tsv
## 4096	PSJ-for/geom99_segment_coverage_chunk2.tsv
## 4096	PSJ-for/geom99_segment_coverage_chunk3.tsv
## 4096	PSJ-for/geom99_segment_coverage_chunk4.tsv
## 4096	PSJ-for/geom9_segment_coverage_chunk1.tsv
## 4096	PSJ-for/geom9_segment_coverage_chunk2.tsv
## 4096	PSJ-for/geom9_segment_coverage_chunk3.tsv
## 4096	PSJ-for/geom9_segment_coverage_chunk4.tsv
## 7471104	total

## ~/R/animint-examples/examples $ du --block-size=1 -cs --apparent-size PSJ-for/*.tsv|tail
## 480	PSJ-for/geom98_segment_coverage_chunk1.tsv
## 154	PSJ-for/geom99_segment_coverage_chunk1.tsv
## 261	PSJ-for/geom99_segment_coverage_chunk2.tsv
## 368	PSJ-for/geom99_segment_coverage_chunk3.tsv
## 475	PSJ-for/geom99_segment_coverage_chunk4.tsv
## 153	PSJ-for/geom9_segment_coverage_chunk1.tsv
## 259	PSJ-for/geom9_segment_coverage_chunk2.tsv
## 365	PSJ-for/geom9_segment_coverage_chunk3.tsv
## 471	PSJ-for/geom9_segment_coverage_chunk4.tsv
## 961044	total

## On my filesystem each tsv file takes 4096 bytes even though its
## apparent size is less. The compiler should not create files less
## than 4KB.

## After gzipping these files, most stay at 4096 bytes, but the bigger
## ones compress:

## thocking@silene:~/R/animint-examples/examples/PSJ-for(master)$ du --block-size=1 -c *.tsv|grep -v 4096
## 8192	geom3_segment_coverage_chunk1.tsv
## 20480	geom465_text_modelSelection_chunk1.tsv
## 245760	geom5_line_coverage_chunk1.tsv
## 7471104	total

## thocking@silene:~/R/animint-examples/examples/PSJ-for(master)$ gzip *.tsv
## thocking@silene:~/R/animint-examples/examples/PSJ-for(master)$ du --block-size=1 -c geom3_segment_coverage_chunk1.tsv.gz geom465_text_modelSelection_chunk1.tsv.gz geom5_line_coverage_chunk1.tsv.gz 
## 4096	geom3_segment_coverage_chunk1.tsv.gz
## 4096	geom465_text_modelSelection_chunk1.tsv.gz
## 61440	geom5_line_coverage_chunk1.tsv.gz
## 69632	total

## thocking@silene:~/R/animint-examples/examples/PSJ-for(master)$ du --block-size=1 -c *.gz|tail
## 4096	geom98_segment_coverage_chunk1.tsv.gz
## 4096	geom99_segment_coverage_chunk1.tsv.gz
## 4096	geom99_segment_coverage_chunk2.tsv.gz
## 4096	geom99_segment_coverage_chunk3.tsv.gz
## 4096	geom99_segment_coverage_chunk4.tsv.gz
## 4096	geom9_segment_coverage_chunk1.tsv.gz
## 4096	geom9_segment_coverage_chunk2.tsv.gz
## 4096	geom9_segment_coverage_chunk3.tsv.gz
## 4096	geom9_segment_coverage_chunk4.tsv.gz
## 7266304	total

## In this case the savings is only 200KB, but there are likely to be
## more savings if we switch to making a smaller number of larger
## files.

## However, for serving these files we need to keep them as
## uncompressed tsv. For file transfer purposes the web server (for
## example gist.github.com) will automatically compress the file and
## serve it with headers
## Content-Encoding: gzip
## Content-Type: text/tab-separated-values

cat("compiling data viz\n")
print(system.time({
  animint2dir(viz, out.dir="PSJ")
}))

not.problem.name <- names(PSJ$first)[names(PSJ$first) != "problem.name"]

problems.by.res <- split(PSJ$problems, PSJ$problems$bases.per.problem)

## get_ functions are passed all of the current values of the selector
## variables, and must return a data.frame to display.
L <- lapply(PSJ$first, paste) # for example, the first selection
get_tallrect <- function(...){
  L <- list(...)
  res.problems <- problems.by.res[[L$bases.per.problem]]
  problem.name.vec <- paste(res.problems$problem.name)
  peaks.by.problem <- list()
  for(problem.name in problem.name.vec){
    selector.name <- peakvar(problem.name)
    prob.peaks <- PSJ$peaks.by.problem[[selector.name]]
    peaks.by.problem[[problem.name]] <-
      subset(prob.peaks, peaks == L[[selector.name]])
  }
  show.peaks <- do.call(rbind, peaks.by.problem)
  error.regions <- PeakErrorSamples(show.peaks, PSJ$filled.regions)
  error.regions
}

## For example, the second problem in the selected resolution.
problem.name <- "chr11:118174946-118177139"
bases.per.problem <- "6516"
get_segment <- function(problem.name, bases.per.problem, ...){
  L <- list(...)
  res.problems <- problems.by.res[[bases.per.problem]]
  is.selected <- res.problems$problem.name == problem.name
  other.problems <- res.problems[!is.selected, ]
  peaks.by.problem <- list()
  for(other.name in other.problems$problem.name){
    selector.name <- peakvar(other.name)
    prob.peaks <- PSJ$peaks.by.problem[[selector.name]]
    peaks.by.problem[[other.name]] <-
      subset(prob.peaks, peaks == L[[selector.name]])
  }
  other.peaks <- do.call(rbind, peaks.by.problem)
  
  selector.name <- peakvar(problem.name)
  prob.peaks <- PSJ$peaks.by.problem[[selector.name]]
  modelSelection <- PSJ$modelSelection.by.problem[[selector.name]]
  modelSelection$errors <- NA
  for(model.i in 1:nrow(modelSelection)){
    model.row <- modelSelection[model.i, ]
    model.peaks <- subset(prob.peaks, peaks == model.row$peaks)
    all.peaks <- rbind(other.peaks, model.peaks)
    ## TODO: this is the same computation as in get_tallrect, so why
    ## don't we just have it used the cached value from get_tallrect?
    error.regions <- PeakErrorSamples(all.peaks, PSJ$filled.regions)
    modelSelection$errors[model.i] <- with(error.regions, sum(fp+fn))
  }
  modelSelection
}

## Cache functions are passed all the current values of the selector
## variables, and return a character vector of selector variable names
## whose values will be used for saving the computed data.
cache_tallrect <- function(bases.per.problem){
  res.problems <- problems.by.res[[bases.per.problem]]
  peakvar(res.problems$problem.name)
}
cache_segment <- function(bases.per.problem, ...){
  res.problems <- problems.by.res[[bases.per.problem]]
  is.selected <- res.problems$problem.name == problem.name
  other.problems <- res.problems[!is.selected, ]
  peakvar(other.problems$problem.name)
}

cat("constructing data viz with lookup function\n")
print(system.time({
  viz <-
    list(coverage=ggplot()+
           geom_segment(aes(chromStart/1e3, problem.i,
                            xend=chromEnd/1e3, yend=problem.i,
                            showSelected=bases.per.problem,
                            clickSelects=problem.name),
                        data=prob.regions)+
           ggtitle("select problem")+
           geom_text(aes(chromStart/1e3, problem.i,
                         showSelected=bases.per.problem,
                         label=sprintf("%d problems mean size %.1f kb",
                           problems, mean.bases/1e3)),
                     data=PSJ$problem.labels,
                     hjust=0)+
           geom_segment(aes(problemStart/1e3, problem.i,
                            showSelected=bases.per.problem,
                            clickSelects=problem.name,
                            xend=problemEnd/1e3, yend=problem.i),
                        size=5,
                        data=PSJ$problems)+
           scale_y_continuous("aligned read coverage",
                              breaks=function(limits){
                                floor(limits[2])
                              })+
           scale_linetype_manual("error type",
                                 limits=c("correct", 
                                   "false negative",
                                   "false positive"
                                          ),
                                 values=c(correct=0,
                                   "false negative"=3,
                                   "false positive"=1))+
           scale_x_continuous(paste("position on chr11",
                                    "(kilo bases = kb)"))+
           coord_cartesian(xlim=c(118167.406, 118238.833))+
           geom_tallrect(aes(xmin=chromStart/1e3, xmax=chromEnd/1e3,
                             fill=annotation),
                         alpha=0.5,
                         color="grey",
                         data=PSJ$filled.regions)+
           scale_fill_manual(values=ann.colors)+
           theme_bw()+
           theme_animint(width=1500, height=facet.rows*100)+
           theme(panel.margin=grid::unit(0, "cm"))+
           facet_grid(sample.id ~ ., labeller=function(df){
             df$sample.id <- sub("McGill0", "", sub(" ", "\n", df$sample.id))
             df
           }, scales="free")+
           geom_line(aes(base/1e3, count),
                     data=PSJ$coverage,
                     color="grey50")+
           geom_tallrect(aes(xmin=chromStart/1e3,
                             xmax=chromEnd/1e3,
                             linetype=status,
                             showSelected.value=peaks,
                             showSelected.variable=peakvar(problem.name),
                             showSelected2=bases.per.problem),
                         updateWhenChanged=not.problem.name,
                         cache=cache_tallrect,
                         data=data.frame(get_tallrect),
                         fill=NA,
                         color="black")+
           geom_segment(aes(chromStart/1e3, 0,
                            xend=chromEnd/1e3, yend=0,
                            clickSelects=problem.name,
                            showSelected.variable=peakvar(problem.name),
                            showSelected.value=peaks,
                            showSelected2=bases.per.problem),
                        data=sample.peaks, size=7, color="deepskyblue")+
           geom_segment(aes(chromStart/1e3, problem.i,
                            xend=chromEnd/1e3, yend=problem.i,
                            clickSelects=problem.name,
                            showSelected.variable=peakvar(problem.name),
                            showSelected.value=peaks,
                            showSelected2=bases.per.problem),
                        data=problem.peaks, size=7, color="deepskyblue"),

         resError=ggplot()+
           ggtitle("select problem size")+
           ylab("minimum percent incorrect regions")+
           geom_tallrect(aes(xmin=min.bases.per.problem,
                             xmax=max.bases.per.problem,
                             clickSelects=bases.per.problem),
                         alpha=0.5,
                         data=res.error)+
           scale_x_log10()+
           geom_line(aes(bases.per.problem, errors/regions*100,
                         color=chunks, size=chunks),
                     data=data.frame(res.error, chunks="this"))+
           geom_line(aes(bases.per.problem, errors/regions*100,
                         color=chunks, size=chunks),
                     data=data.frame(PSJ$error.total.all, chunks="all")),

         modelSelection=ggplot()+
           geom_segment(aes(min.log.lambda, peaks,
                            xend=max.log.lambda, yend=peaks,
                            showSelected=problem.name,
                            showSelected2=bases.per.problem),
                        data=data.frame(all.modelSelection, what="peaks"),
                        size=5)+
           geom_text(aes(min.log.lambda, peaks,
                         showSelected=problem.name,
                         showSelected2=bases.per.problem,
                         label=sprintf("%.1f kb in problem %s",
                           (problemEnd-problemStart)/1e3, problem.name)),
                     data=data.frame(modelSelection.labels, what="peaks"))+
           geom_segment(aes(min.log.lambda, as.integer(errors),
                            xend=max.log.lambda, yend=as.integer(errors),
                            showSelected=problem.name,
                            showSelected2=bases.per.problem),
                        data=data.frame(get_segment),
                        cache=cache_segment,
                        updateWhenChanged=
                          c("problem.name", "bases.per.problem"),
                        size=5)+
           ggtitle("select number of samples with 1 peak")+
           ylab("")+
           geom_tallrect(aes(xmin=min.log.lambda, 
                             xmax=max.log.lambda, 
                             clickSelects.variable=
                               peakvar(problem.name),
                             clickSelects.value=peaks,
                             showSelected=problem.name,
                             showSelected2=bases.per.problem),
                         data=all.modelSelection, alpha=0.5)+
           facet_grid(what ~ ., scales="free"),
         
         title="Animint compiler with .variable .value aesthetics",

         first=PSJ$first)
}))

cat("compiling data viz\n")
print(system.time({
  animint2dir(viz, out.dir="PSJ")
}))

