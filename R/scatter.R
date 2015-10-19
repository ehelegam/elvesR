#' Scatter plots and median or mean bars
#'
#' This function creats a scatter plot and add groups bar representing median or mean values
#' @author Elves H Duarte
#' @param x The values to be plotted
#' @param y Categorical variable containing the different groups
#' @param fun Choose which statistic (medium or mean) will be added
#' @return A plot based in 1-D scatter plot of the base package
#' @export


scatter<- function(x, data, fun="median",
                   main="", sub="",
                   xlim=NULL, ylim=NULL,
                   axes=TRUE, ann=par("ann"),
                   col=par("col"),
                   ...) 
{
  if(class(x)!="formula")
  {
    if(fun!="median" & fun!="mean")
      stop("Only 'median' or 'mean' are currently available.")
    if(fun=="median")
      fun=median(x, na.rm=T)
    if(fun=="mean")
      fun=mean(x, na.rm=T)
    if (is.null(xlim))
      xlim <- range(0.5:1.5)
    if (is.null(ylim))
      ylim <- range(x, na.rm=T)
    par(mar=c(5,4.5,4,2)+0.1)
    #plot.new()
    plot.window(xlim, ylim, ...)
    stripchart(x, vertical=T, method="jitter", jitter=0.075, col=col, xlim=xlim, ylim=ylim, ...)
    segments(1:length(x)-0.15, fun, 1:length(x)+0.15, fun, col="black", lwd=4, xlim=NULL)
    if (ann)
      title(main=main, sub=sub, ...)
  }
  if(class(x)=="formula")
  {
    new.data = model.frame(formula=x, data=data)
    #return(new.data)
    if(fun!="median" & fun!="mean")
      stop("Only 'median' or 'mean' are currently available.")
    else
      fun=tapply(new.data[,1], list(new.data[,2]), fun, na.rm=T)
    if (is.null(xlim))
      xlim <- range(c(0:length(levels(new.data[,2])))+0.5)
    if (is.null(ylim))
      ylim <- range(new.data[,1], na.rm=T)
    par(mar=c(5,4.5,4,2)+0.1)
    #plot.new()
    plot.window(xlim, ylim, ...)
    if(length(new.data)>2)
    {
      stripchart(new.data[,1]~new.data[,2], vertical=T, method="jitter", jitter=0.075, col=col, xlim=xlim, ylim=ylim, subset=new.data[,3]==levels(new.data[,3][1]),...)
      segments(1:length(new.data[,2])-0.15, fun, 1:length(new.data[,2])+0.15, fun, col="black", lwd=4, xlim=NULL)
      if(ann)
        title(main=main, sub=sub, ...) 
    }
    else
    {
      stripchart(new.data[,1]~new.data[,2], vertical=T, method="jitter", jitter=0.075, col=col, xlim=xlim, ylim=ylim, ...)
      segments(1:length(new.data[,2])-0.15, fun, 1:length(new.data[,2])+0.15, fun, col="black", lwd=4, xlim=NULL)
      if(ann)
        title(main=main, sub=sub, ...) 
    }
  }
}
