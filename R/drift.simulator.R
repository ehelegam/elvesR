#' Simulation drift with a Wolbachia population
#'
#' Drift simulation taking in account a population with to different phenotypes
#' @author Elves H Duarte
#' @param p Initial frequency of the population to be modeled
#' @param N Population size per generation
#' @param s Number of individuals to be sampled at each generation
#' @param Total generation to be modeled
#' @param fit Fitness of the population to be modeled. Zero (0) gives no fitness advantage.
#' @param r Number of replicates
#' @return A figure with different results.
#' @export
drift.sim=function(p, N, s, g, fit, r, get="fig"){
  drift=function(p, N, s, g, fit){
    prob=NULL
    if(p==0)
      stop("The initial frequency of wMelCS_b must be greather than 0.")
    prob[1]=p
    f=c(rep("wMel", p*N), rep("wMelCS_b", (1-p)*N*(1-fit)))
    #f=c(rep("wMel", p*N+N*fit), rep("wMelCS_b", (1-p)*N))
    f.sample=sample(f, s, replace=F)
    f.freq=table(f.sample)
    f.freq.fix=f.freq/sum(f.freq)
    f.percentage=ifelse(!("wMel" %in% names(f.freq)), 0, f.freq.fix[[1]])
    for(i in 2:(g-1))
    {
      prob[2]=f.percentage[[1]]
      f2=c(rep("wMel", prob[i]*N), rep("wMelCS_b", (1-prob[i])*N*(1-fit)))
      #f2=c(rep("wMel", prob[i]*N+N*fit), rep("wMelCS_b", (1-prob[i])*N))
      f.sample2=sample(f2, s, replace=F)
      f.freq2=table(f.sample2)
      f.percentage2=f.freq2/sum(f.freq2)
      prob[i+1]=ifelse(!("wMel" %in% names(f.freq2)), 0, f.percentage2[[1]])
    }
    prob
  }
  rep=replicate(r, drift(p, N, s, g, fit))
  {
  if(get=="fig")
  {
    tiff(paste("SIM", format(Sys.time(), "%b%d%Y-%H:%M:%S"), ".tiff", sep=""), width = 1000, height=600)
    def.par <- par(no.readonly = TRUE)
    par(mar=c(4, 5, 4, 2))
    layout(matrix(c(1, 1, 2, 2,  3 , 3, 3, 4), byrow=T, ncol=4), c(3,2), c(3,1), FALSE)
    plot(1:g, rep[,1], ylim=c(0, 1), type="l", col="black", xlim=c(1, g), main="Random drift model", xlab="Generations", ylab="Frequency", las=1, cex.axis=1.5, cex.lab=1.5, cex.main=1.5)
    points(1:g, rep[,1], pch=16, col="black", cex=0.5)
    grid(nx=NULL, ny=NA, lty=2, col="gray")
    for(i in 1:r-1)
    {
      colour="gray"
      lines(1:g, rep[,i+1], col=colour)
      points(1:g, rep[,i+1], col=colour, pch=16, cex=0.5)
    }
    abline(h=p, col="black", lty=2)
    freq=NULL
    prob.fix=NULL
    prob.lost=NULL
    for(i in 1:r)
    {
      prob.fix=c(which(rep[,i]==1)[1])
      prob.fix=as.numeric(prob.fix)
      prob.lost=c(which(rep[,i]==0)[1])
    }
    fix=sum(rep[g,]==1)
    lost=sum(rep[g,]==0)
    hist(rep[g,], col="gray", main="wMelCS_b distribution", ylab="n", xlab="% wMelCS_b", las=1, xlim=c(0,1), cex.axis=1.5, cex.lab=1.5, cex.main=1.5)
    grid(nx=NULL, ny=NULL, lty=2, col="gray")
    abline(v=p, col="red", cex=2)
    abline(v=mean(rep[g,]), lty=2, cex=2)
    #par(mar = c(3,3,2,1))
    plot(1:g, rowMeans(rep), ylim=c(0, 1), type="n", col="black", pch=16, xlim=c(1, g), main=paste("Average wMelCS_b=",max(round(rowMeans(rep), 2), sep="\t")), xlab="Generations", ylab="%", las=1, cex.axis=1.5, cex.lab=1.5, cex.main=1.5)
    stand=apply(rep, 1, sd)
    range1=apply(rep, 1, min)
    range2=apply(rep, 1, max)
    polygon(c(1:g, g:1), c(range1, rev(range2)), col=gray(0.8), border="gray")
    polygon(c(1:g, g:1), c(c(rowMeans(rep)+stand), rev(c(rowMeans(rep)-stand))), col=gray(0.4), border="gray")
    lines(rowMeans(rep), col="black", lty=3)
    points(rowMeans(rep), col="black", pch=16)
    barplot(c(lost/r*100, fix/r*100), names=c("0", 1), ylim=c(0,100), las=1, ylab="%", cex.axis=1.5, cex.names=1.5, cex.lab=1.5)
    par(def.par)
    dev.off()
    sink(paste("SIM", format(Sys.time(), "%b%d%Y-%H:%M:%S"), ".txt", sep=""))
    set.seed(12345)
    cat("=============================\n")
    cat("Parameters of the model\n")
    cat("=============================\n")
    cat("Initial frequency: ", p, "\n")
    cat("Population size: ", N, "\n")
    cat("Population effective size: ", s, "\n")
    cat("Generations: ", g, "\n")
    cat("Relative fitness: ", fit, "\n")
    cat("Independent replicates: ", r, "\n")
    cat("\n")
    cat("=============================\n")
    cat("Results\n")
    cat("=============================\n")
    cat(sprintf("Average frequency of wMelCS_b after %d generations: ", g), mean(rep[g,]), "\n")
    cat(sprintf("Proportion of times wMelCS_b were lost: %d%s\n",lost/r*100, "%"))
    cat(sprintf("Proportion of times wMelCS_b were fixed: %d%s\n",fix/r*100, "%"))
    #cat("Time to be lost: ", match(0, t(rep)), "\n")
    sink()
    file.show(paste("SIM", format(Sys.time(), "%b%d%Y-%H:%M:%S"), ".txt", sep=""))
  }
  else
  {
    return(rep)
  }
  }
}
