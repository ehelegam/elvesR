#' Analyze qPCR data in R
#'
#' It takes the your date and by computing accordingly to Pfaffl method it return the fold changes
#' @author Elves H Duarte
#' @param samples A column containing your samples name
#' @param genes A column with your reference and control genes
#' @param Cts Your Ct values
#' @param standev The maximum standard deviation among your technical replicates
#' @param sd.rm A logical parameter to decide whether the  outlayer sample will be removed
#' @param ref The reference samples
#' @param data the name of your data frame
#' @return Return a five-column data frame
#' @export

fold.calc=function(samples, genes, Cts, standev=standev, sd.rm=TRUE, ref, data)
{
  sd.values1=tapply(Cts[genes==paste(levels(genes)[1], sep="")], samples[genes==paste(levels(genes)[1], sep="")], list)
  sd.values2=tapply(Cts[genes==paste(levels(genes)[2], sep="")], samples[genes==paste(levels(genes)[2], sep="")], list)
  sd1=NULL
  sd2=NULL
  rm.sd1=NULL
  rm.sd2=NULL
  dt2.sd1=NULL
  dt2.sd2=NULL
  dtf=NULL
  for(i in 1:length(levels(samples)))
    {
    sd1[i]=sd(sd.values1[[i]])
    sd2[i]=sd(sd.values2[[i]])
  }
  {
  if(sd.rm==TRUE)
  {
    if(max(c(sd1, sd2), na.rm=T)>=standev)
    {
      dtf=NULL
      dt2=NULL
      if(max(sd1, na.rm=T)>=standev)
      {
        rm.sd1=which(sd1>=standev)
        for(i in rm.sd1)
        {
          mat1=combn(sd.values1[[i]], 2)
          val.sd1=which(apply(mat1, 2, sd)==min(apply(mat1, 2, sd)))
          val.rm=as.vector(mat1[, val.sd1])
          val.rm2=c(mat1)
          val.rm3=setdiff(val.rm2, val.rm)
          dt2.sd1[i]=match(val.rm3, data$Cq)
        }
        dt2.sd1=dt2.sd1[!is.na(dt2.sd1)]
      }
      if(max(sd2, na.rm=T)>=standev)
      {
        rm.sd2=which(sd2>=standev)
        for(i in rm.sd2)
        {
          mat2=combn(sd.values2[[i]], 2)
          val.sd2=which(apply(mat2, 2, sd)==min(apply(mat2, 2, sd)))
          val.rm=as.vector(mat2[, val.sd2])
          val.rm2=c(mat2)
          val.rm3=setdiff(val.rm2, val.rm)
          dt2.sd2[i]=match(val.rm3, data$Cq)
        }
        dt2.sd2=dt2.sd2[!is.na(dt2.sd2)]
      }
      dtf=c(dt2.sd1, dt2.sd2)
      dt2=data[c(-dtf),]
      write.table(dt2[c(dtf),], file="comments.txt", sep="\t", row.names=F, quote=F)
      write.table(dt2, file="data_corrected.txt", sep="\t", row.names=F, quote=F)
    }
    else
    {
      dt2=data
      write("No value was removed.", file="comments.txt")
    }
  }
    else
    {
  dt2=data
  write("You didn't test for variation among your samples.", file="comments.txt", sep="\t")
  }
  }
  result=matrix(1, ncol=5, nrow=length(levels(dt2$sample))+1)
  colnames(result)=c(paste(levels(dt2$gene)[2], c(".Mean", ".SD"), sep=""), paste(levels(dt2$gene)[1], c(".Mean", ".SD"), sep=""), "Fold change")
  rownames(result)=c(ref, levels(dt2$sample))
  ctr.gen1=tapply(dt2$Cq[dt2$gene==paste(levels(dt2$gene)[2], sep="")], dt2$sample[dt2$gene==paste(levels(dt2$gene)[2], sep="")], mean)
  ctr.gen2=tapply(dt2$Cq[dt2$gene==paste(levels(dt2$gene)[2], sep="")], dt2$sample[dt2$gene==paste(levels(dt2$gene)[2], sep="")], sd)
  tar.gen1=tapply(dt2$Cq[dt2$gene==paste(levels(dt2$gene)[1], sep="")], dt2$sample[dt2$gene==paste(levels(dt2$gene)[1], sep="")], mean)
  tar.gen2=tapply(dt2$Cq[dt2$gene==paste(levels(dt2$gene)[1], sep="")], dt2$sample[dt2$gene==paste(levels(dt2$gene)[1], sep="")], sd)
  dt3=subset(dt2, grepl(ref, sample))
  ref2=with(dt3, tapply(Cq, gene, mean))[[1]]
  ref1=with(dt3, tapply(Cq, gene, mean))[[2]]
  result[,1]=c(ref1, ctr.gen1)
  result[,2]=c(NA, ctr.gen2)
  result[,3]=c(ref2, tar.gen1)
  result[,4]=c(NA, tar.gen2)
  DCt.cal=result[-1,1]-result[1,1]
  DCt.test=result[-1,3]-result[1,3]
  DDCt=DCt.test-DCt.cal
  fold.dif=2^(-DDCt)
  result[,5]=c(1, fold.dif)
  write.table(round(result, 2), file="results.txt", se="\t", quote=F)
  return(round(result, 2))
}
