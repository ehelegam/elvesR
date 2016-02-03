#' Annotate Wolbachia positions
#'
#' By using a bed file each position of Wolbahcia is annotated
#' @author Elves H Duarte
#' @param x A list of positions to be annotated.
#' @return Retur a list of genes (or not) of Wolbachia locus.
#' @export

pos2gene=function(x)
{
  data("wolb_genes")
  .chr=1:1268079
  {
    if(length(x)==1)
    {
      ##
      if(x %in% .chr)
      {
        for(i in 1:length(.genes$V3))
        {
          rg=.genes[i,4]:.genes[i,5]
          if(length(x)==1)
          {
            if(x %in% rg)
            {
              gene_list=as.character(.genes[i,6])
              break
            }
            else
            {
              gene_list="non_coding"
            }
      }
      }
      }
      else
      {
        gene_list=NA
      }
    }
    else
    {
      gene_list=NULL
      for(n in 1:length(x))
      {
        if(x[n] %in% .chr)
        {
          for(i in 1:length(.genes$V3))
          {
            rg=.genes[i,4]:.genes[i,5]
            if(x[n] %in% rg)
            {
              gene_list[n]=as.character(.genes[i,6])
              break
            }
            else
            {
              gene_list[n]="non_coding"
            }
          }
        }
        else
        {
          gene_list[n]=NA
        }
      }
    }
  }
  return(gene_list)
}
