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
  {
    if(length(x)==1)
    {
      for(i in 1:length(genes$V4))
      {
        rg=genes[i,2]:genes[i,3]
        {
          if(length(x)==1)
          {
            if(x %in% rg)
            {
              gene_list=as.character(genes[i,4])
              break
            }
            else
            {
              gene_list="non_coding"
            }
          }
        }
      }
    }
    else
    {
      gene_list=NULL
      for(n in 1:length(x))
      {
        for(i in 1:length(genes$V4))
        {
          rg=genes[i,2]:genes[i,3]
          {
            if(x[n] %in% rg)
            {
              gene_list[n]=as.character(genes[i,4])
              break
            }
            else
            {
              gene_list[n]="non_coding"
            }
          }
        }
      }
    }
  }
  return(gene_list)
}
