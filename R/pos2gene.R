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
      gen=subset(.genes, .genes$V5>=x)[1,]
      exon=gen$V4:gen$V5
      if(x %in% exon)
        gene_list=as.character(gen[1,6])
      else
        gene_list="non-coding"
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
        gen=subset(.genes, .genes$V5>=x[n])[1,]
        exon=gen$V4:gen$V5
        if(x[n] %in% exon)
          gene_list[n]=as.character(gen[1,6])
        else
          gene_list[n]="non-coding"
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
