standardize_celltypes <- function(celltype_vector,
                                  allow_subtypes=F){
  if(allow_subtypes){
    Cell_group_dict <- c("astrocytes"="astrocytes",
                         "microglia"="microglia",
                         "oligo"="oligo",
                         "neurons"="neurons",
                         "brain"="brain",
                         ExcitatoryNeurons="neurons (+)",
                         InhibitoryNeurons="neurons (-)",
                         NigralNeurons="neurons (nigral)",
                         Microglia="microglia",
                         Oligodendrocytes="oligo",
                         Astrocytes="astrocytes",
                         OPCs="OPCs")
  }else {
    Cell_group_dict <- c("astrocytes"="astrocytes",
                         "microglia"="microglia",
                         "oligo"="oligo",
                         "OPCs"="oligo",
                         "neurons"="neurons",
                         "neurons (+)"="neurons",
                         "neurons (-)"="neurons",
                         "neurons (nigral)"="neurons",
                         "brain"="brain",
                         ExcitatoryNeurons="neurons",
                         InhibitoryNeurons="neurons",
                         NigralNeurons="neurons",
                         Microglia="microglia",
                         Oligodendrocytes="oligo",
                         Astrocytes="astrocytes")
  }
  celltypes <- Cell_group_dict[celltype_vector]
  return(celltypes)
}

