dep_graph_add_meta <- function(g2,
                               pkg_name,
                               deps_metadata){
  igraph::V(g2)$URL <- deps_metadata[names(igraph::V(g2)),]$URL
  igraph::V(g2)$Version <- deps_metadata[names(igraph::V(g2)),]$Version
  igraph::V(g2)$value <- ifelse(names(igraph::V(g2))==pkg_name, 40, 30)
  igraph::V(g2)$group <- ifelse(igraph::V(g2)==pkg_name, 'y', 'n')
  igraph::V(g2)$outputs <- igraph::degree(g2, mode = "out")
  igraph::V(g2)$inputs <- igraph::degree(g2, mode = "in")
  igraph::V(g2)$title <- paste(
    paste0("<strong>Package</strong>: ",
           deps_metadata[names(igraph::V(g2)),]$Package
    ),
    paste0(
      "<strong>GitHub</strong>:<br> ","<a href='",igraph::V(g2)$URL,"'",
      " target='_blank'>",igraph::V(g2)$URL,"</a>"
    ),
    paste0("<strong>Title</strong>: ",
           deps_metadata[names(igraph::V(g2)),]$Title
    ),
    paste0("<strong>Description</strong>: ",
           deps_metadata[names(igraph::V(g2)),]$Description
    ),
    paste0("<strong>Version</strong>: ",
           deps_metadata[names(igraph::V(g2)),]$Version
    ),
    paste0("<strong>Depends</strong>: ",
           deps_metadata[names(igraph::V(g2)),]$Depends
    ),
    paste0("<strong>Imports</strong>: ",
           lapply(deps_metadata[names(igraph::V(g2)),]$Imports,length),
           " packages"
    ),
    paste0("<strong>Suggests</strong>: ",
           lapply(deps_metadata[names(igraph::V(g2)),]$Suggests,length),
           " packages"
    ),
    paste0("<strong>Remotes</strong>: ",
           lapply(deps_metadata[names(igraph::V(g2)),]$Remotes,length),
           " packages"
    ),
    paste0("<strong>SystemRequirements</strong>: ",
           deps_metadata[names(igraph::V(g2)),]$SystemRequirements
    ),
    sep="<br>"
  )
  return(g2)
}
