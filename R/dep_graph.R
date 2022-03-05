#' Dependency graph plot
#'
#' Create a dependency graph between a set of R packages and plot them
#' as an interactive network. By default,
#' plots only packages within the
#'  \href{https://github.com/topics/echoverse}{\code{echoverse}}.
#' @param pkg_name Package to search dependencies for.
#' @param deps A subset of of the main package's (\code{pkg_name} )
#' dependencies to include in the plot visualization.
#' @param show_plot Whether to print the plot.
#' @param save_path Path to save the plot to, as an interactive,
#'  self-container HTML file.
#' @param verbose Print messages.
#' @inheritParams visNetwork::visNodes
#' @inheritParams visNetwork::visSave
#' @export
#' @importFrom dplyr %>%
#' @examples
#' vis <- echolocatoR::dep_graph()
dep_graph <- function(pkg_name="echolocatoR",
                      deps = "echoverse",
                      shape = c("image", "hexagon"),
                      image =
                        file.path(
                          "https://github.com/RajLabMSSM",
                          "Fine_Mapping/blob/master/echolocatoR",
                          "images/bat_silhouette.png?raw=true"),
                      show_plot = TRUE,
                      save_path = file.path(tempdir(),
                                            paste0(pkg_name,".dep_graph.html")),
                      background = "#25355c",
                      verbose = TRUE){
  requireNamespace("pkgnet")
  requireNamespace("igraph")

  ### Gather metadata ####
  deps_metadata <- package_metadata(pkgs = deps)
  all_pkgs <- unique(c(pkg_name,unlist(deps_metadata$Package)))
  #### Gather dependency graph data ####
  report <- pkgnet::CreatePackageReport(pkg_name = pkg_name)
  g <- report$DependencyReporter$pkg_graph$igraph
  #### Subset to only echoverse modules ####
  all_pkgs <- all_pkgs[all_pkgs %in% names(igraph::V(g))]
  g2 <- igraph::induced_subgraph(
    graph =  g,
    vids = igraph::V(g)[all_pkgs]
  )
  #### Add graph metadata ####
  g2 <- dep_graph_add_meta(g2 = g2,
                           pkg_name = pkg_name,
                           deps_metadata = deps_metadata)
  #### Create interactive plot ####
  vis <- dep_graph_plot(g2 = g2,
                        shape = shape,
                        image = image,
                        pkg_name = pkg_name,
                        show_plot = show_plot,
                        save_path = save_path,
                        background = background)
  return(vis)
}
