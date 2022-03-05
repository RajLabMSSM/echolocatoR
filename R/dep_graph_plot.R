dep_graph_plot <- function(g2,
                           pkg_name,
                           shape = c("image", "hexagon"),
                           image =
                             file.path(
                               "https://github.com/RajLabMSSM",
                               "Fine_Mapping/blob/master/echolocatoR",
                               "images/bat_silhouette.png?raw=true"),
                           show_plot = TRUE,
                           save_path = NULL,
                           width = NULL,
                           height = NULL,
                           background = "white"){
  requireNamespace("visNetwork")

  vis <- visNetwork::visIgraph(g2) %>%
    visNetwork::visIgraphLayout(layout = "layout_as_star",
                                center=pkg_name,
                                randomSeed = 11) %>%
    visNetwork::visNodes(
      shape = tolower(shape[1]),
      borderWidth = 2,
      image = image,
      labelHighlightBold = TRUE,
      color = list(
        background =  "#25355c",
        border = "#41c6c8",
        highlight = "#56ffff"
      ),
      font = list(color="white",
                  size=20,
                  face="Aller_Rg",
                  strokeWidth=10,
                  strokeColor="#25355c"),
      shadow = list(enabled = TRUE,
                    size = 40,
                    color="#537bcb") # "#03b1f0"
    ) %>%
    visNetwork::visEdges(
      arrows = 'from',
      shadow = list(enabled=TRUE,color="#686ea6"),
      smooth = TRUE,dashes =FALSE,
      width = 2,
      color = list(color = "#56ffff",
                   opacity=.75,
                   highlight = "#686ea6"),
    ) %>%
    visNetwork::visOptions(nodesIdSelection = list(enabled = FALSE,
                                                   selected=pkg_name,
                                                   main="select package"),
                           highlightNearest=TRUE,
                           width = width,
                           height = height) %>%
    visNetwork::visLayout(randomSeed = 11)

  if(show_plot) print(vis)

  if(!is.null(save_path)) {
    message("Saving dependency graph plot ==> ",save_path)
    visNetwork::visSave(graph = vis,
                        file = save_path,
                        background = background,
                        selfcontained = TRUE)
  }
  return(vis)
}
