makerandomgraph <- function(numvertices,numedges)
  {
     require("graph")
     #vector with the names of the vertices
     mynames <- sapply(seq(1,numvertices),toString)
     myrandomgraph <- randomEGraph(mynames, edges = numedges)
     return(myrandomgraph)
  }

myrandomgraph <- makerandomgraph(15, 43)

myrandomgraphplot <- layoutGraph(myrandomgraph, layoutType="neato");renderGraph(myrandomgraphplot)
