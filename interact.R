require("graph");library("Rgraphviz");library("RBGL")

mytable <- read.table("binary_interaction_data.txt") # Store the data in a data frame

proteins1 <- mytable$V1

proteins2 <- mytable$V2

protnames <- c(levels(proteins1),levels(proteins2))

# Find out how many pairs of proteins there are

numpairs <- length(proteins1)

# Find the unique protein names:
uniquenames <-  unique(protnames)
     
# Make a graph for these proteins with no edges:
mygraph <- new("graphNEL", nodes = uniquenames)

# Add edges to the graph:See http://rss.acs.unt.edu/Rdoc/library/graph/doc/graph.pdf for more examples

weights <- rep(1,numpairs)

mygraph2 <- addEdge(as.vector(proteins1),as.vector(proteins2),mygraph,weights) 

##mygraphplot <- layoutGraph(mygraph2, layoutType="neato");renderGraph(mygraphplot)

mynodes <- nodes(mygraph2) 

mydegrees <- graph::degree(mygraph2);mydegrees<-sort(mydegrees)

Tbledegrees<- table(mydegrees) ;hist(mydegrees, col="red")

myconnectedcomponents <- connectedComp(mygraph2)

componentsizes <- numeric(length(myconnectedcomponents))

for (i in 1:length(myconnectedcomponents)) {
   component <- myconnectedcomponents[[i]] # Store the connected component in a vector "component"
   componentsize <- length(component)      # Find the number of vertices in this connected component
   componentsizes[i] <- componentsize      # Store the size of this component
}

findcomponent <- function(graph,vertex)
  {
     # Function to find the connected component that contains a particular vertex
     require("RBGL")
     found <- 0
     myconnectedcomponents <- connectedComp(graph)
     numconnectedcomponents <- length(myconnectedcomponents)
     for (i in 1:numconnectedcomponents)
     {
        componenti <- myconnectedcomponents[[i]]
        numvertices <- length(componenti)
        for (j in 1:numvertices)
        {
           vertexj <- componenti[j]
           if (vertexj == vertex)
           {
              found <- 1
              return(componenti)
           }
        }
     }
     print("ERROR: did not find vertex in the graph")
  }

mycomponent <- findcomponent(mygraph2, "YBR009C")

component3 <- myconnectedcomponents[[3]]

mysubgraph <- subGraph(component3, mygraph2) ;
mygraphplot <- layoutGraph(mysubgraph, layoutType="neato");renderGraph(mygraphplot)

findcommunities <- function(mygraph,minsize)
  {
     # Function to find network communities in a graph
     # Load up the igraph library:
     require("igraph")
     # Set the counter for the number of communities:
     cnt <- 0
     # First find the connected components in the graph:
     myconnectedcomponents <- connectedComp(mygraph)
     # For each connected component, find the communities within that connected component:
     numconnectedcomponents <- length(myconnectedcomponents)
     for (i in 1:numconnectedcomponents)
     {
        component <- myconnectedcomponents[[i]]
        # Find the number of nodes in this connected component:
        numnodes <- length(component)
        if (numnodes > 1) # We can only find communities if there is more than one node
        {
           mysubgraph <- subGraph(component, mygraph)
           # Find the communities within this connected component:
           # print(component)
           myvector <- vector()
           mylist <- findcommunities2(mysubgraph,cnt,"FALSE",myvector,minsize)
           cnt <- mylist[[1]]
           myvector <- mylist[[2]]
        }
     }
     print(paste("There were",cnt,"communities in the input graph"))
  }
findcommunities2 <- function(mygraph,cnt,plot,myvector,minsize)
  {
     # Function to find network communities in a connected component of a graph
     # Find the number of nodes in the input graph
     nodes <- nodes(mygraph)
     numnodes <- length(nodes)
     # Record the vertex number for each vertex name
     myvector <- vector()
     for (i in 1:numnodes)
     {
        node <- nodes[i] # "node" is the vertex name, i is the vertex number
        myvector[`node`] <- i  # Add named element to myvector
     }
     # Create a graph in the "igraph" library format, with numnodes nodes:
     newgraph <- graph.empty(n=numnodes,directed=FALSE)
     # First record which edges we have seen already in the "mymatrix" matrix,
     # so that we don't add any edge twice:
     mymatrix <- matrix(nrow=numnodes,ncol=numnodes)
     for (i in 1:numnodes)
     {
        for (j in 1:numnodes)
        {
           mymatrix[i,j] = 0
           mymatrix[j,i] = 0
        }
     }
     # Now add edges to the graph "newgraph":
     for (i in 1:numnodes)
     {
        node <- nodes[i] # "node" is the vertex name, i is the vertex number
        # Find the nodes that this node is joined to:
        neighbours <- adj(mygraph, node)
        neighbours <- neighbours[[1]] # Get the list of neighbours
        numneighbours <- length(neighbours)
        if (numneighbours >= 1) # If this node "node" has some edges to other nodes
        {
           for (j in 1:numneighbours)
           {
              neighbour <- neighbours[j]
              # Get the vertex number
              neighbourindex <- myvector[neighbour]
              neighbourindex <- neighbourindex[[1]]
              # Add a node in the new graph "newgraph" between vertices i and neighbourindex
              # In graph "newgraph", the vertices are counted from 0 upwards.
              indexi <- i
              indexj <- neighbourindex
              # If we have not seen this edge already:
              if (mymatrix[indexi,indexj] == 0 && mymatrix[indexj,indexi] == 0)
              {
                 mymatrix[indexi,indexj] <- 1
                 mymatrix[indexj,indexi] <- 1
                 # Add edges to the graph "newgraph"
                 newgraph <- add.edges(newgraph, c(i, neighbourindex))
              }
           }
        }
     }
     # Set the names of the vertices in graph "newgraph":
     newgraph <- set.vertex.attribute(newgraph, "name", value=nodes)
     # Now find communities in the graph:
     communities <- spinglass.community(newgraph)
     # Find how many communities there are:
     sizecommunities <- communities$csize
     numcommunities <- length(sizecommunities)
     # Find which vertices belong to which communities:
     membership <- communities$membership
     # Get the names of vertices in the graph "newgraph":
     vertexnames <- get.vertex.attribute(newgraph, "name")
     # Print out the vertices belonging to each community:
     for (i in 1:numcommunities)
     {
        cnt <- cnt + 1
        nummembers <- 0
        printout <- paste("Community",cnt,":")
        for (j in 1:length(membership))
        {
           community <- membership[j]
           if (community == i) # If vertex j belongs to the ith community
           {
              vertexname <- vertexnames[j]
              if (plot == FALSE)
              {
                 nummembers <- nummembers + 1
                 # Print out the vertices belonging to the community
                 printout <- paste(printout,vertexname)
              }
              else
              {
                 # Colour in the vertices belonging to the community
                 myvector[`vertexname`] <- cnt
              }
           }
         }
         if (plot == FALSE && nummembers >= minsize)
         {
            print(printout)
         }
      }
      return(list(cnt,myvector))
   }

findcommunities(mysubgraph, 1)
plotcommunities <- function(mygraph)
  {
     # Function to plot network communities in a graph
     # Load the "igraph" package:
     require("igraph")
     # Make a plot of the graph
     graphplot <- layoutGraph(mygraph, layoutType="neato")
     renderGraph(graphplot)
     # Get the names of the nodes in the graph:
     vertices <- nodes(mygraph)
     numvertices <- length(vertices)
     # Now record the colour of each vertex in a vector "myvector":
     myvector <- vector()
     colour <- "red"
     for (i in 1:numvertices)
     {
        vertex <- vertices[i]
        myvector[`vertex`] <- colour   # Add named element to myvector
     }
     # Set the counter for the number of communities:
     cnt <- 0
     # First find the connected components in the graph:
     myconnectedcomponents <- connectedComp(mygraph)
     # For each connected component, find the communities within that connected component:
     numconnectedcomponents <- length(myconnectedcomponents)
     for (i in 1:numconnectedcomponents)
     {
        component <- myconnectedcomponents[[i]]
        # Find the number of nodes in this connected component:
        numnodes <- length(component)
        if (numnodes > 1) # We can only find communities if there is more than one node
        {
           mysubgraph <- subGraph(component, mygraph)
           # Find the communities within this connected component:
           mylist <- findcommunities2(mysubgraph,cnt,"TRUE",myvector,0)
           cnt <- mylist[[1]]
           myvector <- mylist[[2]]
        }
      }
      # Get a set of cnt colours, where cnt is equal to the number of communities found:
      mycolours <- rainbow(cnt)
      # Set the colour of the vertices, so that vertices in each community are of the same colour,
      # and vertices in different communities are different colours:
      myvector2 <- vector()
      for (i in 1:numvertices)
      {
         vertex <- vertices[i]
         community <- myvector[vertex]
         mycolour <- mycolours[community]
         myvector2[`vertex`] <- mycolour
     }
     nodeRenderInfo(graphplot) = list(fill=myvector2)
     renderGraph(graphplot)
 }
plotcommunities(mysubgraph)



