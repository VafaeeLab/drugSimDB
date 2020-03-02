library(igraph)
i2d.db = read.delim("../data/external/i2d.2_9.Public.HUMAN.tab", sep="\t", header = T)
drugRepo.db = read.delim("../data/new_net_info_V3_pval.csv", sep=",", header = T)
drugTarget.db = read.delim("../data/external/uniprot links.csv", sep=",", header = T)
drugTarget.db = drugTarget.db[which(drugTarget.db$Type == "SmallMoleculeDrug"),]

ppi.net.all = simplify(graph.data.frame(i2d.db[1:100,2:3], directed = F))

# Testing zone: START
v = as.vector(drugTarget.db$UniProt.ID[which(drugTarget.db$DrugBank.ID == "DB00014")])
s.paths = get.all.shortest.paths(ppi.net, unique(v), to = unique(v), mode = "all")
g = make_empty_graph() + vertices(s.paths) + path(s.paths$res) + path(s.paths$res[[2]]) + + path(s.paths$res[[3]]) + + path(s.paths$res[[4]]) + + path(s.paths$res[[5]]) + + path(s.paths$res[[6]])
tkplot(make_empty_graph() + path(s.paths))

# Create a subgraph from the neighborhood of specific vertices in igraph
## small example (i.e. part of the main PPI net)
v =  c("P63104", "Q9H2F3", "A0AV96", "A0AVT1", "Q15717")
g = subgraph(ppi.net.all, v)
## try on the whole PPI
g = ppi.net.all

V(g)$name <- as_ids(V(g))
nodes_of_interest = c("P63104", "Q9H2F3", "A0AVT1")
# nodes_of_interest = c("P63104", "Q9H2F3", "Q15717")
selnodes <- V(g)[name %in% nodes_of_interest]
selegoV <- ego(g, order=3, nodes = selnodes, mode = "all", mindist = 0)
selegoG <- induced_subgraph(g,unlist(selegoV))
allV <- as.character(V(selegoG)$name)
plot(selegoG)
list_of_edges <- E(selegoG)[from(nodes_of_interest) | to(nodes_of_interest)]
your_subgraph <- subgraph.edges(selegoG, list_of_edges)
tkplot(your_subgraph)

# ego_graph = make_ego_graph(g, order = 2, nodes = selnodes, mode = "all", mindist = 1)


## try on the whole PPI

all.sp = all_shortest_paths(small.net, from = V(small.net), to = V(small.net))

list_of_edges <- E(ppi.net.new)[from(v) | to(v)]
your_subgraph <- subgraph.edges(ppi.net.new, list_of_edges)
# Testing zone: END
