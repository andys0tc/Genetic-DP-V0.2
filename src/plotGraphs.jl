#using LightGraphs 
using GraphPlot
using Graphs

using Cairo, Compose
#g = SimpleGraph(3)
g = loadgraph("C:/Users/atorr/OneDrive - Instituto Politecnico Nacional/Documents/Julia Projects/Genetic Dinamic Programming/test/graphs/gtest80.lgz")
#g = SimpleGraph(10)
#add_edge!(g,1,2);
#add_edge!(g,1,3);
#add_edge!(g,2,4);
#add_edge!(g,2,5);
#add_edge!(g,5,7);
#add_edge!(g,5,8);
#add_edge!(g,3,6);
#add_edge!(g,6,9);
#add_edge!(g,6,10);

nodelabel = 1:nv(g)
gplot(g)
draw(PNG("tree80.png", 30cm, 30cm), gplot(g, nodelabel=nodelabel))
#plot()