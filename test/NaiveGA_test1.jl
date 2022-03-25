include("../src/GeneticAlgorithms.jl")
import .GeneticAlgorithms as ga

using LightGraphs, DataFrames, DelimitedFiles

#g = ga.random_tree(300)
#g = ga.random_tree(300)
#g = ga.rand_binary_tree(30)
g = ga.random_tree(80)
savegraph("C:/Users/atorr/OneDrive - Instituto Politecnico Nacional/Documents/Julia Projects/Genetic Dinamic Programming/test/graphs/gtest80.lgz", g)
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


w = ones(nv(g))
#w = rand(nv(g))

#opt = ga.max_ind_set_on_trees(g, w)
#mwis = ga.MaxWeightIndependentSet(penalty=1.0)
#print("Max Independent set is $opt \n")


#N = Int64(floor(3*nv(g)))
#N = 10
#popu = ga.GeneticGraphData{Bool}(N, [true, false], g, w, mwis)
#pm = 0.5/nv(g)
##
##  loggerNaive = open("test/graphs/loggerNaiveGA.txt", "w+")
##  print(loggerNaive, "Parent1 Parent2 Child1\n")
##  naive_ga = ga.NaiveGA{Bool}(pm, popu, loggerNaive)
##  print("Running Naive GA... \n")
##  nga = ga.run(naive_ga, 5)
##  print("$nga \n")
##  close(loggerNaive)
##  ##
##  loggerGen = open("test/graphs/loggerGlobalGenDP.txt", "w+")
##  print(loggerGen, "Parent1 Parent2 Child1\n")
##  gendp_global = ga.GenDP_Global{Bool}(pm, popu, ga.random_cut_tree, loggerGen)
##  print("Running local GenDP... \n")
##  lgdp = ga.run(gendp_global, 5)
##  print("GenDp:$lgdp Naive:$nga")
##  close(loggerGen)
##  ##
##  loggerGenDP = open("test/graphs/loggerLocalGenDP.txt", "w+")
##  print(loggerGenDP, "Parent1 Parent2 Child1\n")
##  gendpL = ga.GenDP{Bool}(pm, popu, ga.random_cut_tree, loggerGenDP)
##  print("\nRunning GenDP... \n")
##  gdp = ga.run(gendpL, 5)
##  print("GenDp:$lgdp Naive:$nga GenDpLocal:$gdp\n")
##  close(loggerGenDP)
##  


###Min Vertex Cover###

opt_mvc = ga.min_vertex_cover_on_trees(g, w,1)

mvc = ga.MinWeightVertexCover(penalty=1.0)
print("Min Vertex cover is $opt_mvc \n")


N = Int64(floor(3*nv(g)))
#N = 30
popu = ga.GeneticGraphData{Bool}(N, [true, false], g, w, mvc)
pm = 0.5/nv(g)
#loggerNaive = open("C:/Users/atorr/OneDrive - Instituto Politecnico Nacional/Documents/Julia Projects/Genetic Dinamic Programming V0.2/test/history/NaiveGA80.txt", "w+")
#
##print(loggerNaive, "Parent1 Parent2 Child1\n")
#naive_ga = ga.NaiveGA{Bool}(pm, popu, loggerNaive)
#print("Running Naive GA... \n")
#nga = -1*ga.run(naive_ga, 100)
#print("$nga \n")
#close(loggerNaive)
#
#
#loggerGen = open("C:/Users/atorr/OneDrive - Instituto Politecnico Nacional/Documents/Julia Projects/Genetic Dinamic Programming V0.2/test/history/GlobalGenDP80.txt", "w+")
##print(loggerGen, "Parent1 Parent2 Child1\n")
#gendp_global = ga.GenDP_Global{Bool}(pm, popu, ga.random_cut_tree, loggerGen)
##gendp_global = ga.GenDP_Global{Bool}(pm, popu, ga.random_cut_tree)
#print("Running local GenDP... \n")
#lgdp = -1*ga.run(gendp_global, 100)
#print("GenDp:$lgdp Naive:$nga")
#close(loggerGen)

loggerGenDP = open("C:/Users/atorr/OneDrive - Instituto Politecnico Nacional/Documents/Julia Projects/Genetic Dinamic Programming V0.2/test/history/GenDP80.txt", "w+")
#history_file = 
#print(loggerGenDP, "Parent1 Parent2 Child1\n")
gendpL = ga.GenDP{Bool}(pm, popu, ga.random_cut_tree, loggerGenDP)
#gendpL = ga.GenDP{Bool}(pm, popu, ga.random_cut_tree)
print("\nRunning GenDP... \n")
gdp = -1*ga.run(gendpL, 100)
print("GenDp:$gdp \n")
#print("GenDp:$lgdp Naive:$nga GenDpLocal:$gdp\n")
close(loggerGenDP)