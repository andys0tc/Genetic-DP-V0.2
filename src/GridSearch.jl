include("GeneticAlgorithms.jl")

import .GeneticAlgorithms as ga

using LightGraphs, DataFrames,CSV, StatsPlots,Statistics
println("Antes de funcion")
function GridNaiveGAVertexCover(tree_sizes::Vector{Int64},max_weight,no_trials)
    #create DataFrame
    #Columns:
    #   tree_sizes
    #   Value Min VertexCove GenDP
    #   pm
    #   Num populations
    #   Num generations
    #thread
    #lk = ReentrantLock()
    #df = DataFrame(TreeSize=Int64[], MinVertexCover = Float64[], pm=Float64[], NumPopu = Float64[],NumGeneration = Float64[], NaiveGA = Float64[],HalfGenDP = Float64[], GenDP=Float64[])
    df = DataFrame(TreeSize=Int64[], MinVertexCover = Float64[], pm=Float64[], NumPopu = Float64[],NumGeneration = Float64[], Value = Float64[], alg=String[])  
    results = zeros(no_trials,3)
    #opt_VC_DP = zeros(size(tree_sizes)[1])
    
    for i in 1:size(tree_sizes)[1]
        #Graph settings
        nodes = tree_sizes[i]
        g = ga.random_tree(nodes)
        w = float(rand(1:max_weight,nv(g)))
        
        #Naive Genetic Algorithm settings
        
        #Probability of mutation
        pm=(range(0.5,3,step = 0.5)|>collect)/nodes
        #pm=[0.1,0.17]
        
        #Population sizes
        #population_sizes = ((range(0.5,3,step = 0.5)|>collect)*sqrt(nodes))
        population_sizes = ((range(0.5,3,step = 0.5)|>collect)*(nodes))
        population_sizes = round.(population_sizes)
        #population_sizes = [16]
        population_sizes = convert(Array{Int64},population_sizes)
        
        #Generation_sizes 
        #generation_sizes = (range(0.5,3,step = 0.5)|>collect)*sqrt(nodes)log(nodes)
        #generation_sizes = round.(generation_sizes)
        generation_sizes=[100]
        generation_sizes = convert(Array{Int64},generation_sizes)##
        #Get opt with Dynammic Programming
        opt_VC_DP = ga.min_vertex_cover_on_trees(g,float(w)) #Programacion Dinamica
        
        #ga_results = zeros(size(pm)[1],size(population_sizes)[1],size(generation_sizes)[1])
        #mwis = ga.MinWeightIndependentSet()
        mwvc = ga.MinWeightVertexCover() #
        #Grid search pm|population_sizes|generation_sizes
        for j in 1:size(pm)[1]
            for k in 1:size(population_sizes)[1]
                popu_size = min(population_sizes[k],100)
                popu = ga.GeneticGraphData{Bool}(popu_size, [true,false],g,w,mwvc)
                #loggerNaive = open("test/graphs/loggerNaiveGA.txt", "w+")


                naive_ga = ga.NaiveGA{Bool}(pm[j], popu)#Naive
                gendp_global = ga.GenDP_Global{Bool}(pm[j], popu,ga.random_cut_tree)#halfGenDP
                gendp = ga.GenDP{Bool}(pm[j], popu,ga.random_cut_tree) #halfGenDP
                
                for l in 1:size(generation_sizes)[1]
                    gen = min(generation_sizes[l],100)
                    #Threads.@threads 
                    for m in 1:no_trials
                    
                        results[m,1] = -1*ga.run(naive_ga, gen)
                        results[m,2] = -1*ga.run(gendp_global, gen)
                        results[m,3] = -1*ga.run(gendp,gen)
                        
                    end

                    nga, lgdp, gdp = mean(results,dims = 1)
                    #push!(df, [nodes opt_VC_DP pm[j] popu_size gen nga lgdp gdp])
                    push!(df, [nodes, opt_VC_DP, pm[j], popu_size, gen, nga, "naive"])
                    push!(df, (nodes, opt_VC_DP, pm[j], popu_size, gen, lgdp, "halfGenDP"))

                    push!(df, (nodes, opt_VC_DP, pm[j], popu_size, gen, gdp, "GenDP"))
                    #df[size_tree,val_opt_VC,pm,N,T]
                    #push!(df, [nodes opt_VC_DP pm[j] popu_size gen nga lgdp gdp])
                end
            end
        end
        println("Nuevo arbol $i")
        
    end
    #df.NaiveGAAppRatio = df.NaiveGA ./ df.MinVertexCover
    #df.HalfGenDPAppRatio = df.HalfGenDP ./ df.MinVertexCover
    #df.GenDPAppRatio = df.GenDP ./ df.MinVertexCover
    return df
end

println("Empezando ...")

#df = GridNaiveGAVertexCover([10],1,2)
df = GridNaiveGAVertexCover([10, 30, 50, 80, 110, 140, 170],1,3)



#CSV.write("/home/andrea/Julia/GenDP-main v02/test/raw_results_VC_13-11-21.csv", df)
CSV.write("C:/Users/atorr/OneDrive - Instituto Politecnico Nacional/Documents/Julia Projects/GenDP-main-GITHUB-Ric-MVC/test/raw_results_VC_23022022.csv",df)
print("Fin de proceso")

# groupby number_nodes
#df_group = groupby(df,:tree_sizes)
#CSV.write("/home/andrea/Julia/GenDP-main/test/raw_results_VC_Group.csv", df_group)

#csv_file = CSV.read("/home/andrea/Julia/GenDP-main/test/raw_results_VC.csv")
#df = DataFrame(csv_file)
#csv_group = CSV.read("/home/andrea/Julia/GenDP-main/test/raw_results_VC_Group.csv")
#dfg = DataFrame(csv_group)

#dfc = combine(dfg,:naiveGA => mean,:MinVertexCover =>mean)


#pyplot()
#mycolours = [:green :blue :orange]
#@df dfc plot(:tree_sizes, :MinVertexCover, group=:tree_sizes, color=mycolours)