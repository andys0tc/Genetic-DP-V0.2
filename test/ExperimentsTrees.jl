include("../src/GeneticAlgorithms.jl")
import .GeneticAlgorithms as ga

using LightGraphs, DataFrames, CSV, Statistics

#, 110, 140, 170, 200

avg_m = [0.5, 1.0, 1.5]
pop_size_factor = [0.5, 1.0, 2.0]
gens_factor = [1.0, 2.0, 3.0]

df = ga.random_tree_experiments([200, 250, 300], avg_m, pop_size_factor, gens_factor)
CSV.write("raw_results.csv", df)
df_group = groupby(df, [:tree_size, :instance, :p_mutation, :population, :generations])

dfc = combine(df_group, :GAAppRatio => mean, :GAAppRatio => std, :GenDPHalfAppRatio => mean, :GenDPHalfAppRatio => std, :GenDPAppRatio => mean, :GenDPAppRatio => std)
CSV.write("results_stats.csv", dfc)
print(dfc)
