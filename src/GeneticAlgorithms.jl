module GeneticAlgorithms

export Individual
export GeneticGraphData

mutable struct Individual{T}
    chromosome::Union{Vector{T}, Nothing}
    domain::Union{Vector{T}, Nothing}
    fitness::Float64

    function Individual{T}() where {T}
        chromosome = nothing
        domain = nothing
        fitness = 0.0
        return new(chromosome, domain, fitness)
    end
end

"""
    rand_init!(ind, n)

Create a new chromosome (array) of size `n` for individual `ind` initialized with random values from the `domain` of the individual.
"""
function rand_init!(ind::Individual{T}, n::Int64) where {T}
    ind.chromosome = rand(ind.domain, n)
end

"""
    mutate!(ind, pm::Float64)

With probability `pm` change independently the value of each entry in chromosome to a random value in the `domain`.
"""
function mutate!(ind::Individual{T}, pm::Float64) where {T}
    for i in 1:size(ind.chromosome, 1)
        if rand() <= pm
            #ind.chromosome[i] = rand(ind.domain)
            ind.chromosome[i] = !ind.chromosome[i]
        end        
    end 
end

"""
    naive_crossover!(p1, p2, ch1, ch2)

Generate new children using the standard crossover operation defined by Genetic Algorithms.

Pick random index `cp` between 2 and ``n-1`` and generate new children by taking the first subsequence `chromosome[1:cp]` of chromosome form parent `p1` and the second subsequence `chromosome[cp+1:end]` form the second parent `p2` and viceversa.
"""
function naive_crossover!(p1::Individual{T}, p2::Individual{T}, ch1::Individual{T}, ch2::Individual{T}) where {T}   
    n = size(p1.chromosome, 1) 
    cp = rand(2:n-1)

    ch1.chromosome[1:cp] = p1.chromosome[1:cp]
    ch1.chromosome[cp+1:end] = p2.chromosome[cp+1:end]

    ch2.chromosome[1:cp] = p2.chromosome[1:cp]
    ch2.chromosome[cp+1:end] = p1.chromosome[cp+1:end]
end

function update_parents!(p1::Individual{T}, p2::Individual{T}, ch1::Individual{T}, ch2::Individual{T}) where {T}
    p1.chromosome = ch1.chromosome
    p2.chromosome = ch2.chromosome
end

"""
    partition_crossover!(parents, partition, child)

Create a new `child` by taking elements from the `parents[i]` specified by the `partition[i]`
"""
function partition_crossover!(parents::Vector{Individual{T}}, partition::Vector{Vector{Int64}},child::Individual{T}) where {T}   
    
    child.chromosome[partition[1]] = parents[1].chromosome[partition[1]]
    child.chromosome[partition[2]] = parents[2].chromosome[partition[2]]
end

abstract type GeneticAlgorithm end

"""
    GlobalGA

Genetic Algorithms that select individuals using global fitness values.
"""
abstract type GlobalGA <: GeneticAlgorithm end

"""
    LocalGA

Genetic Algorithms that select individuals using local fittness values.
"""
abstract type LocalGA <: GeneticAlgorithm end

include("FitnessFunctions.jl")
include("GeneticGraphAlgorithms.jl")
include("GeneticDPonTrees.jl")
include("GraphGenerators.jl")
include("GenDPExperiments.jl")

end