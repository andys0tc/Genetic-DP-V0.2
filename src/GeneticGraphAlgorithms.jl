using LightGraphs, Random, Distributions, StatsBase



mutable struct GeneticGraphData{T}
    population::Vector{Individual{T}}
    new_population::Vector{Individual{T}}
    g::SimpleGraph{Int64}
    w::Vector{Float64}
    domain::Vector{T}
    fit::FitnessFunction

    function GeneticGraphData{T}(N::Int64, domain::Vector{T}, g::SimpleGraph{Int64}, w::Vector{Float64}, fit::FitnessFunction) where {T}
        population = Vector{Individual{T}}(undef, N)
        new_population = Vector{Individual{T}}(undef, N)

        for i in 1:N
            population[i] = Individual{T}()
            new_population[i] = Individual{T}()
        end
        return new{T}(population, new_population, g, w, domain, fit)
    end
end

mutable struct NaiveGA{T} <: GlobalGA
    pm::Float64
    popu::GeneticGraphData{T}
    logger::Union{IOStream, Nothing}
    
    function NaiveGA{T}(pm::Float64, popu::GeneticGraphData{T}) where {T}
        return new{T}(pm, popu, nothing)
    end

    function NaiveGA{T}(pm::Float64, popu::GeneticGraphData{T}, logger::IOStream) where {T}
        return new{T}(pm, popu, logger)
    end
end

function initialize!(popu::GeneticGraphData{T}) where {T}
    for i in 1:size(popu.population, 1)
        popu.population[i].domain = popu.domain
        popu.new_population[i].domain = popu.domain
        rand_init!(popu.population[i], nv(popu.g))
        rand_init!(popu.new_population[i], nv(popu.g))
    end
end

function get_fitness!(popu::GeneticGraphData{T}) where {T}
    fitness_values = zeros(size(popu.population, 1))
    for i in 1:size(popu.population, 1)
        fitness_values[i] = fitness!(popu.fit, popu.g, popu.w, popu.population[i])
    end
    return fitness_values
end

function create_distribution_linear!(values::Vector{Float64})
    min_v = minimum(values)
    max_v = maximum(values)

    for i in 1:size(values, 1)
        values[i] = values[i]-min_v + 1.0
    end
    s = sum(values)
    for i in 1:size(values, 1)
        values[i] = values[i]/s
    end

    return Categorical(values)
end

function create_distribution_exp!(values::Vector{Float64}, alpha=1.0)
    
    #values = scaleVector!(values)
    for i in 1:size(values, 1)
        values[i] = exp(alpha*values[i])+floatmin()
    end
    s = sum(values)
    for i in 1:size(values, 1)
        values[i] = values[i]/s
    end

    return Categorical(values)
end

function scaleVector!(values::Vector{Float64})
    
    min = minimum(values)
    max = maximum(values)
    values = (values.-min)./(max-min)
    return values
end

function getDist(values::Vector{Float64},num_exp::Int64)
    count = zeros(size(values,1))
    round = Categorical(values)
    for i in 1:num_exp
        r = rand(round)
        count[r] = count[r]+1
    end
       
    return count
    
end

function tournament!(popu::GeneticGraphData, fitness_values::Vector{Float64}, alpha=1.0)     
    
    
    for i in 1:size(popu.population, 1)
        #Select 2 individuals
        selected = rand(1:size(popu.population, 1),2)
        #Exponentiate the fitness values
    
        f_i = exp(alpha*fitness_values[selected[1]])
        f_j = exp(alpha*fitness_values[selected[2]])
        winner = sample(selected, Weights([f_i, f_j]))
        popu.new_population[i] = popu.population[winner] 
    end
    

end
    

function new_generation!(global_ga::GlobalGA) 
    popu = global_ga.popu
    fitness_values = get_fitness!(popu)
    ft=copy(fitness_values)
    roulette = create_distribution_exp!(fitness_values)
    # roulette = create_distribution_linear!(fit)
    n = 1
    N = size(popu.population, 1)
    while n < N
        p1 = rand(roulette) # naive selection
        p2 = rand(roulette) # naive selection

        crossover!(global_ga, popu.population[p1], popu.population[p2], popu.new_population[n], popu.new_population[n+1])

        @debug begin
            save_iteration(popu, global_ga.logger, p1, p2, n)
            save_iteration(popu, global_ga.logger, p1, p2, n+1)
        end

        mutate!(popu.new_population[n], global_ga.pm)
        mutate!(popu.new_population[n+1], global_ga.pm)
        n += 2

    end

    tmp = popu.population
    popu.population = popu.new_population
    popu.new_population = tmp
end

function new_generation!(global_ga::GlobalGA,alpha=0.3) 
    popu = global_ga.popu
    fitness_values = get_fitness!(popu)
    ft=copy(fitness_values)

    n = 1
    N = size(popu.population, 1)
    roulette = create_distribution_exp!(fitness_values)
    #ex = getDist(roulette.p,1000)
    #tournament!(popu, fitness_values, alpha)
    while n < N
       
        #p1, p2 = rand(1:size(popu.population, 1),2)
        p1 = rand(roulette) # naive selection
        p2 = rand(roulette) # naive selection

        #after tournament the new values are in popu.new_population
        #parents are in popu.new_population and kids are popu.population
        #crossover!(global_ga, popu.new_population[p1], popu.new_population[p2], popu.population[n], popu.population[n+1])
        crossover!(global_ga, popu.population[p1], popu.population[p2], popu.new_population[n], popu.new_population[n+1])
        
        mutate!(popu.new_population[n], global_ga.pm)
        mutate!(popu.new_population[n+1], global_ga.pm)
        save_history(popu, global_ga.logger,n)
        save_history(popu, global_ga.logger,n+1)
        n += 2
    end

    tmp = popu.population
    popu.population = popu.new_population
    popu.new_population = tmp
end
function save_history(popu::GeneticGraphData, io::IOStream , n::Int64)

    print(io, ga.popu.population[n].chromosome)
end

function save_iteration(popu::GeneticGraphData, io::IOStream, parent1, parent2, child1, child2)
    p1_fitness = fitness!(popu.fit, popu.g, popu.w, popu.population[parent1])
    p2_fitness = fitness!(popu.fit, popu.g, popu.w, popu.population[parent2])
    ch1_fitness = fitness!(popu.fit, popu.g, popu.w, popu.new_population[child1])
    ch2_fitness = fitness!(popu.fit, popu.g, popu.w, popu.new_population[child2])
    print("$p1_fitness $p2_fitness $ch1_fitness $ch2_fitness\n")
    print(io, "$p1_fitness $p2_fitness $ch1_fitness $ch2_fitness\n")
    # println(io, popu.population[parent1])
    # println(io, popu.population[parent2])
    # println(io, popu.population[child1])
    # println(io, popu.population[child2])
end

function run(ga::GeneticAlgorithm, gen::Int64)

    initialize!(ga.popu)
    max_val = -floatmax(Float64)
    #record = zeros(gen)
    for i in 1:gen
        fit = get_fitness!(ga.popu)
        max_gen = maximum(fit)
        #print("$max_gen ")
        if max_gen >= max_val
            max_val = max_gen
        end
        new_generation!(ga)
    end

    return max_val
end

function crossover!(naive_ga::NaiveGA, p1::Individual{T}, p2::Individual{T}, ch1::Individual{T}, ch2::Individual{T}) where {T}
    naive_crossover!(p1, p2, ch1, ch2)
end




include("GeneticDPonTrees.jl")