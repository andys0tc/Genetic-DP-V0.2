using LightGraphs, Random, Distributions, JLD, DelimitedFiles

mutable struct GenDP{T} <: LocalGA
    pm::Float64
    popu::GeneticGraphData{T}
    partition::Function
    logger::Union{IOStream, Nothing}
    
    function GenDP{T}(pm::Float64, popu::GeneticGraphData{T}, partition::Function) where {T}
        return new{T}(pm, popu, partition, nothing)
    end

    function GenDP{T}(pm::Float64, popu::GeneticGraphData{T}, partition::Function, logger::IOStream) where {T}
        return new{T}(pm, popu, partition, logger)
    end
end

mutable struct GenDP_Global{T} <: GlobalGA
    pm::Float64
    popu::GeneticGraphData{T}
    partition::Function
    logger::Union{IOStream, Nothing}
    
    function GenDP_Global{T}(pm::Float64, popu::GeneticGraphData{T}, partition::Function) where {T}
        return new{T}(pm, popu, partition, nothing)
    end

    function GenDP_Global{T}(pm::Float64, popu::GeneticGraphData{T}, partition::Function, logger::IOStream) where {T}
        return new{T}(pm, popu, partition, logger)
    end
end

function crossover!(gendp_global::GenDP_Global, p1::Individual{T}, p2::Individual{T}, ch1::Individual{T}, ch2::Individual{T}) where {T}
    popu = gendp_global.popu
    comps = Vector{Vector{Int64}}(undef, 2)
    interface = Vector{Vector{Int64}}(undef, 2)
    gendp_global.partition(popu.g, comps, interface)
    partition_crossover!([p1, p2], comps, ch1)
    partition_crossover!([p2, p1], comps, ch2)
end

function random_cut_tree(t::SimpleGraph{Int64}, comps::Vector{Vector{Int64}}, interface::Vector{Vector{Int64}})

    n = nv(t)
    """
    edges_iter = collect(edges(t))
    num_edges = size(edges_iter, 1)

    re = rand(edges_iter)
    ed_src = src(re)
    ed_dst = dst(re)
    """
    v₀ = rand(1:n)
    walk = randomwalk(t, v₀, n)
    ed_src = walk[end-1]
    ed_dst = walk[end]
    
    rem_edge!(t, ed_src, ed_dst)

    components = connected_components(t)
    comps[1] = components[1]
    comps[2] = components[2]

    interface[1] = Vector{Int64}(undef, 1)
    interface[2] = Vector{Int64}(undef, 1)
    if ed_src ∈ comps[1]
        interface[1][1] = ed_src
        interface[2][1] = ed_dst
    else
        interface[2][1] = ed_src
        interface[1][1] = ed_dst
    end

    """ Check if t is not a tree
    if size(edges) >= n-1
        throw Exception
    end
    """
    add_edge!(t, ed_src, ed_dst)

    return (ed_src, ed_dst)
end

function gendp_get_fitness!(popu::GeneticGraphData{T}, g::SimpleGraph{Int64}, vmap::Vector{Int64}) where {T}
    fit = zeros(size(popu.population, 1))
    for i in 1:size(popu.population, 1)
        fit[i] = fitness!(popu.fit, g, popu.w, popu.population[i], vmap)
    end
    return fit 
end

function gendp_get_fitness!(popu::GeneticGraphData{T}, g::SimpleGraph{Int64}, vmap::Vector{Int64}, frozen::Vector{Int64}, values::Vector{Bool}) where {T}
    fit = zeros(size(popu.population, 1))
    old = Vector{Bool}(undef, size(frozen, 1))

    for i in 1:size(popu.population, 1)
        old = popu.population[i].chromosome[frozen]
        popu.population[i].chromosome[frozen] = values
        fit[i] = fitness!(popu.fit, g, popu.w, popu.population[i], vmap)
        popu.population[i].chromosome[frozen] = old
    end
    return fit 
end

function gendp_get_fitness!(popu::GeneticGraphData{T}, g::SimpleGraph{Int64}, vmap::Vector{Int64}, frozen::Vector{Int64}, values::Vector{Bool},res::Vector{Int64}) where {T}
    fit = zeros(size(popu.population, 1),2)
    #fit_zeros = zeros(size(popu.population, 1))
    old = Vector{Bool}(undef, size(frozen, 1))
    old2 = Vector{Bool}(undef, size(res, 1))
    #flag = Vector{Bool}(undef, size(frozen, 1))
    val_interface = Vector{Bool}(undef, size(res, 1))
    for j in 1:2
        if j == 1
            val_interface[1] = false
        else
            val_interface[1] = true
        end
        for i in 1:size(popu.population, 1)
            old = popu.population[i].chromosome[frozen]
            old2 = popu.population[i].chromosome[res]
            popu.population[i].chromosome[frozen] = values
            popu.population[i].chromosome[res] = val_interface
            fit[i,j] = fitness!(popu.fit, g, popu.w, popu.population[i], vmap)
            popu.population[i].chromosome[frozen] = old
            popu.population[i].chromosome[res] = old2
        end
    end
    
    return fit
end


function new_generation!(local_ga::LocalGA)
    popu = local_ga.popu
    n = 1
    j = 1
    l = nv(popu.g)
    comps = Vector{Vector{Int64}}(undef, 2)
    interface = Vector{Vector{Int64}}(undef, 2)

    while n < size(popu.population, 1)
        # Select parents based on random cut
        edge = local_ga.partition(popu.g, comps, interface)

        dom1, res1 = sample_dominant_recesive(popu, comps, interface)
        #pd_dom = min_vertex_cover_on_trees(comps[1],ones(size(comps[1], 1), ))
        #popu.population[res1].chromosome[res]
        partition_crossover!([popu.population[dom1], popu.population[res1]], comps, popu.new_population[n])
        
        fix_interface!(popu.fit, popu.g, popu.w, popu.new_population[n], interface)
        save_history(popu ,local_ga.logger , n)
        #@debug begin
        
            #save_iteration(popu, local_ga.logger, dom1, res1, n)
        #end
        mutate!(popu.new_population[n], local_ga.pm)
        
        n += 1
        reverse!(comps)
        reverse!(interface)
        dom2, res2 = sample_dominant_recesive(popu, comps, interface)
        partition_crossover!([popu.population[dom2], popu.population[res2]], comps, popu.new_population[n])
        fix_interface!(popu.fit, popu.g, popu.w, popu.new_population[n], interface)
        #@debug begin
        #save_iteration(popu, local_ga.logger, dom2, res2, n)
        
        #end
        #save_iteration(popu, local_ga.logger    , p1, p2, n)

        mutate!(popu.new_population[n], local_ga.pm)
        save_history(popu ,local_ga.logger , n)
        n += 1
        j += 1
    end    

    #historico 

    tmp = popu.population
    popu.population = popu.new_population
    popu.new_population = tmp
end

function sample_dominant_recesive(popu::GeneticGraphData{T}, comps::Vector{Vector{Int64}}, interface::Vector{Vector{Int64}}) where {T}
    t_dom, vmap_dom = induced_subgraph(popu.g, comps[1])
    
    val = Vector{Bool}(undef, size(interface[2], 1))
    #root_dom = findall(x->x==interface[1][1], vmap_dom)
    #pd_dom = min_vertex_cover_on_trees(t_dom,ones(size(comps[1], 1), ),root_dom[1])
    fit_dom = gendp_get_fitness!(popu, t_dom, vmap_dom)
    roulette_dom = create_distribution_exp!(fit_dom)
    dom = rand(roulette_dom)

    t_res, vmap_res = induced_subgraph(popu.g, append!(comps[2], interface[1]))
    frozen = popu.population[dom].chromosome[interface[1]]
    #get best recesive given all the posibles values in interface[2]
    #res_ind = popu.population[dom].chromosome[interface[1]]

    #root_res = findall(x->x==interface[2][1], vmap_res)
    #pd_res = min_vertex_cover_on_trees(t_res,ones(size(comps[2], 1), ),root_res[1])
    fit_res = gendp_get_fitness!(popu, t_res, vmap_res, interface[1], frozen,interface[2])
    roulette_res_zeros = create_distribution_exp!(fit_res[:,1])
    roulette_res_ones = create_distribution_exp!(fit_res[:,2])

    #roulette_res = create_distribution_exp!(fit_res)
    res_zeros = rand(roulette_res_zeros)
    res_ones = rand(roulette_res_ones)

    if fit_res[res_zeros,1]>fit_res[res_ones,2]
        res = res_zeros
        val[1] = false
        
        #flag[1] = false
    else
        res = res_ones
        val[1] = true
        #popu.population[dom].chromosome[interface[2]] = true
        #flag[1] = true
    end
    popu.population[dom].chromosome[interface[2]] = val
    
    for i in 1:size(interface[1], 1)
        pop!(comps[2])
    end

    return (dom, res)
end


function run_gendp(popu::GeneticGraphData{T}, gen::Int64) where {T}
    initialize_ga!(popu)
    max_val = typemin(Float64)
    for i in 1:gen
        fit = get_fitness!(popu)
        max_gen = maximum(fit)
        
        if max_gen >= max_val
            max_val = max_gen
        end
        gendp!(popu)
    end 
    return max_val
end

function save_history(popu::GeneticGraphData, io::IOStream , n::Int64)
    #print(popu.new_population[n].chromosome,"\n")
    #write(io, popu.new_population[n].chromosome,"\n")
    print(io, popu.new_population[n].chromosome,"\n")
    #writedlm(io, popu.new_population[n].chromosome)
end

function save_iteration(popu::GeneticGraphData, io::IOStream, parent1, parent2, child1)
    p1_fitness = fitness!(popu.fit, popu.g, popu.w, popu.population[parent1])
    p2_fitness = fitness!(popu.fit, popu.g, popu.w, popu.population[parent2])
    ch1_fitness = fitness!(popu.fit, popu.g, popu.w, popu.new_population[child1])
    print("$p1_fitness $p2_fitness $ch1_fitness\n")
    print(io, "$p1_fitness $p2_fitness $ch1_fitness\n")
    # println(io, popu.population[parent1])
    # println(io, popu.population[parent2])
    # println(io, popu.population[child1])
    # println(io, popu.population[child2])
end

include("FitnessFunctions.jl")