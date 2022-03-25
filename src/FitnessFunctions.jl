using LightGraphs


abstract type FitnessFunction end

struct MaxWeightIndependentSet <: FitnessFunction 
    penalty::Float64
    function MaxWeightIndependentSet(;penalty=1.0)
        return new(penalty)
    end
end

struct MinWeightVertexCover <: FitnessFunction 
    penalty::Float64
    function MinWeightVertexCover(;penalty=1.0)
        return new(penalty)
    end
end


function fitness_max_ind_set!(g::SimpleGraph{Int64}, w::Vector{Float64}, ind::Individual{Bool}, vmap::Vector{Int64}, penalty::Float64=10.0)
    n = nv(g)
    fitness = 0.0
    conflict = 0.0

    for i in 1:n
        if ind.chromosome[vmap[i]] == true
            neighbors = outneighbors(g, i)
            for node in neighbors
                if ind.chromosome[vmap[node]] == true
                    conflict += w[vmap[i]]
                    conflict += w[node]
                end
            end

            if isapprox(conflict, 0.0)
                fitness += w[vmap[i]]
            else
                fitness -= (penalty * conflict)/2
            end
            conflict = 0.0
        end
    end
    ind.fitness = fitness
    return fitness
end

function fitness!(mwis::MaxWeightIndependentSet, g::SimpleGraph{Int64}, w::Vector{Float64}, ind::Individual{Bool}, vmap::Vector{Int64})
    return fitness_max_ind_set!(g, w, ind, vmap, mwis.penalty)
end

function fitness_max_ind_set!(g::SimpleGraph{Int64}, w::Vector{Float64}, ind::Individual{Bool}, penalty::Float64=10.0)
    vmap = collect(1:size(ind.chromosome, 1))
    return fitness_max_ind_set!(g, w, ind, vmap, penalty)
end

function fitness!(mwis::MaxWeightIndependentSet, g::SimpleGraph{Int64}, w::Vector{Float64}, ind::Individual{Bool})
    return fitness_max_ind_set!(g, w, ind, mwis.penalty)
end

function fix_interface!(mwis::MaxWeightIndependentSet, g::SimpleGraph{Int64}, w::Vector{Float64}, ind::Individual{Bool}, interface::Vector{Vector{Int64}})
    w_in_left = 0.0
    for node in outneighbors(g, interface[1][1])
        if node != interface[2][1] && ind.chromosome[node] == true
            w_in_left += 10
        end
    end

    w_in_right = 0.0
    for node in outneighbors(g, interface[2][1])
        if node != interface[1][1] && ind.chromosome[node] == true
            w_in_right += 10
        end
    end

    if ind.chromosome[interface[1][1]] == false && ind.chromosome[interface[2][1]] == false
        if w_in_left == 0.0
            ind.chromosome[interface[1][1]]=true
        elseif w_in_right == 0.0
            ind.chromosome[interface[2][1]]=true
        end
    end

    if ind.chromosome[interface[1][1]] == true && ind.chromosome[interface[2][1]] == true
        #print("Conflict in interface!\n")
        if w_in_left >= w_in_right
            ind.chromosome[interface[1][1]]=false
        else
            ind.chromosome[interface[2][1]]=false
        end
    end
end

"""
Computes the Maximum Weighted Indepndent Set over a tree using dynamic programing
"""
function max_ind_set_on_trees(g::SimpleGraph{Int64}, w::Vector{Float64})
    visited = zeros(Bool, nv(g))
    
    function max_ind_set_on_trees(t::Int64)
        visited[t] = true
        num_children = 0
        neighbors = outneighbors(g, t)

        max_in = w[t]
        max_out = 0.0

        for node in neighbors
            if visited[node] == false
                max_in_ch, max_out_ch = max_ind_set_on_trees(node)
                max_in += max_out_ch 
                max_out += max(max_in_ch, max_out_ch)
            end            
        end
        
        return max_in, max_out
    end

    max_in_f, max_out_f = max_ind_set_on_trees(1)
    return max(max_in_f, max_out_f)
end


######Min Vertex Cover############


function fitness_min_vertex_cover!(g::SimpleGraph{Int64}, w::Vector{Float64}, ind::Individual{Bool}, vmap::Vector{Int64}, penalty::Float64=10.0)
    n = nv(g)
    fitness = 0.0
    penalization = 0.0
    missing_edges = 0
    is_VC= true
    for e in (edges(g))
        s = src(e)
        d = dst(e)
        if (ind.chromosome[vmap[s]] == true) || (ind.chromosome[vmap[d]] == true)
            continue
        else
            missing_edges += 1    
            is_VC = false
            penalization = penalization +penalty+ 2*(w[vmap[s]]+ w[vmap[d]])
            
        end
    end
    if is_VC
        #fitness = sum(ind.chromosome.*w)
        fitness = sum((w.*ind.chromosome)[vmap])
    else
        
        #fitness = penalization + sum(ind.chromosome.*w)
        fitness = penalization + sum((w.*ind.chromosome)[vmap])
    end
    ind.fitness = -fitness
    #ind.computed = true
    
    return -fitness 
  
end

function fitness!(mwvc::MinWeightVertexCover, g::SimpleGraph{Int64}, w::Vector{Float64}, ind::Individual{Bool}, vmap::Vector{Int64}, penalty::Float64=10.0)
    return fitness_min_vertex_cover!(g, w, ind, vmap, penalty)
end

function fitness_min_vertex_cover!(g::SimpleGraph{Int64}, w::Vector{Float64}, ind::Individual{Bool}, penalty::Float64=10.0)
    vmap = collect(1:size(ind.chromosome, 1))
    return fitness_min_vertex_cover!(g, w, ind, vmap, penalty)
end

function fitness!(mwvc::MinWeightVertexCover, g::SimpleGraph{Int64}, w::Vector{Float64}, ind::Individual{Bool}, penalty::Float64=10.0)
    return fitness_min_vertex_cover!(g, w, ind, penalty)
end

function fix_interface_or!(mwvc::MinWeightVertexCover, g::SimpleGraph{Int64}, w::Vector{Float64}, ind::Individual{Bool}, interface::Vector{Vector{Int64}})

    w_in_neighbors = zeros(Float64, size(interface,1))
    n_kids = zeros(Int64, size(interface,1))
    all_kids_are_one = zeros(Bool, size(interface,1))
    for i in 1:size(interface,1)
        neighbors = outneighbors(g, interface[i][1])
        n_kids[i] = size(neighbors,1)
        for node in neighbors
            if ind.chromosome[node] == true && node != interface[3-i][1]
                w_in_neighbors[i] += 1.0
            end
        end
        if w_in_neighbors[i] == n_kids[i]-1 
            all_kids_are_one[i] = true
        end
    end
    #CHECK IF ALL KIDS ARE 1
    #if w_in_neighbors[1] == n_kids[1]-1 
    #    all_kids_are_one[1] = true
    #end
    #if w_in_neighbors[2] == n_kids[2]-1
    #    all_kids_are_one[2] = true
    #end

    
    #Both in interface are zero
    if ind.chromosome[interface[1][1]] == false && ind.chromosome[interface[2][1]] == false
        
        for i in size(interface,1)
            if  w_in_neighbors[i] != 0.0
                if all_kids_are_one[i] == true 
                    #all kids are one
                    #get sure just one in interface is true
                    #case 1.1
                    ind.chromosome[interface[i][1]] = false
                    ind.chromosome[interface[3-i][1]] = true
                else
                    #not all kids are one
                    #case 1.2
                    ind.chromosome[interface[i][1]] = true
                end
            else
                #all kids are 0
                #case 1.3
                ind.chromosome[interface[i][1]]=true
            end

        end

    elseif ind.chromosome[interface[1][1]] == true && ind.chromosome[interface[2][1]] == true
        for  i in size(interface,1)
            if  w_in_neighbors[i] != 0.0 #some kids aren't zero 
                if all_kids_are_one[i] == true 
                    #all kids are one
                    #get sure jus one in interface is true
                    #r = rand(1:2)
                    #case 2.1
                    ind.chromosome[interface[i][1]] = false
                    ind.chromosome[interface[3-i][1]] = true
                else

                    #case 2.2
                    #not all kids are one
                    ind.chromosome[interface[i][1]] = true
                end
            else
                #all kids are 0
                #case 2.3
                ind.chromosome[interface[i][1]]=true
            end
        end

    else #one of the nodes of interface is 0 and the other is 1
        #get the index of the one that is 1
        if ind.chromosome[interface[1][1]]==true
            ind_t=1
        else
            ind_t=2
        end
    
        #case 3.1 do nothing
        if w_in_neighbors[1]==0 && w_in_neighbors[2]==0 #all kids are zero
            #case 3.2
            ind.chromosome[interface[3-ind_t][1]]=true
        end
        if w_in_neighbors[3-ind_t]==0 #kids of the interface node of zero is also zeros
           #case 3.3.a
            ind.chromosome[interface[3-ind_t][1]]=true
            #case 3.3.b
            if all_kids_are_one[ind_t]==true # 
                ind.chromosome[interface[ind_t][1]]=false
            end
        end
    end
end



function fix_interface!(mwvc::MinWeightVertexCover, g::SimpleGraph{Int64}, w::Vector{Float64}, ind::Individual{Bool}, interface::Vector{Vector{Int64}})
    #print(ind.chromosome)
    w_in_neighbors = zeros(Float64, size(interface,1))
    n_kids = zeros(Int64, size(interface,1))
    all_kids_are_one = zeros(Bool, size(interface,1))
    for i in 1:size(interface,1)
        neighbors = outneighbors(g, interface[i][1])
        n_kids[i] = size(neighbors,1)-1
        for node in neighbors
            if ind.chromosome[node] == true && node != interface[3-i][1]
                w_in_neighbors[i] += 1.0
            end
        end
        if w_in_neighbors[i] == n_kids[i] 
            all_kids_are_one[i] = true
        end
    end
    #CHECK IF ALL KIDS ARE 1
    #if w_in_neighbors[1] == n_kids[1]-1 
    #    all_kids_are_one[1] = true
    #end
    #if w_in_neighbors[2] == n_kids[2]-1
    #    all_kids_are_one[2] = true
    #end

    
    #Both in interface are zero
    if ind.chromosome[interface[1][1]] == false && ind.chromosome[interface[2][1]] == false
        #print(ind.chromosome)
        if all_kids_are_one[1] == true && all_kids_are_one[2] == true
            #case1.1
            #all kids are one, choose random
            r = rand(1:2)
            ind.chromosome[interface[r][1]] = true
        elseif w_in_neighbors[1] == 0 && w_in_neighbors[2] == 0
            #all kids are zero
            #case1.3
            ind.chromosome[interface[1][1]] = true
            ind.chromosome[interface[2][1]] = true
        else
            #not all kids are one
            #case 1.2
            if all_kids_are_one[2]==true
                #first not one, second is one
                ind.chromosome[interface[1][1]] = true
            elseif all_kids_are_one[1]==true
                #second not one, first is one
                ind.chromosome[interface[2][1]] = true
            else
                #some of the kids of both may not be one
                ind.chromosome[interface[1][1]] = true
                ind.chromosome[interface[2][1]] = true

            end
        end

    elseif ind.chromosome[interface[1][1]] == true && ind.chromosome[interface[2][1]] == true
        
        if all_kids_are_one[1] == true && all_kids_are_one[2] == true
            #case 2.1
            r = rand(1:2)
            ind.chromosome[interface[r][1]] = false
            #ind.chromosome[interface[3-r][1]] = true
        
        elseif w_in_neighbors[1] == 0 && w_in_neighbors[2] == 0
            #case 2.3
            #do nothing
            ind.chromosome[interface[1][1]] = true
            ind.chromosome[interface[2][1]] = true
        else
            if all_kids_are_one[2] == true
                #case 2.2
                #ind.chromosome[interface[1][1]] = true
                ind.chromosome[interface[2][1]] = false
            elseif all_kids_are_one[1] == true
                #case 2.3
                ind.chromosome[interface[1][1]] = false
                #ind.chromosome[interface[2][1]] = true
            end
        
        end

    else #one of the nodes of interface is 0 and the other is 1
        #get the index of the one that is 1
        if ind.chromosome[interface[1][1]]==true
            ind_t=1
        else
            ind_t=2
        end

        #case 3.1 do nothing
        if w_in_neighbors[1]==0 && w_in_neighbors[2]==0 #all kids are zero
            #case 3.2
            ind.chromosome[interface[3-ind_t][1]]=true
        
        elseif all_kids_are_one[3-ind_t]==false#kids of the interface node of zero is also zeros
           #case 3.3.a
            ind.chromosome[interface[3-ind_t][1]]=true
            #case 3.3.b
            if all_kids_are_one[ind_t]==true # 
                
                ind.chromosome[interface[ind_t][1]]=false
            end
        end
    end

    #print(ind.chromosome)
end




"""
Calculate Min Vertex Cover with DP on trees
"""


function min_vertex_cover_on_trees(g::SimpleGraph{Int64}, w::Vector{Float64},root::Int64)
    visited = zeros(Bool, nv(g))
    
    function min_vertex_cover_on_trees(t::Int64)
        visited[t] = true
        num_children = 0
        neighbors = outneighbors(g, t)

        min_in = w[t]
        min_out = 0.0

        for node in neighbors
            if visited[node] == false
                min_in_ch, min_out_ch = min_vertex_cover_on_trees(node)
                min_in +=  min(min_in_ch, min_out_ch)
                min_out += min_in_ch
            end            
        end
        
        return min_in, min_out
    end

    min_in_f, min_out_f = min_vertex_cover_on_trees(root)
    return min(min_in_f, min_out_f)
end

