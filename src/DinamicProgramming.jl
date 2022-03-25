using LightGraphs 
using DelimitedFiles
using Plots


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
                
                #dp_val_ch[node,1]  = min_in_ch
                #dp_val_ch[node,2]  = min_out_ch
                
            end   
            
                   
        end
        #dp_val[t,1] = min_in
        #dp_val[t,2] = min_out
        return min_in, min_out
    end

    min_in_f, min_out_f = min_vertex_cover_on_trees(root)
    #dp_val[root,1] = min_in_f
    #dp_val[root,2] = min_out_f
    return (min_in_f, min_out_f)
end

function fitness_min_vertex_cover!(g::SimpleGraph{Int64}, w::Vector{Float64}, ind::Vector{Float64}, vmap::Vector{Int64}, penalty::Float64=10.0)
    n = nv(g)
    fitness = 0.0
    penalization = 0.0
    missing_edges = 0
    is_VC= true
    for e in (edges(g))
        s = src(e)
        d = dst(e)
        if (ind[vmap[s]] == true) || (ind[vmap[d]] == true)
            continue
        else
            missing_edges += 1    
            is_VC = false
            penalization = penalization +penalty+ 2*(w[vmap[s]]+ w[vmap[d]])
            
        end
    end
    if is_VC == true
        #fitness = sum(ind.chromosome.*w)
        fitness = sum((w.*ind)[vmap])      
        
    else
                
        #fitness = penalization + sum(ind.chromosome.*w)
        fitness = penalization + sum((w.*ind)[vmap])
    end
    #ind.fitness = -fitness
    #ind.computed = true
    #if fitness == 0.0
    #    fitness = penalty
    #end
    
    return fitness 
  
end

function opt_tree(g::SimpleGraph{Int64})
    opt = zeros(ne(g),2)
    for (i,e) in enumerate(edges(g))
       
        u, v = src(e), dst(e)

        rem_edge!(g, u, v)
        components = connected_components(g)
        #comps[1] = components[1]
        #comps[2] = components[2]
        interface[1] = Vector{Int64}(undef, 1)
        interface[2] = Vector{Int64}(undef, 1)
        t_1, vmap_t1 = induced_subgraph(g, components[1])
        t_2, vmap_t2 = induced_subgraph(g, components[2])
        
        if u ∈ components[1]
            interface[1][1] = u
            interface[2][1] = v
            comps[1] = components[1]
            comps[2] = components[2]
            g_1 = t_1
            g_2 = t_2
        else
            interface[2][1] = v
            interface[1][1] = u
            comps[1] = components[2]
            comps[2] = components[1]
            g_1=t_2
            g_2=t_1
        end

        
        root_u = findall(x->x==interface[1][1], comps[1])
        #if isempty(root_u)
        #    root_u = findall(x->x==u, vmap_t2)
        #end

        root_v = findall(x->x==interface[2][1], comps[2]) 
        
        #if isempty(root_v)
        #    root_v = findall(x->x==v, vmap_t1)
        #end

        
        
        
        add_edge!(g, u, v)
        #i+=1
    end
    return opt
end 





#for i in 1:size(dp_val, 1)
#
#    print("min edge $i ",minimum(@view dp_val[i,:]))
#    print("\n")
#    print("\n")
#end

#data = readdlm("C:/Users/atorr/OneDrive - Instituto Politecnico Nacional/Documents/Julia Projects/Genetic Dinamic Programming/test/GenDP.txt", ',', Int)
#print(data)
#number of popu
function open_data(file_name::String,n_popu::Int64,n_chrom::Int64)
    
    data = readdlm(file_name, ',', Int)
    history = zeros(n_popu,n_chrom,Int64(size(data, 1)/n_popu))
    println("Reading data")
    j=0
    for i in 0:(size(data, 1)+1)
        if(i==size(data, 1))
            break
        end
        x= (i%n_popu)+1
        z = div(j,n_popu)+1 #generacion

        #if j%n_popu==0
        #    print("\n")
        #    print("Generation ", z)
        #    print("\n")
        #end
        history[x,:,z] = data[i+1,:]
        #println(data[i+1,:,])
        j+=1   
    end
    #print(history)
    return history
end

function check_solutions(g::SimpleGraph{Int64},data::Array{Float64})
    #missing_edges = 0
    #popu,2particiones,aristas,generaciones
    #30*2*5*9
    opt = zeros(ne(g),2,2)
    w = ones(nv(g))
    fitness=zeros(size(data,1),2,size(data,3),ne(g),2) # 1 or 0 fixed for last dimension
    #i=1
    for (i,e) in enumerate(edges(g))
       
        u, v = src(e), dst(e)

        rem_edge!(g, u, v)
        interface[1] = Vector{Int64}(undef, 1)
        interface[2] = Vector{Int64}(undef, 1)
        components = connected_components(g)
        comps[1] = components[1]
        comps[2] = components[2]
        t_1, vmap_t1 = induced_subgraph(g, comps[1])
        t_2, vmap_t2 = induced_subgraph(g, comps[2])
        
        if u ∈ components[1]
            interface[1][1] = u
            interface[2][1] = v
            comps[1] = components[1]
            comps[2] = components[2]
            g_1 = t_1
            g_2 = t_2
            uv_edges[i,1] = u
            uv_edges[i,2] = v
        else
            interface[2][1] = v
            interface[1][1] = u
            comps[1] = components[2]
            comps[2] = components[1]
            g_1=t_2
            g_2=t_1
            uv_edges[i,2] = v
            uv_edges[i,1] = u
        end
        root_u = findall(x->x==interface[1][1], comps[1])
        root_v = findall(x->x==interface[2][1], comps[2]) 

        opt[i,1,1],opt[i,1,2] = min_vertex_cover_on_trees(g_1, ones(size(g_1, 1)), root_u[1])
        opt[i,2,1],opt[i,2,2] = min_vertex_cover_on_trees(g_2, ones(size(g_2, 1)), root_v[1])

        #opt[i,1] = min_vertex_cover_on_trees(g_1, ones(size(g_1, 1)), comps[1][root_u[1]])
        #opt[i,2] = min_vertex_cover_on_trees(g_2, ones(size(g_2, 1)), comps[2][root_v[1]])

        for j in 1:size(fitness,3)#generation            
            for l in 1:size(fitness,1)#popu
                old_1 = data[:,interface[1],j]
                #fixing with zero
                data[:,interface[1],j] .= 0   
                fitness[l,1,j,i,1]=fitness_min_vertex_cover!(g_1, w,data[l,:,j], comps[1])
                #fixing with one
                data[:,interface[1],j] .= 1   
                fitness[l,1,j,i,2]=fitness_min_vertex_cover!(g_1, w,data[l,:,j], comps[1])

                data[:,interface[1],j] = old_1
                #if (fitness[l,1,j,i]>=5)
                #    println(comps[1])
                #    println(data[l,:,j])
                #    println(fitness[l,1,j,i])
                #end
                old_2 = data[:,interface[2],j]
                data[:,interface[2],j] .= 0 
                fitness[l,2,j,i,1]=fitness_min_vertex_cover!(g_2, w,data[l,:,j], comps[2])

                data[:,interface[2],j] .= 1
                fitness[l,2,j,i,2]=fitness_min_vertex_cover!(g_2, w,data[l,:,j], comps[2])
                data[:,interface[2],j] = old_2

                #if (fitness[l,2,j,i]>=6)
                #    println(comps[2])
                #    println(data[l,:,j])
                #    println(fitness[l,2,j,i])
                #end
            end
        end
        add_edge!(g, u, v)
        #i+=1
    end
    return fitness,opt
end

function get_best_per_generation_working(g::SimpleGraph{Int64},fitness::Array{Float64},data::Array{Float64})
    best_per_gen = zeros(2,2,size(data,3),ne(g))
    n_edges=1
    for e in edges(g)#edge
        u, v = src(e), dst(e)

        for i in 1:size(data,3)#generation
            
            #se puede hacer mas eficiente
            index_ones_u =  findall(x->x==1.0, data[:,u,i])
            index_zeros_u =  findall(x->x==0.0, data[:,u,i])
            index_ones_v =  findall(x->x==1.0, data[:,v,i])
            index_zeros_v =  findall(x->x==0.0, data[:,v,i])
            best_in_cover_u = fitness[:,1,i,n_edges][index_ones_u]
            best_out_cover_u = fitness[:,1,i,n_edges][index_zeros_u]
            best_in_cover_v = fitness[:,2,i,n_edges][index_ones_v]
            best_out_cover_v = fitness[:,2,i,n_edges][index_zeros_v]
            try
                best_per_gen[1,1,i,n_edges] = minimum(best_in_cover_u)
            catch e
                best_per_gen[1,1,i,n_edges] = 0
            end
            try
                best_per_gen[1,2,i,n_edges] = minimum(best_out_cover_u)
            catch e
                best_per_gen[1,2,i,n_edges] = 0
            end
            try
                best_per_gen[2,1,i,n_edges] = minimum(best_in_cover_v)
            catch e
                best_per_gen[2,1,i,n_edges] = 0
            end
            try
                best_per_gen[2,2,i,n_edges] = minimum(best_out_cover_v)
            catch e
                best_per_gen[2,2,i,n_edges] = 0
                      end
            #best_per_gen[1,2,i,n_edges] = minimum(best_out_cover_u)
            #best_per_gen[2,1,i,n_edges] = minimum(best_in_cover_v)
            #best_per_gen[2,2,i,n_edges] = minimum(best_out_cover_v)
            #best_in_index = argmin(best_in_cover)
            #for k in 1:size(fitness,2)#partition
            #    best_per_gen[k,j,i,:]=data[j,:,i]
            #end
            
        end
        n_edges+=1
    end
    return best_per_gen
end


function get_best_per_generation(g::SimpleGraph{Int64},fitness::Array{Float64},data::Array{Float64})
    #best_per_gen = zeros(2,2,size(data,3),ne(g))
    best_per_gen = Array{Any}(undef,2,2,size(data,3),ne(g))
    #n_edges=1
    for (n_edges,e) in enumerate(eachrow(uv_edges))#edge
        u, v = e[1], e[2]

        for i in 1:size(data,3)#generation
            
            #se puede hacer mas eficiente
            
            #index_ones_u =  findall(x->x==1.0, data[:,u,i])
            #index_zeros_u =  findall(x->x==0.0, data[:,u,i])
            #index_ones_v =  findall(x->x==1.0, data[:,v,i])
            #index_zeros_v =  findall(x->x==0.0, data[:,v,i])
            #best_in_cover_u = fitness[:,1,i,n_edges][index_ones_u]
            #best_out_cover_u = fitness[:,1,i,n_edges][index_zeros_u]
            #best_in_cover_v = fitness[:,2,i,n_edges][index_ones_v]
            #best_out_cover_v = fitness[:,2,i,n_edges][index_zeros_v]

            best_in_cover_u = fitness[:,1,i,n_edges,2]
            best_out_cover_u = fitness[:,1,i,n_edges,1]
            best_in_cover_v = fitness[:,2,i,n_edges,2]
            best_out_cover_v = fitness[:,2,i,n_edges,1]

            try
                best_per_gen[1,1,i,n_edges] = minimum(best_in_cover_u)
            catch e
                #best_per_gen[1,1,i,n_edges] = missing
                best_per_gen[1,1,i,n_edges] = size(uv_edges,1)*10
            end
            try
                best_per_gen[1,2,i,n_edges] = minimum(best_out_cover_u)
            catch e
                #best_per_gen[1,2,i,n_edges] = missing
                best_per_gen[1,2,i,n_edges] = size(uv_edges,1)*10
            end
            try
                best_per_gen[2,1,i,n_edges] = minimum(best_in_cover_v)
            catch e
                #best_per_gen[2,1,i,n_edges] = missing
                best_per_gen[2,1,i,n_edges] = size(uv_edges,1)*10
            end
            try
                best_per_gen[2,2,i,n_edges] = minimum(best_out_cover_v)
            catch e
                #best_per_gen[2,2,i,n_edges] = missing
                best_per_gen[2,2,i,n_edges] = size(uv_edges,1)*10
            end
            #best_per_gen[1,2,i,n_edges] = minimum(best_out_cover_u)
            #best_per_gen[2,1,i,n_edges] = minimum(best_in_cover_v)
            #best_per_gen[2,2,i,n_edges] = minimum(best_out_cover_v)
            #best_in_index = argmin(best_in_cover)
            #for k in 1:size(fitness,2)#partition
            #    best_per_gen[k,j,i,:]=data[j,:,i]
            #end
            
        end
        #n_edges+=1
    end
    return best_per_gen
end

#######################START##############3

g = loadgraph("C:/Users/atorr/OneDrive - Instituto Politecnico Nacional/Documents/Julia Projects/Genetic Dinamic Programming/test/graphs/gtest80.lgz")

#g = SimpleGraph(10)
#
#add_edge!(g,1,2);
#add_edge!(g,1,3);
#add_edge!(g,2,4);
#add_edge!(g,2,5);
#add_edge!(g,5,7);
#add_edge!(g,5,8);
#add_edge!(g,3,6);
#add_edge!(g,6,9);
#add_edge!(g,6,10);



print("Dinamic Programming")
print("\n")

w = ones(nv(g))
dp_val= zeros(nv(g),2)
opt_val = zeros(nv(g),2)
#dp_val_ch = zeros(nv(g),2)
comps = Vector{Vector{Int64}}(undef, 2)
interface = Vector{Vector{Int64}}(undef, 2)
uv_edges = zeros(Int64,ne(g),2)

#min_vertex_cover_on_trees(g, ones(size(g, 1), ),4)
#opt_val=opt_tree(g)


print(dp_val)
print("\n")
print("\n")
#print(dp_val_ch)
print("\n")


history  = open_data("C:/Users/atorr/OneDrive - Instituto Politecnico Nacional/Documents/Julia Projects/Genetic Dinamic Programming V0.2/test/history/GenDP80.txt",nv(g)*3,nv(g))
fitness_solutions,opt_val = check_solutions(g,history)

best_gen = get_best_per_generation(g,fitness_solutions,history)

x = collect(1:1:size(history,3))
#
#for e in edges(g)
#for i in 1:ne(g)
for (i,e) in enumerate(edges(g))
    u, v = src(e), dst(e)
    y1 = zeros(1,size(best_gen,3))
    y2 = zeros(1,size(best_gen,3))
    y3 = zeros(1,size(best_gen,3))
    y4 = zeros(1,size(best_gen,3))

    try
        #y1=opt_val[i,1] ./ best_gen[1,1,:,i]
        y1=best_gen[1,1,:,i]
        replace!(y1, NaN=>1)
        replace!(y1, Inf=>0)
        #y1=best_gen[1,1,:,i]./opt_val[i,1]#in u
    catch
       
    end
    try
        #y2=opt_val[i,1] ./ best_gen[1,2,:,i]
        y2 = best_gen[1,2,:,i]
        replace!(y2, NaN=>1)
        replace!(y2, Inf=>0)
        #y2=best_gen[1,2,:,i]./opt_val[i,1] #out u
    catch
        
        replace!(y2, NaN=>1)
    end
    try
        #y3=opt_val[i,2] ./ best_gen[2,1,:,i]
        y3 = best_gen[2,1,:,i]
        replace!(y3, NaN=>1)
        replace!(y3, Inf=>0)
        #y3=best_gen[2,1,:,i]./opt_val[i,2] #in v
    catch
        
        
    end
    try
        #y4=opt_val[i,2] ./ best_gen[2,2,:,i]
        y4=best_gen[2,2,:,i]
        replace!(y4, NaN=>1)
        replace!(y4, Inf=>0)
        #y4=best_gen[2,2,:,i]./opt_val[i,2] #out v
    catch
        
        
    end
    opt_u_in= opt_val[i,1,1]
    opt_u_out = opt_val[i,1,2]
    opt_v_in =opt_val[i,2,1]
    opt_v_out = opt_val[i,2,2]

    #y2=best_gen[1,2,:,i]#out u
    #y3=best_gen[2,1,:,i]#in v
    #y4=best_gen[2,2,:,i]#out v
    (plot(x,[y1,y2], title = "Results GenDP over Edge $u-$v, Node $u ",label=["In $opt_u_in" "Out $opt_u_out" ]))
    savefig("C:/Users/atorr/OneDrive - Instituto Politecnico Nacional/Documents/Julia Projects/Genetic Dinamic Programming V0.2/test/figures/Tree80GenDPUV2/myplot$u-$v 1.png")
    #y4=best_gen[2,2,:,i]#out v
    (plot(x,[y3,y4], title = "Results GenDP over Edges $v-$u,  Node $v ",label=["In $opt_v_in" "Out $opt_v_out" ]))
    savefig("C:/Users/atorr/OneDrive - Instituto Politecnico Nacional/Documents/Julia Projects/Genetic Dinamic Programming V0.2/test/figures/Tree80GenDPUV2/myplot$u-$v 2.png")
    #i=+1
end


#println(size(fitness_solutions))



#open("C:/Users/atorr/OneDrive - Instituto Politecnico Nacional/Documents/Julia Projects/Genetic Dinamic Programming/test/GenDP.txt") do file
#    gen=30
#    j=0
#    for l in eachline(file)
#        if j%gen==0
#            print("\n")
#            print("Generation ", div(j,gen)+1)
#            print("\n")
#        end
#        
#        println(l)
#        j+=1
#    end
#    # do stuff with the open file
#end

#for e in edges(g);
#    u, v = src(e), dst(e)
#    #println("edge $u - $v")
#    rem_edge!(g, u, v)
#
#    components = connected_components(g)
#    comps[1] = components[1]
#    comps[2] = components[2]
#    t_1, vmap_t1 = induced_subgraph(g, comps[1])
#    t_2, vmap_t2 = induced_subgraph(g, comps[2])
#    #Check if t is not a tree
#
#    root_u = findall(x->x==u, vmap_t1)
#    opt_1=min_vertex_cover_on_trees(t_1, ones(size(comps[1], 1), ),root_u[1])
#    root_v = findall(x->x==v, vmap_t2)
#    opt_2=min_vertex_cover_on_trees(t_2, ones(size(comps[2], 1), ),root_v[1])
#
#    #if size(comps[1], 1) >= 2
#    #    root_u = findall(x->x==u, vmap_t1)
##
#    #    opt_1=min_vertex_cover_on_trees(t_1, ones(size(comps[1], 1), ),root_u[1])
#    #else
#    #    opt_1=0
#    #end
#    #
#    #if size(comps[2], 1)>=2
#    #    root_v = findall(x->x==v, vmap_t2)
#    #    opt_2=min_vertex_cover_on_trees(t_2, ones(size(comps[2], 1), ),root_v[1])
#    #else
#    #    opt_2=0
#    #end
#    
#    #opt_1=min_vertex_cover_on_trees(t_1, ones(size(comps[1], 1), ),u)
#    #opt_2=min_vertex_cover_on_trees(t_2, ones(size(comps[2], 1), ),v)
#    print("opt_1 $u = $opt_1, opt_2 $v = $opt_2 \n")
#    #dp_val[e,1] = opt_1
#    #dp_val[e,2] = opt_2
#    add_edge!(g, u, v)
#end
#
#include("FitnessFunctions.jl")