using DataFrames

function random_tree_experiments(tree_sizes::Vector{Int64}, avg_m::Vector{Float64}, fₙ::Vector{Float64}, gen::Vector{Float64}, trials=6, num_per_size=5)
    
    max_dif = typemin(Float64)
    max_dif_g = SimpleGraph{Int64}()
    min_dif = typemax(Float64)
    min_dif_g = SimpleGraph{Int64}()
    lk = ReentrantLock()
    df = DataFrame(tree_size = Int64[], instance = Int64[], MIS = Float64[], p_mutation = Float64[], population = Float64[], generations = Float64[], trial = Int64[], NaiveGA = Float64[], GenDPHalf = Float64[], GenDP = Float64[])

    for t_size in tree_sizes
        Threads.@threads for j in 1:num_per_size
            for am in avg_m
                for f in fₙ
                    for ng in gen 
                        g = random_tree(t_size)
                        w = ones(t_size)
                        opt = max_ind_set_on_trees(g, w)
                        N = floor(f * t_size)
                        pm = am / t_size
                        gens = Int(floor(t_size * log(10, t_size) * ng))
                        for i in 1:trials
                            g_copy = copy(g)       
                            mwis = MaxWeightIndependentSet()
                            popu = GeneticGraphData{Bool}(Int(N), [true, false], g_copy, w, mwis)         
                            
                            nga = NaiveGA{Bool}(pm, popu)
                            naive_ga = run(nga, gens)

                            ggdp = GenDP_Global{Bool}(pm, popu, random_cut_tree)
                            gendp_global = run(ggdp, gens)
                            
                            lgdp = GenDP{Bool}(pm, popu, random_cut_tree)
                            gendp = run(lgdp, gens)                      

                            begin
                                lock(lk)
                                try
                                    push!(df, [t_size j opt pm N gens i naive_ga gendp_global gendp])
                                    print("$([t_size j opt pm N gens i naive_ga gendp_global gendp])\n")
                                finally
                                    unlock(lk)
                                end
                            end
                            
                        end
                    end
                end
            end
        end

        # savegraph("test/graphs/g$(t_size)maxdif.lgz", max_dif_g)
        # savegraph("test/graphs/g$(t_size)mindif.lgz", min_dif_g)
    end
    df.GAAppRatio = df.NaiveGA ./ df.MIS
    df.GenDPHalfAppRatio = df.GenDPHalf ./ df.MIS
    df.GenDPAppRatio = df.GenDP ./ df.MIS
    return df

end