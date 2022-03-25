using LightGraphs, Random, DataStructures

function rand_binary_tree(h::Int64)
    n = 2^h-1
    labels = collect(1:n)
    shuffle!(labels)
    g = SimpleGraph(n)
    k = 0

    function rec_binary_tree(lvl::Int64)
        if lvl == 1
            k += 1
            return k 
        end

        left = rec_binary_tree(lvl-1)
        right = rec_binary_tree(lvl-1)
        k += 1
        add_edge!(g, labels[k], labels[left])
        add_edge!(g, labels[k], labels[right])
        return k
    end
    
    rec_binary_tree(h)
    return g
end

""" Generates a random tree on [n] nodes uniformly distributed among
all possible n^(n-2) possible trees.
"""
function random_tree(n::Int64)
    g = SimpleGraph(n)
    uf = IntDisjointSets(n)

    i = 0

    while i < n - 1
        a = rand(1:n)
        b = rand(1:n)
        if a == b
            continue
        end
        if find_root!(uf, a) != find_root!(uf, b)
            add_edge!(g, a, b)
            i += 1
            union!(uf, a, b)
        end

    end

    return g

end

