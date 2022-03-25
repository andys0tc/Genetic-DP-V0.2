using Plots
using Random, Distributions, StatsBase
#using PyPlot
import Distributions: Uniform 

function exponential_dist(values::Vector{Float64})
    val = zeros(size(values,1))
    #local val =values
    for i in 1:size(val, 1)
        val[i] = exp(values[i])+floatmin()
    end
    s = sum(val)
    for i in 1:size(val, 1)
        val[i] = val[i]/s
    end
    return val
end

function exponential(values::Vector{Float64})
    val = zeros(size(values,1))
    for i in 1:size(values, 1)
        val[i] = exp(values[i])+floatmin()
    end
    return val
end

function scaleVector(values::Vector{Float64})
    
    min = minimum(values)
    max = maximum(values)
    val = (values.-min)./(max-min)
    return val
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





x_100 = Array{Float64}(collect(1:1:100));
#x_1000 = Array{Float64}(collect(1:1:1000));

vec_100 = rand(Float64(1):Float64(100),(100));
#vec_1000 = rand(Float64(1):Float64(100),(1000));

vec_interval_100 = rand(Uniform(-1,1), 100);
#vec_interval_1000 = rand(Uniform(-1,1), 1000); 

vec_rand_100_pos=rand(100);
vec_rand_100_neg = -1*rand(100);

#vec_rand_1000_pos=rand(1000);
#vec_rand_1000_neg = -1*rand(1000);


vec_100_norm=scaleVector(vec_100);
#vec_1000_norm=scaleVector!(vec_1000);
vec_interval_100_norm=scaleVector(vec_interval_100);
#vec_interval_1000_norm=scaleVector!(vec_interval_1000);


vec_100_exp = exponential_dist(vec_100);
#vec_1000_exp = exponential_dist(vec_1000);
vec_norm_100 = exponential_dist(scaleVector(vec_100));
#vec_norm_1000 = exponential_dist(scaleVector!(vec_1000));

vec_interval_100_exp = exponential_dist(vec_interval_100);
#vec_interval_1000_exp = exponential_dist(vec_interval_1000);
vec_interval_norm_100 = exponential_dist(scaleVector(vec_interval_100));
#vec_interval_norm_1000 = exponential_dist(scaleVector!(vec_interval_1000));



#####DIST count


count_100 = getDist(vec_100_exp,100);
print(count_100)
histogram(randn(1000), bins = :scott, weights = repeat(1:5, outer = 200))


#plt.hist(vec_100_exp)

#plt.show()
#plt.hist(vec_1000_exp)
#plt.show()
#plt.hist(vec_norm_100)
#plt.show()
#plt.hist(vec_norm_1000)
#plt.show()
#plt.hist(vec_interval_100_exp)
#plt.show()
#plt.hist(vec_interval_1000_exp)
#plt.show()
#plt.hist(vec_interval_norm_100)
#plt.show()
#plt.hist(vec_interval_norm_1000)
#plt.show()
#"""

