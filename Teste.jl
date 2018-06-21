using ProgressMeter
using PmapProgressMeter
using Parameters
using DataArrays,DataFrames
using QuadGK
using Distributions
using StatsBase
using ParallelDataTransfer
using Match
using Lumberjack
using FileIO

include("parameters.jl")
include("PopStruc.jl")
include("mutation.jl")
include("functions.jl")


P = InfluenzaParameters(
    mutation_rate = 0.00416,
    matrix_strain_lines = 10,
    grid_size_human = 1760

)

#Here I am creating a random strain with sequence_size sites and 1 to 20 states per site

Original_Strain = rand(1:P.number_of_states,P.sequence_size)
p = 1 - exp(-P.mutation_rate/365)

Vector_Prob = zeros(Float64,P.max_infectious_period)

for i = 1:P.max_infectious_period
    Vector_Prob[i] = p*((1-p)^(i-1))
end

CumProb = cumsum(Vector_Prob)

humans = Array{Human}(P.grid_size_human)
setup_human(humans)

###### Let's test the functions and how the mutations is going

TransmitingStrain = Original_Strain
t = 7

for i = 1:P.grid_size_human
    humans[i].strains_matrix,humans[i].Vector_time,humans[i].NumberStrains = mutation(TransmitingStrain,P,t,CumProb,1) ### t must be changed by h[i].infectious period
    transm = rand(1:humans[i].NumberStrains)
    TransmitingStrain = humans[i].strains_matrix[transm,:]
end

DistMatrix = Calculating_Distance(humans,P)

DistMatrix/P.sequence_size
B = DistMatrix[1,:]/P.sequence_size

find(x-> x.NumberStrains>2,humans)