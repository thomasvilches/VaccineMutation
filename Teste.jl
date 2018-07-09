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
    grid_size_human = 7100

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
Vaccine_Strain = Original_Strain
TransmitingStrain = Original_Strain
t = 7

for i = 1:P.grid_size_human
    humans[i].strains_matrix,humans[i].Vector_time,humans[i].NumberStrains = mutation(TransmitingStrain,P,t,CumProb,1) ### t must be changed by humans[i].infectious period
    available_strains = 1:humans[i].NumberStrains
    VaccineEfVector = Calculating_Efficacy(humans[i].strains_matrix[available_strains,:],length(available_strains),Vaccine_Strain,0.0,P)
    transm = Which_One_Will_Transmit(VaccineEfVector,humans[i].Vector_time[available_strains],8,2)#rand(1:humans[i].NumberStrains)
    TransmitingStrain = humans[i].strains_matrix[transm,:]
end

DistMatrix = Calculating_Distance(humans,P)

DistMatrix/P.sequence_size
B = DistMatrix[1,:]/P.sequence_size

find(x-> x.NumberStrains>2,humans)



########################################################3

Matrix = zeros(Int64,P.matrix_strain_lines,P.sequence_size)
Matrix[1,:] = rand(1:P.number_of_states,P.sequence_size)

Matrix[2,:] = Matrix[1,:]

for i = 1:30
    aux = rand(1:P.sequence_size)
    Matrix[2,aux] = rand(1:P.number_of_states)
end


Matrix[3,:] = Matrix[2,:]

for i = 1:15
    aux = rand(1:P.sequence_size)
    Matrix[3,aux] = rand(1:P.number_of_states)
end

Calculating_Distance_Two_Strains(Matrix[1,:],Matrix[3,:])

VaccineEfVector = Calculating_Efficacy(Matrix,3,Matrix[1,:],0.8,P)
Vector_time = [0; 5; 8]
Which_One_Will_Transmit(VaccineEfVector,Vector_time,11,2)