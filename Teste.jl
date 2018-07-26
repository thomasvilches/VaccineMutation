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

Vector_Prob = zeros(Float64,365)

for i = 1:365
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

A = Vector{Int64}(100000)
for i = 1:100000
A[i] = Which_One_Will_Transmit(VaccineEfVector,Vector_time,9,2)
end
find(x-> x == 2,A)
r = rand()
retorno = findfirst(x->x>r,probs)



function main(simulationNumber::Int64,P::InfluenzaParameters)
  
    humans = Array{Human}(P.grid_size_human)
    
   # srand(100*simulationNumber)

    setup_human(humans)
    setup_demographic(humans,P)
    Vaccine_Strain = Vector{Int8}(P.sequence_size)
    Creating_Vaccine_Vector(Vaccine_Strain,P)

    if P.GeneralCoverage == 1
        vaccination(humans,P)
    end
    latent_ctr = zeros(Int64,P.sim_time)##vector of results latent
    symp_ctr = zeros(Int64,P.sim_time) #vector for results symp 
    asymp_ctr = zeros(Int64,P.sim_time) #vector for results asymp 
    gd_ctr = zeros(Float64,P.sim_time) #vector for results genetic distance

    initial = setup_rand_initial_latent(humans,P,Vaccine_Strain)### for now, we are using only 1
    Number_in_age_group = zeros(Int64,15)
    Age_group_Matrix = Matrix{Int64}(15,P.grid_size_human)
    for i = 1:P.grid_size_human
        Age_group_Matrix[humans[i].contact_group,(Number_in_age_group[humans[i].contact_group]+1)] = humans[i].index
        Number_in_age_group[humans[i].contact_group] += 1
    end

    for t=1:P.sim_time
        
        contact_dynamic2(humans,P,Age_group_Matrix,Number_in_age_group,Vaccine_Strain)
        
        for i=1:P.grid_size_human
            increase_timestate(humans[i],P)
        end
        latent_ctr[t],symp_ctr[t],asymp_ctr[t],gd_ctr[t]=update_human(humans,P,Vaccine_Strain)
        
    end
    first_inf = find(x-> x.WhoInf == initial && x.WentTo == SYMP,humans)

    numb_first_inf = length(first_inf)

    ## Calculating the proportion of infected people in function of hamming distance

    #number_of_infected = sum(latent_ctr)
    p = zeros(Int64,P.matrix_strain_lines)
    Ef = zeros(Float64,P.grid_size_human)
    count::Int64 = 0
    for i = 1:P.grid_size_human
        if humans[i].WhoInf > 0
            count += 1
            p[Int64(Calculating_Distance_Two_Strains(Vaccine_Strain,humans[i].strains_matrix[1,:]))+1]+=1
            auxMatrix = Matrix{Int8}(1,P.sequence_size)
            auxMatrix[1,:] = humans[i].strains_matrix[1,:]
            Ef[count] = Calculating_Efficacy(auxMatrix,1,Vaccine_Strain,humans[i].vaccineEfficacy,P)[1]
        end
    end
    
    #return latent_ctr,symp_ctr,asymp_ctr,numb_first_inf,numb_symp_inf,numb_asymp_inf
    return latent_ctr,symp_ctr,asymp_ctr,numb_first_inf,p,gd_ctr,Ef

end

