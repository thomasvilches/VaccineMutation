using ProgressMeter
using PmapProgressMeter
using DataArrays, DataFrames
using Match
using ParallelDataTransfer
using QuadGK
using Parameters #module
using Distributions
using StatsBase

include("parameters.jl")
include("PopStruc.jl")
include("functions.jl")
include("mutation.jl")

function main(cb,simulationNumber::Int64,P::InfluenzaParameters)
  
    humans = Array{Human}(P.grid_size_human)
    
   # srand(100*simulationNumber)

    setup_human(humans)
    setup_demographic(humans,P)
    Vaccine_Strain = Vector{Int8}(P.sequence_size)
    Mod_Strain = Vector{Int8}(P.sequence_size)
    Mod_Strain2 = Vector{Int8}(P.sequence_size)
    Creating_Vaccine_Vector(Vaccine_Strain,P)

    if P.GeneralCoverage == 1
        vaccination(humans,P)
    end
    latent_ctr = zeros(Int64,P.sim_time)##vector of results latent
    symp_ctr = zeros(Int64,P.sim_time) #vector for results symp 
    asymp_ctr = zeros(Int64,P.sim_time) #vector for results asymp 
    gd_ctr = zeros(Float64,P.sim_time) #vector for results genetic distance

    initial = setup_rand_initial_latent(humans,P,Vaccine_Strain)### for now, we are using only 1
    #############################################################33
    ###############################################################
   #= if P.start_different == 1
        Mod_Strain = Vaccine_Strain
        for i = 1:1:(round(P.initial_p0*P.sequence_size))
            aux = rand(1:P.sequence_size)
            change::Int8 =  Mod_Strain[aux]
            while change == Mod_Strain[aux]
                change = rand(1:P.number_of_states)
            end
            Mod_Strain[aux] = change
        end
    end
    Mod = setup_rand_initial_latent(humans,P,Mod_Strain)

   if P.start_different == 1
        Mod_Strain2 = Vaccine_Strain
        for i = 1:1:(round(P.initial_p0*P.sequence_size))
            aux = rand(1:P.sequence_size)
            change::Int8 =  Mod_Strain2[aux]
            while change == Mod_Strain2[aux]
                change = rand(1:P.number_of_states)
            end
            Mod_Strain2[aux] = change
        end
    end
    Mod = setup_rand_initial_latent(humans,P,Mod_Strain2)



   if P.start_different == 1
    Mod_Strain2 = Vaccine_Strain
    for i = 1:1:(round(P.initial_p0*P.sequence_size))
        aux = rand(1:P.sequence_size)
        change::Int8 =  Mod_Strain2[aux]
        while change == Mod_Strain2[aux]
            change = rand(1:P.number_of_states)
        end
        Mod_Strain2[aux] = change
    end
end
Mod = setup_rand_initial_latent(humans,P,Mod_Strain2)
=#
####################################################################################################3
##########################################################################################################

   
   
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
        latent_ctr[t],symp_ctr[t],asymp_ctr[t],gd_ctr[t]= update_human(humans,P,Vaccine_Strain)
        cb(1) ## increase the progress metre by 1.. callback function
    end
    first_inf = find(x-> x.WhoInf == initial && x.WentTo == SYMP,humans)

    numb_first_inf = length(first_inf)

    ## Calculating the proportion of infected people in function of hamming distance

    #number_of_infected = sum(latent_ctr)
    pv = zeros(Int64,P.matrix_strain_lines)
    pnv = zeros(Int64,P.matrix_strain_lines)
    Ef = zeros(Float64,P.matrix_strain_lines)
    count::Int64 = 0
    for i = 1:P.grid_size_human
        if humans[i].WhoInf > 0
            count += 1
            auxP = Int64(Calculating_Distance_Two_Strains(Vaccine_Strain,humans[i].strains_matrix[1,:]))+1
            if humans[i].vaccinationStatus == 1
                pv[auxP]+=1
            else
                pnv[auxP]+=1
            end
            auxMatrix = Matrix{Int8}(1,P.sequence_size)
            auxMatrix[1,:] = humans[i].strains_matrix[1,:]
            Ef[auxP] += Calculating_Efficacy(auxMatrix,1,Vaccine_Strain,humans[i].vaccineEfficacy,P)[1]
        end
    end
    
    #return latent_ctr,symp_ctr,asymp_ctr,numb_first_inf,numb_symp_inf,numb_asymp_inf
    return latent_ctr,symp_ctr,asymp_ctr,numb_first_inf,pv,pnv,gd_ctr,Ef

end

