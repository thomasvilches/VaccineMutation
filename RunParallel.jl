addprocs(4)

@everywhere using ProgressMeter
@everywhere using PmapProgressMeter
@everywhere using Parameters
@everywhere using DataArrays,DataFrames
@everywhere using QuadGK
@everywhere using Distributions
@everywhere using StatsBase
@everywhere using ParallelDataTransfer
@everywhere using Match
@everywhere using Lumberjack
@everywhere using FileIO
println("added $(nworkers()) processors")
info("starting @everywhere include process...")
@everywhere include("basicModelSlurm.jl")
################## To run this files, You must check the return of BasicModel.jl
#######################3

function dataprocess(results,P::InfluenzaParameters,numberofsims)

    resultsL = Matrix{Int64}(P.sim_time,numberofsims)
    resultsA = Matrix{Int64}(P.sim_time,numberofsims)
    resultsS = Matrix{Int64}(P.sim_time,numberofsims)
    resultsGD = Matrix{Float64}(P.sim_time,numberofsims)
    resultsR0 = Vector{Int64}(numberofsims)
    resultsSymp = Vector{Int64}(numberofsims)
    resultsAsymp = Vector{Int64}(numberofsims)
    
    resultsP = Matrix{Float64}(P.grid_size_human,numberofsims)

    for i=1:numberofsims
        resultsL[:,i] = results[i][1]
        resultsS[:,i] = results[i][2]
        resultsA[:,i] = results[i][3]
        resultsR0[i] = results[i][4]
        resultsSymp[i] = results[i][5]
        resultsAsymp[i] = results[i][6]
        resultsP[:,i] = results[i][7]
        resultsGD[:,i] = results[i][8]

    end
   
    directory = "July9/"

    writedlm(string("$directory","result","$(P.Prob_transmission)","Mut","$(P.mutation_rate)","Ef","$(P.VaccineEfficacy)","_latent.dat"),resultsL)
    writedlm(string("$directory","result","$(P.Prob_transmission)","Mut","$(P.mutation_rate)","Ef","$(P.VaccineEfficacy)","_symp.dat"),resultsS)
    writedlm(string("$directory","result","$(P.Prob_transmission)","Mut","$(P.mutation_rate)","Ef","$(P.VaccineEfficacy)","_asymp.dat"),resultsA)
    writedlm(string("$directory","result","$(P.Prob_transmission)","Mut","$(P.mutation_rate)","Ef","$(P.VaccineEfficacy)","_R0.dat"),resultsR0)
    writedlm(string("$directory","result","$(P.Prob_transmission)","Mut","$(P.mutation_rate)","Ef","$(P.VaccineEfficacy)","_SympInf.dat"),resultsSymp)
    writedlm(string("$directory","result","$(P.Prob_transmission)","Mut","$(P.mutation_rate)","Ef","$(P.VaccineEfficacy)","_AsympInf.dat"),resultsAsymp)
    writedlm(string("$directory","result","$(P.Prob_transmission)","Mut","$(P.mutation_rate)","Ef","$(P.VaccineEfficacy)","_P.dat"),resultsP)
    writedlm(string("$directory","result","$(P.Prob_transmission)","Mut","$(P.mutation_rate)","Ef","$(P.VaccineEfficacy)","_GD.dat"),resultsGD)
end

function run_main(P::InfluenzaParameters,numberofsims::Int64)
    
    results = pmap((cb, x) -> main(cb, x, P), Progress(numberofsims*P.sim_time), 1:numberofsims, passcallback=true)

    dataprocess(results,P,numberofsims)
end


@everywhere P=InfluenzaParameters(
    VaccineEfficacy = 0.2,
    GeneralCoverage = 1,
    Prob_transmission = 0.079,
    sim_time = 200,
    grid_size_human = 1000,
    mutation_rate = 0.00416
)

run_main(P,100)