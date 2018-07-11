include("SlurmConnect.jl")

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
using SlurmConnect 

add_truck(LumberjackTruck("processrun.log"), "my-file-logger")
remove_truck("console")
info("lumberjack process started up, starting repl")

info("adding procs...")

s = SlurmManager(512)
@eval Base.Distributed import Base.warn_once
addprocs(s, partition="defq", N=16)

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
   
    directory = "July10/"

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

for Vef = 0.2:0.05:0.8
    P=InfluenzaParameters(
        VaccineEfficacy = $Vef,
        GeneralCoverage = 1,
        Prob_transmission = 0.079,
        sim_time = 300,
        grid_size_human = 10000,
        mutation_rate = 0.00416
    )

    run_main(P,1000)

end