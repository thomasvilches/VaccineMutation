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
    
    resultsPV = Matrix{Float64}(P.matrix_strain_lines,numberofsims)
    resultsPNV = Matrix{Float64}(P.matrix_strain_lines,numberofsims)
    resultsEf = Matrix{Float64}(P.matrix_strain_lines,numberofsims)
    for i=1:numberofsims
        resultsL[:,i] = results[i][1]
        resultsS[:,i] = results[i][2]
        resultsA[:,i] = results[i][3]
        resultsR0[i] = results[i][4]
        resultsPV[:,i] = results[i][5]
        resultsPNV[:,i] = results[i][6]
        resultsGD[:,i] = results[i][7]
        resultsEf[:,i] = results[i][8]
    end
   
    directory = "Aug02/results/"

    writedlm(string("$directory","result","$(P.Prob_transmission)","Mut","$(P.mutation_rate)","Ef","$(P.VaccineEfficacy)","_latent.dat"),resultsL)
    writedlm(string("$directory","result","$(P.Prob_transmission)","Mut","$(P.mutation_rate)","Ef","$(P.VaccineEfficacy)","_symp.dat"),resultsS)
    writedlm(string("$directory","result","$(P.Prob_transmission)","Mut","$(P.mutation_rate)","Ef","$(P.VaccineEfficacy)","_asymp.dat"),resultsA)
    writedlm(string("$directory","result","$(P.Prob_transmission)","Mut","$(P.mutation_rate)","Ef","$(P.VaccineEfficacy)","_R0.dat"),resultsR0)
    writedlm(string("$directory","result","$(P.Prob_transmission)","Mut","$(P.mutation_rate)","Ef","$(P.VaccineEfficacy)","_PV.dat"),resultsPV)
    writedlm(string("$directory","result","$(P.Prob_transmission)","Mut","$(P.mutation_rate)","Ef","$(P.VaccineEfficacy)","_PNV.dat"),resultsPNV)
    writedlm(string("$directory","result","$(P.Prob_transmission)","Mut","$(P.mutation_rate)","Ef","$(P.VaccineEfficacy)","_GD.dat"),resultsGD)
    writedlm(string("$directory","result","$(P.Prob_transmission)","Mut","$(P.mutation_rate)","Ef","$(P.VaccineEfficacy)","_Ef.dat"),resultsEf)
end

function run_main(P::InfluenzaParameters,numberofsims::Int64)
    
    results = pmap((cb, x) -> main(cb, x, P), Progress(numberofsims*P.sim_time), 1:numberofsims, passcallback=true)

    dataprocess(results,P,numberofsims)
end


@everywhere P=InfluenzaParameters(
    VaccineEfficacy = 0.8,
    GeneralCoverage = 1,
    Prob_transmission = 0.079,
    sim_time = 200,
    grid_size_human = 1000,
    matrix_strain_lines = 1200,
    mutation_rate = 0.3,
    initial_p = 0.01,
    initial_p2 = 0.005,
    start_different = 1
)

run_main(P,1000)