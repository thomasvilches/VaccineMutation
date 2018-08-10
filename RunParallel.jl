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
    resultsR02 = Vector{Int64}(numberofsims)
    resultsPV = Matrix{Float64}(P.matrix_strain_lines,numberofsims)
    resultsPNV = Matrix{Float64}(P.matrix_strain_lines,numberofsims)
    resultsEf = Matrix{Float64}(P.matrix_strain_lines,numberofsims)
    resultsTimeInf = Matrix{Int64}(P.grid_size_human,numberofsims)
    resultsTimeRec = Matrix{Int64}(P.grid_size_human,numberofsims)
    resultsDistance = Matrix{Int64}(P.grid_size_human,numberofsims)
    for i=1:numberofsims
        resultsL[:,i] = results[i][1]
        resultsS[:,i] = results[i][2]
        resultsA[:,i] = results[i][3]
        resultsR0[i] = results[i][4]
        resultsPV[:,i] = results[i][5]
        resultsPNV[:,i] = results[i][6]
        resultsGD[:,i] = results[i][7]
        resultsEf[:,i] = results[i][8]
        resultsTimeInf[:,i] = results[i][9]
        resultsDistance[:,i] = results[i][10]
        resultsTimeRec[:,i] = results[i][11]
        resultsR02[i] = results[i][12]
    end
   
    directory = "Aug03/results1/"

    writedlm(string("$directory","result","$(P.Prob_transmission)","Mut","$(P.mutation_rate)","Ef","$(P.VaccineEfficacy)","_latent.dat"),resultsL)
    writedlm(string("$directory","result","$(P.Prob_transmission)","Mut","$(P.mutation_rate)","Ef","$(P.VaccineEfficacy)","_symp.dat"),resultsS)
    writedlm(string("$directory","result","$(P.Prob_transmission)","Mut","$(P.mutation_rate)","Ef","$(P.VaccineEfficacy)","_asymp.dat"),resultsA)
    writedlm(string("$directory","result","$(P.Prob_transmission)","Mut","$(P.mutation_rate)","Ef","$(P.VaccineEfficacy)","_R0.dat"),resultsR0)
    writedlm(string("$directory","result","$(P.Prob_transmission)","Mut","$(P.mutation_rate)","Ef","$(P.VaccineEfficacy)","_PV.dat"),resultsPV)
    writedlm(string("$directory","result","$(P.Prob_transmission)","Mut","$(P.mutation_rate)","Ef","$(P.VaccineEfficacy)","_PNV.dat"),resultsPNV)
    writedlm(string("$directory","result","$(P.Prob_transmission)","Mut","$(P.mutation_rate)","Ef","$(P.VaccineEfficacy)","_GD.dat"),resultsGD)
    writedlm(string("$directory","result","$(P.Prob_transmission)","Mut","$(P.mutation_rate)","Ef","$(P.VaccineEfficacy)","_Ef.dat"),resultsEf)
    writedlm(string("$directory","result","$(P.Prob_transmission)","Mut","$(P.mutation_rate)","Ef","$(P.VaccineEfficacy)","_TimeInf.dat"),resultsTimeInf)
    writedlm(string("$directory","result","$(P.Prob_transmission)","Mut","$(P.mutation_rate)","Ef","$(P.VaccineEfficacy)","_Distance.dat"),resultsDistance)
    writedlm(string("$directory","result","$(P.Prob_transmission)","Mut","$(P.mutation_rate)","Ef","$(P.VaccineEfficacy)","_TimeRec.dat"),resultsTimeRec)
    writedlm(string("$directory","result","$(P.Prob_transmission)","Mut","$(P.mutation_rate)","Ef","$(P.VaccineEfficacy)","_R02.dat"),resultsR02)

end

function run_main(P::InfluenzaParameters,numberofsims::Int64)
    
    results = pmap((cb, x) -> main(cb, x, P), Progress(numberofsims*P.sim_time), 1:numberofsims, passcallback=true)

    dataprocess(results,P,numberofsims)
end


@everywhere P=InfluenzaParameters(
    VaccineEfficacy = 0.0,
    GeneralCoverage = 0,
    Prob_transmission = 0.015,
    sim_time = 200,
    grid_size_human = 1000,
    matrix_strain_lines = 1200,
    mutation_rate = 0.3,
    initial_p = 0.01,
    initial_p2 = 0.005,
    initial_p3 = 0.015,
    start_different = 1
)

run_main(P,1000)