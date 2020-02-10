""" Run this file to generate plots from the solutions """

using JLD2, FileIO, GraphIO
using Interpolations



begin
	dir = @__DIR__
	include("$dir/src/experiments.jl")
	include("$dir/src/network_dynamics.jl")
end

begin
		using DifferentialEquations
		using Distributions
		using LightGraphs
		using LinearAlgebra
		using Random
		Random.seed!(42)
end

begin
	N = 24
	#num_days =  50
	obs_days = 10
	#num_monte = 100
	batch_size = 100
end

begin
	freq_threshold = 0.001
	phase_filter = 1:N
	freq_filter = N+1:2N
	control_filter = 2N+1:3N
	energy_filter = 3N+1:4N
	energy_abs_filter = 4N+1:5N
end


using Dates
date = Dates.today() - Dates.Day(1)

@load "$dir/solutions/$(date)/exp0_sol_new.jld2" res_tl0 monte_prob0
@load "$dir/solutions/$(date)/expI_sol_new.jld2" res_tlI monte_probI
@load "$dir/solutions/$(date)/expII_sol_new.jld2" res_tlII monte_probII
@load "$dir/solutions/$(date)/expIII_sol_new.jld2" res_tlIII monte_probIII
@load "$dir/solutions/$(date)/expIV_sol_new.jld2" res_tlIV monte_probIV

#date = Dates.today()

if isdir("$dir/plots/$(date)") == false
	mkdir("$dir/plots/$(date)")
end

using Plots
using DataFrames
using StatsPlots
using LaTeXStrings

Plots.scalefontsizes(2)

######################################
######## BOX #####################


function index_calc(i,batch_size)
	out = []
	if i <= 0 || batch_size <= 0
		println("error: values <= 0 are not allowed.")
	elseif i%batch_size == 0
		out = batch_size
	else
		out = i%batch_size
	end
	out
end

######################################
lambda_lst = 0.2:0.1:1
lam = 0

@time include("exp0-IV_ex_plots.jl")

###############################################
@time include("exp0-IV_max_freq_plots.jl")
###########################################
@time include("exp0-IV_control_energy_plots.jl")

@time include("exp0-IV_cost_plots.jl")

############################################
# By experiment

include("plot_by_exp.jl")
