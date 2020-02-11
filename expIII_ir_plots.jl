using JLD2, FileIO, GraphIO
using Distributed
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
	#num_days =  30
	obs_days=10
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

function index_calc(i,batch_size)
	out = []
	if i <= 0 || batch_size <= 0
		error("Inputs <= 0 are not allowed.")
	elseif i%batch_size == 0
		out = batch_size
	else
		out = i%batch_size
	end
	out
end

using Dates
date = Dates.today()  - Dates.Day(2)

@load "$dir/solutions/$(date)/expIII_sol_ir_lambda_0.8.jld2" res_tlIII_ir monte_probIII_ir
#res_tlIII_ir = res_tlIV_ir

lamda = monte_probIII_ir.prob.p.hl.lambda

date = Dates.today() - Dates.Day(1)
if isdir("$dir/plots/$(date)") == false
	mkdir("$dir/plots/$(date)")
end

using Plots
using DataFrames
using StatsPlots

Plots.scalefontsizes(2)

######################################
######## BOX #####################
######################################
ilc_lst = 1:N

@time include("expIII_ir_ex_plots.jl")
###############################################
@time include("expIII_ir_max_freq_plots.jl")
###########################################
@time include("expIII_ir_control_energy_plots.jl")

@time include("expIII_ir_cost_plots.jl")
