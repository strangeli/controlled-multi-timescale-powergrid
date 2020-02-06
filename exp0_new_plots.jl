using JLD2, FileIO, GraphIO
using Distributed
using Interpolations

_calc = false
_slurm = false

if _calc
    using ClusterManagers
	if length(ARGS) > 0
		N_tasks = parse(Int, ARGS[1])
	else
		N_tasks = 1
	end
    N_worker = N_tasks
	if _slurm
    	addprocs(SlurmManager(N_worker))
	else
		addprocs(N_worker)
	end
	println()
	println(nprocs(), " processes")
	println(length(workers()), " workers")
else
	using Plots
end

# here comes the broadcast
# https://docs.julialang.org/en/v1/stdlib/Distributed/index.html#Distributed.@everywhere
@everywhere begin
	calc = $_calc # if false, only plotting
end

@everywhere begin
	dir = @__DIR__
	#include("$dir/exp_base.jl")
	include("$dir/src/experiments.jl")
#	include("$dir/input_data/demand_curves.jl")
	include("$dir/src/network_dynamics.jl")
end

@everywhere begin
		using DifferentialEquations
		using Distributions
		using LightGraphs
		using LinearAlgebra
		using Random
		Random.seed!(42)
end

@everywhere begin
	N = 24
	num_days =  5
	num_monte = 100
	batch_size = 100
end

@everywhere begin
	freq_threshold = 0.2
	phase_filter = 1:N
	freq_filter = N+1:2N
	control_filter = 2N+1:3N
	energy_filter = 3N+1:4N
	energy_abs_filter = 4N+1:5N
end

@load "$dir/solutions_only_MC_no_time_series/20190909_cluster/exp0_sol_new.jld2" res_tl0 monte_prob0

plot([p[1] for p in res_tl0.u],seriestype=:scatter)
xlabel!("i")
ylabel!("max. freq. deviation")
savefig("$dir/plots/20190909/exp0_new/exp0_new_freq_dev.png")

plot([p[2] for p in res_tl0.u],seriestype=:scatter)
xlabel!("i")
ylabel!("exceedance")
savefig("$dir/plots/20190909/exp0_new/exp0_new_exceedance.png")

plot(hcat(Vector([p[3] for p in res_tl0.u])...)',seriestype=:scatter)
xlabel!("i")
ylabel!("control energy")
savefig("$dir/plots/20190909/exp0_new/exp0_new_control_energy.png")

plot(hcat(Vector([p[4] for p in res_tl0.u])...)',seriestype=:scatter)
xlabel!("i")
ylabel!("frequency variance")
savefig("$dir/plots/20190909/exp0_new/exp0_new_freq_var.png")
