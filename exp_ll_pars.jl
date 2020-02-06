using JLD2, FileIO
using Distributed

_test = false # false for exp run on cluster
_calc = true
_slurm =true

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
	println(nprocs(), " processes")
	println(length(workers()), " workers")
else
	using Plots
end

# here comes the broadcast
# https://docs.julialang.org/en/v1/stdlib/Distributed/index.html#Distributed.@everywhere
@everywhere begin
	calc = $_calc # if false, only plotting
	test = $_test
end

@everywhere begin
	dir = @__DIR__
	include("$dir/src/experiments.jl")
	include("$dir/src/network_dynamics.jl")
end

@everywhere begin
		using DifferentialEquations
		using Distributions
		using LightGraphs
		using LinearAlgebra
		using Random
		using Interpolations
		Random.seed!(42)
end


@everywhere if test
	N = 4 # N must be even for the random regular graph with k=3
	num_days = 1
	num_runs = 2
else
	N = 24
	num_days =  1
	num_runs = 100
end

@everywhere begin
	batch_size = num_runs
	freq_threshold = 0.0005
	phase_filter = 1:N
	freq_filter = N+1:2N
	control_filter = 2N+1:3N
	energy_filter = 3N+1:4N
	energy_abs_filter = 4N+1:5N
end

############################################
# Parameters
############################################

@everywhere begin
	low_layer_control = experiments.LeakyIntegratorPars(M_inv=0.2,kP=525,T_inv=1/0.05,kI=0.005)
	kappa= 0.15/experiments.l_hour
	vc = [] # ilc_nodes (here: without communication)
	cover = Dict([v => [] for v in vc])# ilc_cover
	lambda = 1
	obs_days=1
#	compound_pars = experiments.compound_pars(N, low_layer_control, kappa, vc, cover)
end

############################################
# ODE Problem
############################################

# add modifications from exp_base here
@everywhere begin
	if test
		kP_lst_s = 0:250:750
		kI_lst_s = 0.001:0.25:1
		kP_lst = repeat(kP_lst_s; inner=length(kI_lst_s))
		kI_lst = repeat(kI_lst_s; outer=length(kP_lst_s))
	else
		kP_lst_s = 0:25:1000
		kI_lst_s = 0.001:0.025:1
		kP_lst = repeat(kP_lst_s; inner=length(kI_lst_s))
		kI_lst = repeat(kI_lst_s; outer=length(kP_lst_s))
	end
	num_batches = length(kP_lst)
	num_monte = batch_size * num_batches
end


_graph_lst = []
for i in 1:num_monte
	push!(_graph_lst, random_regular_graph(iseven(3N) ? N : (N-1), 3))
end
@everywhere graph_lst = $_graph_lst

begin
	demand_amp = -rand(N) #.* 100.
	# demand_ramp = rand(N) .* 2.
	periodic_demand = t -> demand_amp .* sin(t*pi/(24*3600))^2
	samples = 24*4
	inter = interpolate([0.2 * randn(N) for i in 1:(num_days * samples + 1)], BSpline(Linear()))
	residual_demand = t -> inter(1. + t / (24*3600) * samples)
end

lambda = 0.8

begin
	_compound_pars = experiments.compound_pars(N, low_layer_control, kappa, vc, cover, lambda)
	_compound_pars.hl.daily_background_power .= 0.
	_compound_pars.hl.current_background_power .= 0.
	_compound_pars.hl.mismatch_yesterday .= 0.
	_compound_pars.periodic_demand = periodic_demand
	_compound_pars.residual_demand = residual_demand
	_compound_pars.graph = graph_lst[1]
end

@everywhere compound_pars = $_compound_pars

@everywhere begin
	ic = zeros(compound_pars.D * compound_pars.N)
	tspan_ll = tspan = (0., 0.1 * num_days * experiments.l_day)
	ode_ll = ODEProblem(network_dynamics.ACtoymodel!, ic, tspan_ll, compound_pars)
end


# @time solve(ode_ll, Tsit5());
# 	176.763785 seconds (44.07 M allocations: 83.515 GiB, 22.40% gc time)

############################################
# Monte Carlo
############################################

if calc
	monte_prob = EnsembleProblem(ode_ll,
	    output_func = (sol,i) -> experiments.observer_basic_types(sol, i, freq_filter, energy_filter, freq_threshold),
	    prob_func = (prob, i, repeat) -> experiments.prob_func(prob, i, repeat, batch_size, kP_lst,kI_lst,graph_lst,num_days),
	    reduction = (u, data, I) -> experiments.reduction(u, data, I, batch_size),
	    u_init = [])

	res = @time solve(monte_prob,
	                     Rodas4P(),
	                     #saveat=1,
	                     #save_idxs=1:nodes,
	                     trajectories=num_monte,
	                     batch_size=batch_size, EnsembleDistributed())

						 jldopen("$dir/solutions/20200120_exp_ll_pars.jld2", true, true, true, IOStream) do file
					 		file["monte_prob"] = monte_prob
					 		file["res"] = res
					 	end
else
	using Plots
	@load "$dir/solutions/cluster/20190604_exp1a_sol_100nodes.jld2" res monte_prob
	#factor = maximum(abs, res.u, dims=1)
	plot(collect(kP_lst), res.u[:,1], label=["omega_max"])
	xlabel!("k_P")
	ylabel!("omega max")
	savefig("$dir/plots/exp1a_freq_dev_100nodes.png")
	plot(collect(kP_lst), res.u[:,2], label=["exceedance"])
	xlabel!("k_P")
	ylabel!("exceedance")
	savefig("$dir/plots/exp1a_exceedance_100nodes.png")
	plot(collect(kP_lst), res.u[:,3], label=["control energy"])
	xlabel!("k_P")
	ylabel!("control energy")
	savefig("$dir/plots/exp1a_control_energy_100nodes.png")
	plot(collect(kP_lst), res.u[:,4], label=["ex_std"])
	xlabel!("k_P")
	ylabel!("exceedance std")
	savefig("$dir/plots/exp1a_ex_std_100nodes.png")
	plot(collect(kP_lst),res.u[:,5], label=["energy_std"])
	xlabel!("k_P")
	ylabel!("energy std")
	savefig("$dir/plots/exp1a_energy_std_100nodes.png")
end
