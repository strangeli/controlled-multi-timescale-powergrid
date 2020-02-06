using JLD2, FileIO, GraphIO
using Distributed
using Interpolations

_calc = true
_slurm = true

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
	num_days =  50
	batch_size = 100
	lambda_lst = 0.2:0.1:1
	num_monte = batch_size * length(lambda_lst)
end

@everywhere begin
	freq_threshold = 0.0005
	obs_days = 10
	phase_filter = 1:N
	freq_filter = N+1:2N
	control_filter = 2N+1:3N
	energy_filter = 3N+1:4N
	energy_abs_filter = 4N+1:5N
end


############################################

@everywhere begin
	l_day = 3600*24 # DemCurve.l_day
	l_hour = 3600 # DemCurve.l_hour
	l_minute = 60 # DemCurve.l_minute
	low_layer_control = experiments.LeakyIntegratorPars(M_inv=0.2,kP=525,T_inv=1/0.05,kI=0.005)
	kappa = 0.15 / l_hour
end

############################################
# this should only run on one process
############################################

_graph_lst = []
for i in 1:num_monte
	push!(_graph_lst, random_regular_graph(iseven(3N) ? N : (N-1), 3))
end
@everywhere graph_lst = $_graph_lst

# graph_lst_orig = Dict("g$i"=> graph_lst[i] for i = 1:num_monte)
# savegraph("$dir/solutions/20191121/graph_lst_orig.lg", graph_lst_orig)


# graph_lst = [collect(values(loadgraphs("$dir/solutions/20191121/graph.lg")))]
# fault_size = 3
# edges_to_delete = sample(edges(graph_lst[1]) |> collect, fault_size, replace=false)
# for e in edges_to_delete
# 	rem_edge!(graph_lst[1], e)
# end
# savegraph("$dir/solutions/20191121/graph_faulty.lg", graph_lst[1])

############################################
#  demand
############################################

demand_amp = rand(N)
periodic_demand = t -> demand_amp .* sin(t*pi/(24*3600))^2
samples = 24*4
inter = interpolate([.2 * randn(N) for i in 1:(num_days * samples + 1)], BSpline(Linear()))
residual_demand = t -> inter(1. + t / (24*3600) * samples) # 1. + is needed to avoid trying to access out of range


#######################################
#            SIM 2                    #
########################################

vc2 = independent_set(graph_lst[1], DegreeIndependentSet()) # ilc_nodes
cover2 = Dict([v => neighbors(graph_lst[1], v) for v in vc2]) # ilc_cover
lambda = lambda_lst[1]

_compound_pars = experiments.compound_pars(N, low_layer_control, kappa, vc2, cover2, lambda)
_compound_pars.hl.daily_background_power .= 0.
_compound_pars.hl.current_background_power .= 0.
_compound_pars.hl.mismatch_yesterday .= 0.
_compound_pars.periodic_demand = periodic_demand
_compound_pars.residual_demand = residual_demand
_compound_pars.graph = graph_lst[1]

@everywhere compound_pars = $_compound_pars


@everywhere begin
	ic = zeros(compound_pars.D * compound_pars.N)
	tspan = (0., num_days * l_day)
	ode_tlII = ODEProblem(network_dynamics.ACtoymodel!, ic, tspan, compound_pars,
		callback=CallbackSet(PeriodicCallback(network_dynamics.HourlyUpdateEcon(), l_hour, initial_affect= false),
						 PeriodicCallback(network_dynamics.DailyUpdate_II, l_day, initial_affect= false)))
end

#sol_test = solve(ode_tlII, Rodas4P())

monte_probII = EnsembleProblem(
	ode_tlII,
	output_func = (sol, i) -> experiments.observer_basic_types_lambda_II(sol, i, freq_filter, energy_filter,  energy_abs_filter, freq_threshold, obs_days),
	prob_func = (prob,i,repeat) -> experiments.prob_func_II(prob,i,repeat, batch_size, graph_lst, num_days, lambda_lst),
	u_init = [])

res_tlII = @time solve(monte_probII,
					 Rodas4P(),
					 trajectories=num_monte,
					 batch_size=batch_size, EnsembleDistributed())

#erg = res_tl.u
#@save "$dir/solutions/exp3_res.jld2" erg
#@save "$dir/solutions/exp3_prob.jld2" monte_prob
using Dates
date = Dates.Date(Dates.now())

if isdir("$dir/solutions/$(date)") == false
 	mkdir("$dir/solutions/$(date)")
end

jldopen("$dir/solutions/$(date)/expII_sol_new.jld2", true, true, true, IOStream) do file
	file["monte_probII"] = monte_probII
	file["res_tlII"] = res_tlII
end

graph_dict = Dict("g$i"=> Graph(res_tlII.u[i][5]) for i = 1:length(res_tlII))
savegraph("$dir/solutions/$(date)/expII_graph_lst.lg", graph_dict)

#@save "$dir/solutions/exp3_sol.jld2" res_tl monte_prob
