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
		using DifferentialEquations
		using Distributions
		using LightGraphs
		using LinearAlgebra
		using Random
		Random.seed!(42)
end

@everywhere begin
	dir = @__DIR__
	#include("$dir/exp_base.jl")
	include("$dir/src/experiments.jl")
#	include("$dir/input_data/demand_curves.jl")
	include("$dir/src/network_dynamics.jl")
end

@everywhere begin
	N = 24
	num_days =  50
	batch_size = 100 #2
	lambda_lst = 0.2:0.1:1 # 0.25:0.25:1#
	num_monte = batch_size * length(lambda_lst)
end

@everywhere begin
	freq_threshold = 0.002
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
#            SIM 3                  #
########################################


vc3 = independent_set(graph_lst[1], DegreeIndependentSet()) # ilc_nodes
cover3 = Dict([v => [] for v in vc3]) # ilc_cover
lambda = lambda_lst[1]

_compound_pars = experiments.compound_pars(N, low_layer_control, kappa, vc3, cover3, lambda)

_compound_pars.hl.daily_background_power .= 0.
_compound_pars.hl.current_background_power .= 0.
_compound_pars.hl.mismatch_yesterday .= 0.
_compound_pars.periodic_demand =  periodic_demand
_compound_pars.residual_demand = residual_demand
_compound_pars.graph = graph_lst[1]

@everywhere compound_pars = $_compound_pars

# @everywhere sign_cb = ContinuousCallback(network_dynamics.SignChangeCondition,network_dynamics.SignChangeSave!)

@everywhere begin
	ic = zeros(compound_pars.D * compound_pars.N)
	tspan = (0., num_days * l_day)
	ode_tlIII = ODEProblem(network_dynamics.ACtoymodel!, ic, tspan, compound_pars,
		callback=CallbackSet(PeriodicCallback(network_dynamics.HourlyUpdateEcon(), l_hour, initial_affect= false),
						 PeriodicCallback(network_dynamics.DailyUpdate_III, l_day, initial_affect= false)))
end

# HourlyUpdateEcon(verbose=true)
# sol_test = solve(ode_tlIII, Rodas4P())
# using Plots
# plot(sol_test, vars=energy_filter)
# println("bla")
# savefig("$dir/plots/energy_kappa_0-075.png")

# plot(sol_test.t,t->periodic_demand(t)[1] + residual_demand(t)[1])

monte_probIII = EnsembleProblem(
	ode_tlIII,
	output_func = (sol, i) -> experiments.observer_basic_types_lambda_III(sol, i, freq_filter, energy_filter,  energy_abs_filter, freq_threshold, obs_days),
	prob_func = (prob,i,repeat) -> experiments.prob_func_III(prob, i, repeat, batch_size, graph_lst, num_days, lambda_lst),
	u_init = [])

res_tlIII = @time solve(monte_probIII,
					 Rodas4P(),
					 trajectories=num_monte,
					 batch_size=batch_size, EnsembleDistributed())

using Dates
date = Dates.Date(Dates.now())

if isdir("$dir/solutions/$(date)") == false
	mkdir("$dir/solutions/$(date)")
end

jldopen("$dir/solutions/$(date)/expIII_sol_new.jld2", true, true, true, IOStream) do file
	file["monte_probIII"] = monte_probIII
	file["res_tlIII"] = res_tlIII
end

graph_dict = Dict("g$i"=> Graph(res_tlIII.u[i][5]) for i = 1:length(res_tlIII))
savegraph("$dir/solutions/$(date)/expIII_graph_lst.lg", graph_dict)



# only debugging hourlyupdate
function fhu(time, lambda, sol)
	sum_lambda = 3600 * (lambda - 1)
	sum_sign = zeros(N)
	for t in time-3600+1:time
		sum_sign += sign.(sol(t)[energy_filter])
	end
	sign_term = lambda .* sum_sign .+ sum_lambda
	my = zeros(24,N)
	my[Int((time-3600)/3600)+1,:] .= ones(N) .*  sign_term # integrator.u[sign_idx]
end
