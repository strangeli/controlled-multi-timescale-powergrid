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
	num_days =  30
	batch_size = 100
	lambda_lst = 0.2:0.1:1
	ilc_lst = 1:24
	num_monte = batch_size * length(ilc_lst)
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


##########################################


ilc_ratio = ilc_lst[1]#Int(N/2)
vc4 = sample(1:N,ilc_ratio,replace = false) # ilc_nodes
cover4 = [] # ilc_covers

# exp IV communication
# kappa_factor = kappa
# kappa4 =  diagm(0 => zeros(N))
#cover = [] # ilc_covers
# for i in 1:ilc_ratio
# 	kappa4[vc4[i],vc4[i]] = 1
# 	#kappa[few[1], vc[i]] = 1 # here multiple entries
# 	#kappa[few[2], vc[i]] = 1
# 	#kappa[few[3], vc[i]] = 1
# 	a = 1:vc4[i]-1
# 	b = vc4[i]+1:N
# 	c = [collect(a); collect(b)]
# 	few = sample(c,3,replace = false)
# 	kappa4[vc4[i], few[1]] = 1
# 	kappa4[vc4[i], few[2]] = 1
# 	kappa4[vc4[i], few[3]] = 1
# 	push!(cover4, Dict([vc4[i] => few]))
# end
# kappa4 = kappa_factor .* kappa4 .* 0.25
lambda = 0.75

_compound_pars = experiments.compound_pars(N, low_layer_control, kappa, vc4, cover4, lambda)

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
	ode_tlIII_ir = ODEProblem(network_dynamics.ACtoymodel!, ic, tspan, compound_pars,
		callback=CallbackSet(PeriodicCallback(network_dynamics.HourlyUpdateEcon(), l_hour, initial_affect= false),
						 PeriodicCallback(network_dynamics.DailyUpdate_III, l_day, initial_affect= false)))
end


# sol2 = solve(ode_tlIII_ir, Rodas4P())
# using Plots
# plot(sol2,vars=energy_filter, legend=false)

monte_probIII_ir = EnsembleProblem(
	ode_tlIII_ir,
	output_func = (sol, i) -> experiments.observer_basic_types_lambda_III(sol, i, freq_filter, energy_filter,  energy_abs_filter, freq_threshold, obs_days),
	prob_func = (prob,i,repeat) -> experiments.prob_func_III_ir(prob,i,repeat, batch_size, graph_lst, num_days, ilc_lst),
	u_init = [])

res_tlIII_ir = @time solve(monte_probIII_ir,
					 Rodas4P(),
					 trajectories=num_monte,
					 batch_size=batch_size, EnsembleDistributed())

using Dates
date = Dates.today()

if isdir("$dir/solutions/$(date)") == false
	mkdir("$dir/solutions/$(date)")
end

jldopen("$dir/solutions/$(date)/expIII_sol_ir_lambda_$(lambda).jld2", true, true, true, IOStream) do file
	file["monte_probIII_ir"] = monte_probIII_ir
	file["res_tlIII_ir"] = res_tlIII_ir
end

graph_dict = Dict("g$i"=> Graph(res_tlIII_ir.u[i][5]) for i = 1:length(res_tlIII_ir))
savegraph("$dir/solutions/$(date)/expIII_ir_graph_lst_$(lambda).lg", graph_dict)
