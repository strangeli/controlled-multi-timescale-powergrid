@doc """
This is a module that contains functions with specified setups for numerical
experiments.
It depends on the modules `network_dynamics.jl` and `observables.jl`

# Examples
```julia-repl
julia> include("src/experiments.jl")
Main.experiments
```
"""
module experiments

using Random # fix random seed
using DifferentialEquations # problems types
using Plots # custom plot functions
using LightGraphs # create incidence matrices
using Parameters
using LinearAlgebra
using Interpolations
using StatsBase
using Statistics
using Distributions

begin
	# import dynamics and observables
	dir = @__DIR__
	include("$dir/network_dynamics.jl")
	include("$dir/observables.jl")
end

l_hour = 60 * 60 # in s
l_day = l_hour * 24 # in s

@doc """
    LeakyIntegratorPars(M_inv, kP, T_inv, kI)
Define a parameter struct.
"""
@with_kw mutable struct LeakyIntegratorPars
    M_inv
    kP
    T_inv
	kI
end

@with_kw mutable struct ILCPars
	kappa
	mismatch_yesterday::Array{Float64, 2}
	daily_background_power::Array{Float64, 2} # 24xN vector with the background power for each hour.
	current_background_power
	ilc_nodes
	ilc_covers
	lambda
	sign_save
end



@doc """
    UlmoPars(N, D, ll. hl, periodic_infeed, periodic_demand, fluctuating_infeed, residual_demand, incidence, coupling)
Define a parameter struct.
"""
@with_kw mutable struct UlMoPars #<: DEParameters # required parameter supertype for MCBB
	N::Int # number of nodes
	D::Int # degrees of freedom of each node
    ll::LeakyIntegratorPars # low level control parameters
	hl:: ILCPars			# high level control parameters
	periodic_infeed
	periodic_demand
	fluctuating_infeed
	residual_demand
	incidence
	coupling
	graph
	"""
	Constructor
	"""
	function UlMoPars(N::Int, #D::Int, # degrees of freedom of each node
    					ll::LeakyIntegratorPars,
						hl:: ILCPars,
						periodic_infeed,
						periodic_demand,
						fluctuating_infeed,
						residual_demand,
						#incidence,
						coupling,
						graph)
			new(N, 5,
			ll,
			hl,
			periodic_infeed,
			periodic_demand,
			fluctuating_infeed,
			residual_demand,
			incidence_matrix(graph,oriented=true),
			coupling,
			graph)
	end
end

modx(x, y) = (mod2pi.(x .+ π/2) .- π/2, y)
mody(x, y) = (x, mod2pi.(y .+ π/2) .- π/2)



@doc """
    default_pars(N)
Setup the system with default parameters.
"""
function default_pars(N)
	low_layer_control = LeakyIntegratorPars(M_inv=60.,kP=1.3,T_inv=1.,kI=0.9)
	g = random_regular_graph(iseven(3N) ? N : (N-1), 3)
	#incidence = incidence_matrix(g, oriented=true)
	coupling= 6. .* diagm(0=>ones(size(incidence_matrix(g, oriented=true),2)))
	vc = independent_set(g, DegreeIndependentSet()) # ilc_nodes
	cover = Dict([v => neighbors(g, v) for v in vc]) # ilc_cover
	lambda = 1
	sign_save = zeros(1,N+1)

	higher_layer_control = ILCPars(kappa=0.35, mismatch_yesterday=zeros(24, N), daily_background_power=zeros(24, N), current_background_power=zeros(N), ilc_nodes=vc, ilc_covers=cover, lambda = lambda, sign_save= sign_save)
	periodic_infeed = t -> zeros(N)
	peak_demand = rand(N)
	periodic_demand= t -> peak_demand .* abs(sin(pi * t/24.))
	fluctuating_infeed = t -> zeros(N)
	residual_demand= t -> zeros(N)



	return UlMoPars(N, low_layer_control,
							higher_layer_control,
							periodic_infeed,
							periodic_demand,
							fluctuating_infeed,
							residual_demand,
							#incidence,
							coupling,
							g)
end

function compound_pars(N, low_layer_control, kappa, ilc_nodes, ilc_covers, lambda)
	sign_save = zeros(1,N+1)
	higher_layer_control = ILCPars(kappa=kappa, mismatch_yesterday=zeros(24, N), daily_background_power=zeros(24, N), current_background_power=zeros(N),ilc_nodes=ilc_nodes, ilc_covers=ilc_covers, lambda = lambda, sign_save=sign_save)

	periodic_infeed = t -> zeros(N)
	periodic_demand= t -> zeros(N)
	fluctuating_infeed = t -> zeros(N)
	residual_demand = t -> zeros(N)

	# make sure N*k is even, otherwise the graph constructor fails
	g = random_regular_graph(iseven(3N) ? N : (N-1), 3)
	#incidence = incidence_matrix(g, oriented=true)
	coupling= 6. .* diagm(0=>ones(ne(g)))

	return UlMoPars(N, low_layer_control,
							higher_layer_control,
							periodic_infeed,
							periodic_demand,
							fluctuating_infeed,
							residual_demand,
							#incidence,
							coupling,
							g)
end

##############################################################################
# Monte Carlo functions

get_run(i, batch_size) = mod(i, batch_size)==0 ? batch_size : mod(i, batch_size)
get_batch(i, batch_size) = 1 + (i - 1) ÷ batch_size



function prob_func(prob, i, repeat, batch_size, kP_lst, kI_lst, graph_lst, num_days)
	println("sim ", i)
    	run = get_run(i, batch_size)
    	batch = get_batch(i, batch_size)
    	demand_amp = - rand(prob.p.N) #.* 100.
	prob.p.periodic_demand = t -> demand_amp .* sin(t*pi/(24*3600))^2
	samples = 24*4
	inter = interpolate([0.2 * randn(prob.p.N) for i in 1:(num_days * samples + 1)], BSpline(Linear()))
	prob.p.residual_demand = t -> inter(1. + t / (24*3600) * samples) # 1. + is needed to avoid trying to access out of range

	prob.p.hl.daily_background_power .= 0.
	prob.p.hl.current_background_power .= 0.
	prob.p.hl.mismatch_yesterday .= 0.
	# Random network in each sim
	prob.p.graph = graph_lst[i]
	prob.p.incidence = incidence_matrix(graph_lst[i], oriented=true)

	prob.p.ll.kP = kP_lst[batch]
	prob.p.ll.kI = kI_lst[batch]
	ODEProblem(network_dynamics.ACtoymodel!, prob.u0, prob.tspan, prob.p)

    prob
end

function prob_func_0(prob, i, repeat, batch_size, graph_lst, num_days, lambda_lst)
	println("sim ", i)
	run = get_run(i, batch_size)
    batch = get_batch(i, batch_size)

	demand_amp = rand(prob.p.N)
	prob.p.periodic_demand = t -> demand_amp .* sin(t*pi/(24*3600))^2
	samples = 24*4
	inter = interpolate([.2 * randn(prob.p.N) for i in 1:(num_days * samples + 1)], BSpline(Linear()))
	prob.p.residual_demand = t -> inter(1. + t / (24*3600) * samples) # 1. + is needed to avoid trying to access out of range

	prob.p.hl.daily_background_power .= 0.
	# T = 3600
	# for i in 1:24
	# 	prob.p.hl.daily_background_power[i,:] .= sum(prob.p.periodic_demand.((i-1)*T:i*T))/T
	# end
	prob.p.hl.current_background_power .= 0.
	prob.p.hl.mismatch_yesterday .= 0.
	# Random network in each sim
	prob.p.graph = graph_lst[i]
	prob.p.incidence = incidence_matrix(graph_lst[i], oriented=true)
	prob.p.hl.lambda = lambda_lst[batch]

	ODEProblem(network_dynamics.ACtoymodel!, prob.u0, prob.tspan, prob.p)
end



function prob_func_I(prob, i, repeat, batch_size, graph_lst, num_days, lambda_lst)
	println("sim ", i)
	run = get_run(i, batch_size)
    batch = get_batch(i, batch_size)

	demand_amp = rand(prob.p.N)
	prob.p.periodic_demand = t -> demand_amp .* sin(t*pi/(24*3600))^2
	samples = 24*4
	inter = interpolate([.2 * randn(prob.p.N) for i in 1:(num_days * samples + 1)], BSpline(Linear()))
	prob.p.residual_demand  =  t -> inter(1. + t / (24*3600) * samples) # 1. + is needed to avoid trying to access out of range

	prob.p.hl.daily_background_power .= 0.
	# T = 3600
	# for i in 1:24
	# 	prob.p.hl.daily_background_power[i,:] .= sum(prob.p.periodic_demand.((i-1)*T:i*T))/T
	# end
	prob.p.hl.current_background_power .= 0.
	prob.p.hl.mismatch_yesterday .= 0.
	# Random network in each sim
	prob.p.graph = graph_lst[i]
	prob.p.incidence = incidence_matrix(graph_lst[i], oriented=true)
	prob.p.hl.lambda = lambda_lst[batch]

	hourly_update = network_dynamics.HourlyUpdateEcon()

	ODEProblem(network_dynamics.ACtoymodel!, prob.u0, prob.tspan, prob.p,
		callback=CallbackSet(PeriodicCallback(hourly_update, 3600, initial_affect= false),
							 PeriodicCallback(network_dynamics.DailyUpdate, 3600*24, initial_affect= false)))
end

function prob_func_II(prob, i, repeat, batch_size, graph_lst, num_days, lambda_lst)
	println("sim ", i)
	run = get_run(i, batch_size)
    batch = get_batch(i, batch_size)

	demand_amp = rand(prob.p.N)
	prob.p.periodic_demand = t -> demand_amp .* sin(t*pi/(24*3600))^2
	samples = 24*4
	inter = interpolate([.2 * randn(prob.p.N) for i in 1:(num_days * samples + 1)], BSpline(Linear()))
	prob.p.residual_demand = t -> inter(1. + t / (24*3600) * samples) # 1. + is needed to avoid trying to access out of range

	prob.p.hl.daily_background_power .= 0.
	# T = 3600
	# for i in 1:24
	# 	prob.p.hl.daily_background_power[i,:] .= sum(prob.p.periodic_demand.((i-1)*T:i*T))/T
	# end

	prob.p.hl.current_background_power .= 0.
	prob.p.hl.mismatch_yesterday .= 0.

	# Random network in each sim
	prob.p.graph = graph_lst[i]
	prob.p.incidence = incidence_matrix(graph_lst[i], oriented=true)
	prob.p.hl.ilc_nodes = independent_set(graph_lst[i], DegreeIndependentSet()) # ilc_nodes
	prob.p.hl.ilc_covers = Dict([v => neighbors(graph_lst[i], v) for v in prob.p.hl.ilc_nodes]) # ilc_covers
	prob.p.hl.lambda = lambda_lst[batch]

	hourly_update = network_dynamics.HourlyUpdateEcon()

	ODEProblem(network_dynamics.ACtoymodel!, prob.u0, prob.tspan, prob.p,
		callback=CallbackSet(PeriodicCallback(hourly_update, 3600, initial_affect= false),
							 PeriodicCallback(network_dynamics.DailyUpdate_II, 3600*24, initial_affect= false)))
end

function prob_func_III(prob, i, repeat, batch_size, graph_lst, num_days, lambda_lst)
	println("sim ", i)
	run = get_run(i, batch_size)
    batch = get_batch(i, batch_size)

	demand_amp = rand(prob.p.N)
	prob.p.periodic_demand = t -> demand_amp .* sin(t*pi/(24*3600))^2
	samples = 24*4
	inter = interpolate([.2 * randn(prob.p.N) for i in 1:(num_days * samples + 1)], BSpline(Linear()))
	prob.p.residual_demand = t -> inter(1. + t / (24*3600) * samples) # 1. + is needed to avoid trying to access out of range

	prob.p.hl.daily_background_power .= 0.
	# T = 3600
	# for i in 1:24
	# 	prob.p.hl.daily_background_power[i,:] .= sum(prob.p.periodic_demand.((i-1)*T:i*T))/T
	# end

	prob.p.hl.current_background_power .= 0.
	prob.p.hl.mismatch_yesterday .= 0.

	# Random network in each sim
	prob.p.graph = graph_lst[i]
	prob.p.incidence = incidence_matrix(graph_lst[i], oriented=true)
	prob.p.hl.ilc_nodes = independent_set(graph_lst[i], DegreeIndependentSet()) # ilc_nodes
	prob.p.hl.lambda = lambda_lst[batch]

	hourly_update = network_dynamics.HourlyUpdateEcon()

	ODEProblem(network_dynamics.ACtoymodel!, prob.u0, prob.tspan, prob.p,
		callback=CallbackSet(PeriodicCallback(hourly_update, 3600, initial_affect= false),
							 PeriodicCallback(network_dynamics.DailyUpdate_III, 3600*24, initial_affect= false)))
end

function prob_func_III_ir(prob, i, repeat, batch_size, graph_lst, num_days, ilc_lst)
	println("sim ", i)
	run = get_run(i, batch_size)
    batch = get_batch(i, batch_size)

	demand_amp = rand(prob.p.N)
	prob.p.periodic_demand = t -> demand_amp .* sin(t*pi/(24*3600))^2
	samples = 24*4
	inter = interpolate([.2 * randn(prob.p.N) for i in 1:(num_days * samples + 1)], BSpline(Linear()))
	prob.p.residual_demand = t -> inter(1. + t / (24*3600) * samples) # 1. + is needed to avoid trying to access out of range

	prob.p.hl.daily_background_power .= 0.
	prob.p.hl.current_background_power .= 0.
	prob.p.hl.mismatch_yesterday .= 0.

	# Random network in each sim
	prob.p.graph = graph_lst[i]
	prob.p.incidence = incidence_matrix(graph_lst[i], oriented=true)
	ilc_ratio = ilc_lst[batch]
	prob.p.hl.ilc_nodes = sample(1:prob.p.N,ilc_ratio,replace = false)
	#prob.p.hl.lambda = lambda_lst[batch]

	hourly_update = network_dynamics.HourlyUpdateEcon()

	ODEProblem(network_dynamics.ACtoymodel!, prob.u0, prob.tspan, prob.p,
		callback=CallbackSet(PeriodicCallback(hourly_update, 3600, initial_affect= false),
							 PeriodicCallback(network_dynamics.DailyUpdate_III, 3600*24, initial_affect= false)))
end

function prob_func_IV(prob, i, repeat, batch_size, graph_lst, num_days, lambda_lst)
	println("sim ", i)
	run = get_run(i, batch_size)
    batch = get_batch(i, batch_size)


	demand_amp = rand(prob.p.N)
	prob.p.periodic_demand = t -> demand_amp .* sin(t*pi/(24*3600))^2
	samples = 24*4
	inter = interpolate([.2 * randn(prob.p.N) for i in 1:(num_days * samples + 1)], BSpline(Linear()))
	prob.p.residual_demand =  t -> inter(1. + t / (24*3600) * samples) # 1. + is needed to avoid trying to access out of range

	prob.p.hl.daily_background_power .= 0.
	# T = 3600
	# for i in 1:24
	# 	prob.p.hl.daily_background_power[i,:] .= sum(prob.p.periodic_demand.((i-1)*T:i*T))/T
	# end

	prob.p.hl.current_background_power .= 0.
	prob.p.hl.mismatch_yesterday .= 0.

	# Random network in each sim
	prob.p.graph = graph_lst[i]
	prob.p.incidence = incidence_matrix(graph_lst[i], oriented=true)
	prob.p.hl.lambda = lambda_lst[batch]

	hourly_update = network_dynamics.HourlyUpdateEcon()

	ODEProblem(network_dynamics.ACtoymodel!, prob.u0, prob.tspan, prob.p,
		callback=CallbackSet(PeriodicCallback(hourly_update, 3600, initial_affect= false),
							 PeriodicCallback(network_dynamics.DailyUpdate_IV, 3600*24, initial_affect= false)))
end

function prob_func_IV_ir(prob, i, repeat, batch_size, graph_lst, num_days, ilc_lst)
	println("sim ", i)
	run = get_run(i, batch_size)
    batch = get_batch(i, batch_size)

	demand_amp = rand(prob.p.N) # .* 100.
	prob.p.periodic_demand = t -> demand_amp .* sin(t*pi/(24*3600))^2
	samples = 24*4
	inter = interpolate([.2 * randn(prob.p.N) for i in 1:(num_days * samples + 1)], BSpline(Linear()))
	prob.p.residual_demand = t -> inter(1. + t / (24*3600) * samples) # 1. + is needed to avoid trying to access out of range

	prob.p.hl.daily_background_power .= 0.
	prob.p.hl.current_background_power .= 0.
	prob.p.hl.mismatch_yesterday .= 0.

	# Random network in each sim
	prob.p.graph = graph_lst[i]
	prob.p.incidence = incidence_matrix(graph_lst[i], oriented=true)
	ilc_ratio = ilc_lst[batch]
	prob.p.hl.ilc_nodes = sample(1:prob.p.N,ilc_ratio,replace = false)
	#prob.p.hl.lambda = lambda_lst[batch]
	kappa_factor =  0.15 / l_hour
	kappa4 =  diagm(0 => zeros(prob.p.N))
	prob.p.hl.ilc_covers = [] # ilc_covers
	for j in 1:ilc_ratio
		kappa4[prob.p.hl.ilc_nodes[j],prob.p.hl.ilc_nodes[j]] = 1
		#kappa[few[1], vc[i]] = 1 # here multiple entries
		#kappa[few[2], vc[i]] = 1
		#kappa[few[3], vc[i]] = 1
		a = 1:prob.p.hl.ilc_nodes[j]-1
		b = prob.p.hl.ilc_nodes[j]+1:prob.p.N
		c = [collect(a); collect(b)]
		few = sample(c,3,replace = false)
		kappa4[prob.p.hl.ilc_nodes[j], few[1]] = 1
		kappa4[prob.p.hl.ilc_nodes[j], few[2]] = 1
		kappa4[prob.p.hl.ilc_nodes[j], few[3]] = 1
		push!(prob.p.hl.ilc_covers, Dict([prob.p.hl.ilc_nodes[j] => few]))
	end
	kappa4 = kappa_factor .* kappa4 .* 0.25
	prob.p.hl.kappa = kappa4

	hourly_update = network_dynamics.HourlyUpdateEcon()

	ODEProblem(network_dynamics.ACtoymodel!, prob.u0, prob.tspan, prob.p,
		callback=CallbackSet(PeriodicCallback(hourly_update, 3600, initial_affect= false),
							 PeriodicCallback(network_dynamics.DailyUpdate_IV, 3600*24, initial_affect= false)))
end


function observer(sol, i, freq_filter, energy_filter, freq_threshold) # what should be extracted from one run
	omega_max = maximum(abs.(sol[freq_filter,:]))
	ex = observables.frequency_exceedance(sol, freq_filter, freq_threshold)
	control_energy = observables.sum_abs_energy_last_days(sol, energy_filter, 2)
	var_omega = var(sol,dims=2)[freq_filter]
	var_ld = observables.var_last_days(sol, freq_filter, 1)
    ((omega_max, ex, control_energy, var_omega, sol.prob.p.graph, sol.prob.p.hl.kappa, sol.prob.p.hl.ilc_nodes, sol.prob.p.hl.ilc_covers, var_ld), false)
end

function observer_basic_types(sol, i, freq_filter, energy_filter, freq_threshold) # what should be extracted from one run
	omega_max = maximum(abs.(sol[freq_filter,:]))
	ex = observables.frequency_exceedance(sol, freq_filter, freq_threshold, sol.prob.tspan[1])
	control_energy = observables.sum_abs_energy_last_days(sol, energy_filter, sol.prob.tspan[2]/(24*3600))
	var_omega = var(sol,dims=2)[freq_filter]
	var_ld = observables.var_last_days(sol, freq_filter, sol.prob.tspan[2]/(24*3600))
	((omega_max, ex, control_energy, var_omega, Array(adjacency_matrix(sol.prob.p.graph)), sol.prob.p.hl.kappa, sol.prob.p.hl.ilc_nodes, sol.prob.p.hl.ilc_covers, var_ld), false)
end


function observer_basic_types_lambda(sol, i, freq_filter, energy_filter, energy_abs_filter, freq_threshold, obs_days) # what should be extracted from one run
	num_days = Int(sol.prob.tspan[2]/(24*3600))
	omega_max = maximum(abs.(sol(sol.prob.tspan[2]-obs_days*3600*24:sol.prob.tspan[2])[freq_filter,:]))
	ex = observables.frequency_exceedance(sol, freq_filter, freq_threshold, obs_days)
	control_energy = observables.sum_abs_energy_last_days(sol, energy_filter, obs_days)
	var_omega = var(sol,dims=2)[freq_filter]
	var_ld = observables.var_last_days(sol, freq_filter, obs_days)
	if isdefined(sol.prob.kwargs.data,:callback)
		cost = observables.energy_cost(sol, energy_filter, energy_abs_filter, num_days, obs_days)
	else
		energy_sum = zeros(sol.prob.p.N)
		for j=1:sol.prob.p.N
			energy_sum[j] = sol(num_days*24*3600.)[energy_abs_filter[j]] .- sol(num_days*24*3600. - obs_days*24*3600+1)[energy_abs_filter[j]]
		end
		cost = sol.prob.p.hl.lambda .* sum(abs.(energy_sum))./3600.
	end
	((omega_max, ex, control_energy, var_omega, Array(adjacency_matrix(sol.prob.p.graph)), sol.prob.p.hl.kappa, sol.prob.p.hl.ilc_nodes, sol.prob.p.hl.ilc_covers, var_ld, cost), false)
end


function observer_basic_types_lambda_II(sol, i, freq_filter, energy_filter, energy_abs_filter, freq_threshold, obs_days) # what should be extracted from one run
	num_days = Int(sol.prob.tspan[2]/(24*3600))
	omega_max = maximum(abs.(sol(sol.prob.tspan[2]-obs_days*3600*24:sol.prob.tspan[2])[freq_filter,:]))
	ex = observables.frequency_exceedance(sol, freq_filter, freq_threshold, obs_days)
	control_energy = observables.sum_abs_energy_last_days(sol, energy_filter, obs_days)
	var_omega = var(sol,dims=2)[freq_filter]
	var_ld = observables.var_last_days(sol, freq_filter, obs_days)
	if isdefined(sol.prob.kwargs.data,:callback)
		cost = observables.energy_cost_II(sol, energy_filter, energy_abs_filter, num_days, obs_days)
	else
		energy_sum = zeros(sol.prob.p.N)
		for j=1:sol.prob.p.N
			energy_sum[j] = sol(num_days*24*3600.)[energy_abs_filter[j]] .- sol(num_days*24*3600. - obs_days*24*3600+1)[energy_abs_filter[j]]
		end
		cost = sol.prob.p.hl.lambda .* sum(abs.(energy_sum))./3600.
	end
	((omega_max, ex, control_energy, var_omega, Array(adjacency_matrix(sol.prob.p.graph)), sol.prob.p.hl.kappa, sol.prob.p.hl.ilc_nodes, sol.prob.p.hl.ilc_covers, var_ld, cost), false)
end

function observer_basic_types_lambda_III(sol, i, freq_filter, energy_filter, energy_abs_filter, freq_threshold, obs_days) # what should be extracted from one run
	num_days = Int(sol.prob.tspan[2]/(24*3600))
	omega_max = maximum(abs.(sol(sol.prob.tspan[2]-obs_days*3600*24:sol.prob.tspan[2])[freq_filter,:]))
	ex = observables.frequency_exceedance(sol, freq_filter, freq_threshold, obs_days)
	control_energy = observables.sum_abs_energy_last_days(sol, energy_filter, obs_days)
	var_omega = var(sol,dims=2)[freq_filter]
	var_ld = observables.var_last_days(sol, freq_filter, obs_days)
	if isdefined(sol.prob.kwargs.data,:callback)
		cost = observables.energy_cost_III(sol, energy_filter, energy_abs_filter, num_days, obs_days)
	else
		energy_sum = zeros(sol.prob.p.N)
		for j=1:sol.prob.p.N
			energy_sum[j] = sol(num_days*24*3600.)[energy_abs_filter[j]] .- sol(num_days*24*3600. - obs_days*24*3600+1)[energy_abs_filter[j]]
		end
		cost = sol.prob.p.hl.lambda .* sum(abs.(energy_sum))./3600.
	end
	((omega_max, ex, control_energy, var_omega, Array(adjacency_matrix(sol.prob.p.graph)), sol.prob.p.hl.kappa, sol.prob.p.hl.ilc_nodes, sol.prob.p.hl.ilc_covers, var_ld, cost), false)
end

function observer_basic_types_lambda_IV(sol, i, freq_filter, energy_filter, energy_abs_filter, freq_threshold, obs_days) # what should be extracted from one run
	num_days = Int(sol.prob.tspan[2]/(24*3600))
	omega_max = maximum(abs.(sol(sol.prob.tspan[2]-obs_days*3600*24:sol.prob.tspan[2])[freq_filter,:]))
	ex = observables.frequency_exceedance(sol, freq_filter, freq_threshold, obs_days)
	control_energy = observables.sum_abs_energy_last_days(sol, energy_filter, obs_days)
	var_omega = var(sol,dims=2)[freq_filter]
	var_ld = observables.var_last_days(sol, freq_filter, obs_days)
	if isdefined(sol.prob.kwargs.data,:callback)
		cost = observables.energy_cost_IV(sol, energy_filter, energy_abs_filter, num_days, obs_days)
	else
		energy_sum = zeros(sol.prob.p.N)
		for j=1:sol.prob.p.N
			energy_sum[j] = sol(num_days*24*3600.)[energy_abs_filter[j]] .- sol(num_days*24*3600. - obs_days*24*3600+1)[energy_abs_filter[j]]
		end
		cost = sol.prob.p.hl.lambda .* sum(abs.(energy_sum))./3600.
	end
	((omega_max, ex, control_energy, var_omega, Array(adjacency_matrix(sol.prob.p.graph)), sol.prob.p.hl.kappa, sol.prob.p.hl.ilc_nodes, sol.prob.p.hl.ilc_covers, var_ld, cost), false)
end

function observer_ll(sol, i, freq_filter, energy_filter, freq_threshold) # what should be extracted from one run
	omega_max = maximum(abs.(sol[freq_filter,:]))
	ex = observables.frequency_exceedance(sol, freq_filter, freq_threshold)
	control_energy = observables.sum_abs_energy_last_days(sol, energy_filter, 1)
    ((omega_max, ex, control_energy), false)#, var_omega_nodemean, var_ld_nodemean), false)
end



function reduction(u, data, I, batch_size) # what should be extracted from one batch
    # u is the solution of previous batches
    # data is the solution of the current batch
    # we obtain:
	omega_max_abs = [dat[1] for dat in data] # max frequency
    	ex = [dat[2] for dat in data] # ex
	energy = [dat[3] for dat in data] # working as array???

	omega_max_max = maximum(omega_max_abs)
	ex_mean = sum(ex)/batch_size
	energy_mean = sum(energy)/length(energy) # summing over all nodes and runs in one batches now
	ex_std = std(ex)
	energy_std = std(energy)
	new_output = [omega_max_max, ex_mean, energy_mean, ex_std, energy_std] #var_omega_max var_omega_ld_mean]
	append!(u, [new_output, ]), false # This way of append! ensures that we get an Array of Arrays
end

##############################################################################
# other



function delete_edges_in_graph!(prob, graph_init, fault_size)
	new_graph = copy(graph_init)
	edges_to_delete = sample(edges(new_graph) |> collect, fault_size, replace=false)
	if length(edges_to_delete) == 1
		rem_edge!(new_graph, edges_to_delete[1])
	else
		for e in edges_to_delete
			rem_edge!(new_graph, e)
		end
	end
	prob.p.graph = new_graph
	prob.p.incidence = incidence_matrix(new_graph, oriented=true)
end



################################################
# plotting help

function movingmean(t,int,a,b)
    idx = findall(x -> t + int > x > t - int, a)
    mean(b[idx])
end

end

#@doc """
#    solver_test(prob::DiffEqBase.DEProblem)
#
#Compare the performance of different solvers for a given problem (`prob`) in DifferentialEquations.jl.
#"""
#function solver_test(prob::DiffEqBase.DEProblem)
#	# The DiffEq documentation discusses a bunch of solver options here:
#	# https://docs.juliadiffeq.org/latest/solvers/ode_solve.html
#
#	# auto-switching
#	println("Test auto-switching methods:")
#	autoswitching = [plot(), plot(), plot()]
#	for (i, solver) in enumerate([AutoTsit5, AutoVern7, AutoVern9])
#		composite = solver(Rodas4())
#		time = minimum([@elapsed solve(prob, composite) for trial in 1:3])
#		name = "$(nameof(solver))(Rodas4)"
#		println("$name best time: $time")
#		sol = solve(prob, composite)
#		plot!(autoswitching[i], sol, vars=prob.p.N+1, tspan=(0., 10.), label=name, alpha=0.25)
#
#		composite = solver(KenCarp4())
#		time = minimum([@elapsed solve(prob, composite) for trial in 1:3])
#		name = "$(nameof(solver))(KenCarp4)"
#		println("$name best time: $time")
#		sol = solve(prob, composite)
#		plot!(autoswitching[i], sol, vars=prob.p.N+1, tspan=(0., 10.), label=name, alpha=0.4)
#
#		composite = solver(Rodas5())
#		time = minimum([@elapsed solve(prob, composite) for trial in 1:3])
#		name = "$(nameof(solver))(Rodas5)"
#		println("$name best time: $time")
#		sol = solve(prob, composite)
#		plot!(autoswitching[i], sol, vars=prob.p.N+1, tspan=(0., 10.), label=name, alpha=0.4)
#	end
#
#	# stiff solvers
#	println("Test recommended stiff methods:")
#	stiff = [plot() for i in 1:8]
#	for (i, solver) in enumerate([Rosenbrock23, TRBDF2, ABDF2, Rodas4P, Rodas5, Kvaerno5, KenCarp4, CVODE_BDF])
#		time = minimum([@elapsed solve(prob, solver()) for trial in 1:3])
#		name = "$(nameof(solver))"
#		println("$name best time: $time")
#		sol = solve(prob, solver())
#		plot!(stiff[i], sol, vars=prob.p.N+1, tspan=(0., 10.), title=name, alpha=0.4)
#	end
#
#	savefig(plot(autoswitching..., legendfontsize=5), "autoswitching.pdf")
#	savefig(plot(stiff..., legend=:none), "recommended_stiff.pdf")
#
#	return nothing
#end
