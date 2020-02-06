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
	freq_threshold = 0.0005
	phase_filter = 1:N
	freq_filter = N+1:2N
	control_filter = 2N+1:3N
	energy_filter = 3N+1:4N
	energy_abs_filter = 4N+1:5N
end

# @load "$dir/solutions_only_MC_no_time_series/20190913/exp0_sol_new.jld2" res_tl0 monte_prob0
# @load "$dir/solutions_only_MC_no_time_series/20190913/expI_sol_new.jld2" res_tlI monte_probI
# @load "$dir/solutions_only_MC_no_time_series/20190913/expII_sol_new.jld2" res_tlII monte_probII
# @load "$dir/solutions_only_MC_no_time_series/20190913/expIII_sol_new.jld2" res_tlIII monte_probIII
# @load "$dir/solutions_only_MC_no_time_series/20190913/expIV_sol_new.jld2" res_tlIV monte_probIV


using Dates
date = Dates.today()

@load "$dir/solutions/$(date)/exp0_sol_new.jld2" res_tl0 monte_prob0
@load "$dir/solutions/$(date)/expI_sol_new.jld2" res_tlI monte_probI
@load "$dir/solutions/$(date)/expII_sol_new.jld2" res_tlII monte_probII
@load "$dir/solutions/$(date)/expIII_sol_new.jld2" res_tlIII monte_probIII
@load "$dir/solutions/$(date)/expIV_sol_new.jld2" res_tlIV monte_probIV

date = Dates.today() + Dates.Day(1)

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

###########################################################
#################### OLD ##################################
##########################################################

# omega_max0 = maximum(abs.([p[1] for p in res_tl0.u]))
# omega_max1 = maximum(abs.([p[1] for p in res_tlI.u]))
# omega_max2 = maximum(abs.([p[1] for p in res_tlII.u]))
# omega_max3 = maximum(abs.([p[1] for p in res_tlIII.u]))
# omega_max4 = maximum(abs.([p[1] for p in res_tlIV.u]))
#
# ex0 = maximum([p[2] for p in res_tl0.u])
# ex1 = maximum([p[2] for p in res_tlI.u])
# ex2 = maximum([p[2] for p in res_tlII.u])
# ex3 = maximum([p[2] for p in res_tlIII.u])
# ex4 = maximum([p[2] for p in res_tlIV.u])
#
#
# seconds = 3600*24*num_days
# plot([0 1 2 3 4], [ex0*100/seconds ex1*100/seconds ex3*100/seconds ex2*100/seconds ex4*100/seconds], seriestype=:bar,  xticks = (0:4, ["no ilc", "local ilc all nodes", "local ilc at vc", "ilc at vc \n+ average neighbors", "ilc at 50% random \n + 3 averag random"]), legend = :none)
# ylabel!("freq exceedance [% of total time]")
# title!("Frequency exceedance")
# savefig("$dir/plots/20190913/frequency_exceedance_percent.png")
#
# plot([0 1 2 3], [ex0 ex1 ex3 ex2], seriestype=:bar,  xticks = (0:3, ["no ilc", "local ilc all nodes", "local ilc at vc", "ilc at vc \n+ average neighbors"]), legend = :none)
# ylabel!("freq exceedance 0-III")
# title!("Frequency exceedance 0-III")
# savefig("$dir/plots/20190913/frequency_exceedance0-III.png")
#
#
# plot([0 1 2 3 4], [omega_max0  omega_max1 omega_max3 omega_max2 omega_max4], seriestype=:bar,  xticks = (0:4, ["no ilc", "local ilc all nodes", "local ilc at vc", "ilc at vc \n+ average neighbors", "ilc at 50% random \n + 3 averag random"]), legend = :none)
# ylabel!("max freq deviation")
# title!("Max frequency deviation")
# savefig("$dir/plots/20190913/max_freq_dev.png")
#
#
# y_f = [ex0 ex1 ex3 ex2 ex4; omega_max0  omega_max1 omega_max3 omega_max2 omega_max4]'
#
# plot(collect(0:4), y_f, seriestype=:scatter,  xticks = (0:4, ["no ilc", "local ilc all nodes", "local ilc at vc", "ilc at vc \n+ average neighbors", "ilc at 50% random \n + 3 averag random"]), label = ["freq dev", "max freq dev"])
# ylabel!("freq exceedance and maximum frequency deviation")
# title!("Frequency exceedance and maximum freq deviation")
# savefig("$dir/plots/20190913/summary_frequency.png")
#
#
#
#
# pnceld0 = mean([p[3] for p in res_tl0.u]) #mean over experiments
# pnceld1 = mean([p[3] for p in res_tlI.u])
# pnceld2 = mean([p[3] for p in res_tlII.u])
# pnceld3 = mean([p[3] for p in res_tlIII.u])
# pnceld4 = mean([p[3] for p in res_tlIV.u])
#
#
#
# plot(1:N, pnceld0,seriestype=:scatter)
# plot!(1:N, pnceld1,seriestype=:scatter)
# plot!(1:N, pnceld2,seriestype=:scatter)
# plot!(1:N, pnceld3,seriestype=:scatter)
# plot!(1:N, pnceld4,seriestype=:scatter)
#
# savefig("$dir/plots/20190913/sim0-IV_pnceld.png")
#
#
# plot(vc2, pnceld0,seriestype=:scatter)
# plot!(vc2, pnceld1,seriestype=:scatter)
# plot!(vc2, pnceld2,seriestype=:scatter)
# plot!(vc3, pnceld3,seriestype=:scatter)
# plot!(vc4, pnceld4,seriestype=:scatter)
#
# savefig("$dir/plots/20190913/sim0-IV_pnceld_vc.png")

#
#
# plot(nonvc2, pnceld0,seriestype=:scatter)
# plot!(nonvc2, pnceld1,seriestype=:scatter)
# plot!(nonvc2, pnceld2,seriestype=:scatter)
# plot!(nonvc3, pnceld3,seriestype=:scatter)
# plot!(nonvc4, pnceld4,seriestype=:scatter)
#
# savefig("$dir/plots/20190913/sim0-IV_pnceld_nonvc.png")
#
#
#
# pnceld_ilc_0 = mean(mean.([[p[3] for p in res_tl0.u][i][vc2[i]] for i in 1:num_monte]))
# pnceld_ilc_1 = mean(mean.([[p[3] for p in res_tlI.u][i][vc2[i]] for i in 1:num_monte]))
# pnceld_ilc_2 = mean(mean.([[p[3] for p in res_tlII.u][i][vc2[i]] for i in 1:num_monte]))
# pnceld_ilc_3 = mean(mean.([[p[3] for p in res_tlIII.u][i][vc3[i]] for i in 1:num_monte]))
# pnceld_ilc_4 = mean(mean.([[p[3] for p in res_tlIV.u][i][vc2[i]] for i in 1:num_monte]))
#
# pnceld_nonilc_0 = mean(mean.([[p[3] for p in res_tl0.u][i][nonvc2[i]] for i in 1:num_monte]))
# pnceld_nonilc_1 = mean(mean.([[p[3] for p in res_tlI.u][i][nonvc2[i]] for i in 1:num_monte]))
# pnceld_nonilc_2 = mean(mean.([[p[3] for p in res_tlII.u][i][nonvc2[i]] for i in 1:num_monte]))
# pnceld_nonilc_3 = mean(mean.([[p[3] for p in res_tlIII.u][i][nonvc3[i]] for i in 1:num_monte]))
# pnceld_nonilc_4 = mean(mean.([[p[3] for p in res_tlIV.u][i][nonvc2[i]] for i in 1:num_monte]))
#
# pnceld_nonilc_0a = mean(mean.([[p[3] for p in res_tl0.u][i][1:N] for i in 1:num_monte]))
# pnceld_ilc_1a = mean(mean.([[p[3] for p in res_tlI.u][i][1:N] for i in 1:num_monte]))
# pnceld_all_2a = mean(mean.([[p[3] for p in res_tlII.u][i][1:N] for i in 1:num_monte]))
# pnceld_all_3a = mean(mean.([[p[3] for p in res_tlIII.u][i][1:N] for i in 1:num_monte]))
# pnceld_all_4a = mean(mean.([[p[3] for p in res_tlIV.u][i][1:N] for i in 1:num_monte]))
#
#
# pnceld_ilc_4a = mean(mean.([[p[3] for p in res_tlIV.u][i][vc4[i]] for i in 1:num_monte]))
# pnceld_nonilc_4b = mean(mean.([[p[3] for p in res_tlIV.u][i][nonvc4[i]] for i in 1:num_monte]))
#
#
#
#
#
# plot([0 1 2 3 4], [pnceld_nonilc_0a pnceld_ilc_1a pnceld_all_3a pnceld_all_2a pnceld_all_4a], seriestype=:bar,  xticks = (0:4, ["no ilc", "local ilc all nodes", "local ilc at vc", "ilc at vc \n+ average neighbors", "ilc at 50% random \n + 3 averag random"]), legend = :none)
# ylabel!("control energy of the last day averaged over all nodes")
# title!("All nodes")
# savefig("$dir/plots/20190913/control_energy_all_nodes_average.png")
#
# plot([0 1 2 3 4], [pnceld_nonilc_0a pnceld_ilc_1a pnceld_ilc_3 pnceld_ilc_2 pnceld_ilc_4a], seriestype=:bar,  xticks = (0:4, ["no ilc \n(benchmark all nodes)", "local ilc all nodes", "local ilc at vc", "ilc at vc \n+ average neighbors", "ilc at 50% random \n + 3 average random"]), legend = :none)
# ylabel!("control energy of the last day averaged over all ILC nodes")
# title!("ILC nodes")
# savefig("$dir/plots/20190913/control_energy_ilc_nodes_average.png")
#
#
# plot([0 1 2 3 4], [pnceld_nonilc_0a pnceld_ilc_1a pnceld_nonilc_3 pnceld_nonilc_2 pnceld_nonilc_4b], seriestype=:bar,  xticks = (0:4, ["no ilc", "local ilc all nodes\n (benchmark all nodes)", "local ilc at vc", "ilc at vc \n+ average neighbors", "ilc at 50% random \n + 3 average random"]), legend = :none)
# ylabel!("control energy last day averaged over all non-ILC nodes")
# title!("Non-ILC nodes")
# savefig("$dir/plots/20190913/control_energy_nonilc_nodes_average.png")
#
#
# y = [pnceld_nonilc_0a pnceld_ilc_1a pnceld_ilc_3 pnceld_ilc_2 pnceld_ilc_4a; pnceld_nonilc_0a pnceld_ilc_1a pnceld_nonilc_3 pnceld_nonilc_2 pnceld_nonilc_4b; pnceld_nonilc_0a pnceld_ilc_1a pnceld_all_3a pnceld_all_2a pnceld_all_4a]'
#
# plot(collect(0:4), y, seriestype=:scatter,  xticks = (0:4, ["no ilc", "local ilc all nodes", "local ilc at vc", "ilc at vc \n+ average neighbors", "ilc at 50% random \n + 3 average random"]), label = ["ILC nodes", "non-ILC nodes", "all nodes"])
# #plot!([0 1 2 3 4], seriestype=:groupedbar)#,  xticks = (0:4, ["no ilc", "local ilc all nodes\n (benchmark all nodes)", "local ilc at vc", "ilc at vc \n+ average neighbors", "ilc at 50% random \n + 3 average random"]), legend = :none)
# #plot!([0 1 2 3 4],, seriestype=:groupedbar)#,  xticks = (0:4, ["no ilc", "local ilc all nodes", "local ilc at vc", "ilc at vc \n+ average neighbors", "ilc at 50% random \n + 3 averag random"]), legend = :none)
# title!("Averaged control energy")
# ylabel!("Control energy last day averaged over certain nodes")
# savefig("$dir/plots/20190913/control_energy_average.png")
#
# yI_IV = [pnceld_ilc_1a pnceld_ilc_3 pnceld_ilc_2 pnceld_ilc_4a; pnceld_ilc_1a pnceld_nonilc_3 pnceld_nonilc_2 pnceld_nonilc_4b; pnceld_ilc_1a pnceld_all_3a pnceld_all_2a pnceld_all_4a]'
#
# plot(collect(1:4), yI_IV, seriestype=:scatter,  xticks = (1:4, ["local ilc all nodes", "local ilc at vc", "ilc at vc \n+ average neighbors", "ilc at 50% random \n + 3 average random"]), label = ["ILC nodes", "non-ILC nodes", "all nodes"])
# #plot!([0 1 2 3 4], seriestype=:groupedbar)#,  xticks = (0:4, ["no ilc", "local ilc all nodes\n (benchmark all nodes)", "local ilc at vc", "ilc at vc \n+ average neighbors", "ilc at 50% random \n + 3 average random"]), legend = :none)
# #plot!([0 1 2 3 4],, seriestype=:groupedbar)#,  xticks = (0:4, ["no ilc", "local ilc all nodes", "local ilc at vc", "ilc at vc \n+ average neighbors", "ilc at 50% random \n + 3 averag random"]), legend = :none)
# title!("Averaged control energy")
# ylabel!("Control energy last day averaged over certain nodes")
# savefig("$dir/plots/20190913/control_energy_averageI-IV.png")
