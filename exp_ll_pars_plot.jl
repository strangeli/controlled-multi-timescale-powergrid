using JLD2, FileIO
using Distributed

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

test = false

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
	freq_threshold = 0.001
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

using Plots
@load "$dir/solutions/20200209_exp_ll_pars.jld2" res monte_prob


#factor = maximum(abs, res.u, dims=1)
plot(collect(kP_lst), [p[1] for p in res.u], label=["omega_max"])
xlabel!("k_P")
ylabel!("omega max")
savefig("$dir/plots/exp_ll_freq_dev.pdf")
plot(collect(kP_lst), [p[2] for p in res.u], label=["exceedance"])
xlabel!("k_P")
ylabel!("exceedance")
savefig("$dir/plots/exp_ll_exceedance.pdf")
plot(collect(kP_lst),mean.([p[3] for p in res.u]), label=["control energy"])
xlabel!("k_P")
ylabel!("averaged control energy")
savefig("$dir/plots/exp_ll_control_energy.pdf")
plot(collect(kP_lst), mean.([p[4] for p in res.u]), label=["freq_var"])
xlabel!("k_P")
ylabel!("freq var")
savefig("$dir/plots/exp_ll_freq_var.pdf")

using LaTeXStrings

p1 = contour(kP_lst_s,kI_lst_s, reshape([p[1] for p in res.u],length(kI_lst_s), length(kP_lst_s)), fill = true,title_location=:right, ylims=(0,0.2), clims=(0,0.005), thickness_scaling=1.1, dpi=150,#colorbar_title=L"\omega_{max}",#,ylims=(-100,2100),
		   ytickfontsize=12,
		   xtickfontsize=12,
		   guidefontsize=12, top_margin=5Plots.mm, left_margin=5Plots.mm, bottom_margin=5Plots.mm, right_margin=8Plots.mm) # ,clims=(0,0.5)xlabel!("k_P")
		  # title!(L"\omega_{max}")
		  xlabel!(L"k_P [s/rad]")
		  ylabel!(L"k_I [rad/s]")
title!("Maximum frequency deviation [rad/s]", titlefontsize=12)
savefig("$dir/plots/contour_kpki_omega_max-02.pdf")

p2 = contour(kP_lst_s,kI_lst_s, reshape([log10(p[2]*100/(24*3600) + 1e-5) for p in res.u],length(kI_lst_s), length(kP_lst_s)), fill = true, thickness_scaling=1.1, ylims=(0,0.2), margin=5Plots.mm, ytickfontsize=12,
		   xtickfontsize=12,
		   guidefontsize=12)
xlabel!(L"k_P [s/rad]")
ylabel!(L"k_I [rad/s]")
title!("log10(exceedance [%])", titlefontsize=12)
savefig("$dir/plots/contour_kpki_ex_log.pdf")

p4 = contour(kP_lst_s,kI_lst_s, reshape([log10(mean(p[4]) + 1E-11) for p in res.u],length(kI_lst_s), length(kP_lst_s)), fill = true, thickness_scaling=1.1, ylims=(0,0.2),margin=5Plots.mm, ytickfontsize=12,
		   xtickfontsize=12,
		   guidefontsize=12)#,clims =(-10,5)) #,xlims = (400,600),ylims=(0,0.1)
		   xlabel!(L"k_P [s/rad]")
		   ylabel!(L"k_I [rad/s]")
title!("log10(frequency variance [rad/s])", titlefontsize=12)
savefig("$dir/plots/contour_kpki_freq_var_log.pdf")

#------------------------------

p3 = contour(kP_lst_s,kI_lst_s, reshape([mean(p[3]) for p in res.u],length(kI_lst_s), length(kP_lst_s)), fill = true, colorbar_title="control energy [s]", thickness_scaling=1.1,ylims=(0,0.1))#,clims=(14000,15000))
xlabel!("k_P")
ylabel!("k_I")
savefig("$dir/plots/contour_kpki_control_energy.pdf")




p2 = contour(kP_lst_s,kI_lst_s, reshape([p[2] for p in res.u],length(kI_lst_s), length(kP_lst_s)), fill = true, colorbar_title="exceedance [s]", thickness_scaling=1.1, ylims=(0,0.2), clims=(0,50))
xlabel!(L"k_P")
ylabel!(L"k_I")
savefig("$dir/plots/contour_kpki_ex.pdf")


###################################################################################

y1 = reshape([p[1] for p in res.u], length(kI_lst_s), length(kP_lst_s))
y2 = reshape([p[2] for p in res.u], length(kI_lst_s), length(kP_lst_s))
y3 = reshape([p[3] for p in res.u], length(kI_lst_s), length(kP_lst_s))
y4 = reshape([p[4] for p in res.u], length(kI_lst_s), length(kP_lst_s))


plot(collect(kI_lst_s), y1, label="omega_max")
xlabel!("k_I")
ylabel!("omega max")
savefig("$dir/plots/exp_ll_freq_dev_kI.pdf")
plot(collect(kI_lst_s),y2,label="exceedance")
xlabel!("k_I")
ylabel!("exceedance")
savefig("$dir/plots/exp_ll_exceedance_kI.pdf")
plot(collect(kI_lst_s),mean.(y3),seriestype=:scatter, label="control energy")
xlabel!("k_I")
ylabel!("averaged control energy")
savefig("$dir/plots/exp_ll_control_energy_kI.pdf")
plot(collect(kI_lst_s), mean.(y4),seriestype=:scatter, label="ex_std")
xlabel!("k_I")
ylabel!("freq var")
savefig("$dir/plots/exp_ll_freq_var_kI.pdf")
