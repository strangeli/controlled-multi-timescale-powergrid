for m = 1:length(lambda_lst) # lambda = 0.2:0.1:1
        println("m ",m)
        lam = lambda_lst[m]
	# [(m-1)*batch_size+1:m*batch_size]
	vc2 = [p[7] for p in res_tlII.u][(m-1)*batch_size+1:m*batch_size]
	vc3 = [p[7] for p in res_tlIII.u][(m-1)*batch_size+1:m*batch_size]
	vc4 = [p[7] for p in res_tlIV.u][(m-1)*batch_size+1:m*batch_size]

	nonvc2 = repeat([Int.(zeros(N))], outer = batch_size)
	for i in 1:batch_size
		noilc = []
		for k in 1:N
				if k ∉ [p[1] for p in vc2[i]]
					append!(noilc, k)
				end
		end
		nonvc2[i] = noilc
	end

	nonvc3 = repeat([Int.(zeros(N))], outer = batch_size)
	for i in 1:batch_size
		noilc = []
		for k in 1:N
				if k ∉ [p[1] for p in vc3[i]]
					append!(noilc, k)
				end
		end
		nonvc3[i] = noilc
	end

	nonvc4 = repeat([Int.(zeros(N))], outer = batch_size)
	for i in 1:batch_size
		noilc = []
		for k in 1:N
				if k ∉ [p[1] for p in vc4[i]]
					append!(noilc, k)
				end
		end
		nonvc4[i] = noilc
	end

	graph_lst = Graph.([p[5] for p in res_tlIV.u][(m-1)*batch_size+1:m*batch_size])


	##############################################
	######### VIOLIN PLOTS #######################
	##############################################

	# [(m-1)*batch_size+1:m*batch_size]

	pnceld_ilc_2 = mean.([[p[3] for p in res_tlII.u][i][vc2[index_calc(i,batch_size)]] for i in (m-1)*batch_size+1:m*batch_size])
	pnceld_ilc_3 = mean.([[p[3] for p in res_tlIII.u][i][vc3[index_calc(i,batch_size)]] for i in (m-1)*batch_size+1:m*batch_size])

	pnceld_nonilc_2 = mean.([[p[3] for p in res_tlII.u][i][nonvc2[index_calc(i,batch_size)]] for i in (m-1)*batch_size+1:m*batch_size])
	pnceld_nonilc_3 = mean.([[p[3] for p in res_tlIII.u][i][nonvc3[index_calc(i,batch_size)]] for i in (m-1)*batch_size+1:m*batch_size])

	pnceld_nonilc_0 = mean.([[p[3] for p in res_tl0.u][i][1:N] for i in (m-1)*batch_size+1:m*batch_size])
	pnceld_ilc_1 = mean.([[p[3] for p in res_tlI.u][i][1:N] for i in (m-1)*batch_size+1:m*batch_size])

	pnceld_ilc_4 = mean.([[p[3] for p in res_tlIV.u][i][vc4[index_calc(i,batch_size)]] for i in (m-1)*batch_size+1:m*batch_size])
	pnceld_nonilc_4 = mean.([[p[3] for p in res_tlIV.u][i][nonvc4[index_calc(i,batch_size)]] for i in (m-1)*batch_size+1:m*batch_size])

	pnceld_all_2 = mean.([[p[3] for p in res_tlII.u][i][1:N] for i in (m-1)*batch_size+1:m*batch_size])
	pnceld_all_3 = mean.([[p[3] for p in res_tlIII.u][i][1:N] for i in (m-1)*batch_size+1:m*batch_size])
	pnceld_all_4 = mean.([[p[3] for p in res_tlIV.u][i][1:N] for i in (m-1)*batch_size+1:m*batch_size])

	################################################################

	# pnceld_ilc_2abs = mean.([[p[3] for p in res_tlIIabs.u][i][vc2[i]] for i in 1:batch_size])
	# pnceld_ilc_3abs = mean.([[p[3] for p in res_tlIIIabs.u][i][vc3[i]] for i in 1:batch_size])
	#
	# pnceld_nonilc_2abs = mean.([[p[3] for p in res_tlIIabs.u][i][nonvc2[i]] for i in 1:batch_size])
	# pnceld_nonilc_3abs = mean.([[p[3] for p in res_tlIIIabs.u][i][nonvc3[i]] for i in 1:batch_size])
	#
	# pnceld_nonilc_0abs = mean.([[p[3] for p in res_tl0abs.u][i][1:N] for i in 1:batch_size])
	# pnceld_ilc_1abs = mean.([[p[3] for p in res_tlIabs.u][i][1:N] for i in 1:batch_size])
	#
	# pnceld_ilc_4abs = mean.([[p[3] for p in res_tlIVabs.u][i][vc4[i]] for i in 1:batch_size])
	# pnceld_nonilc_4abs = mean.([[p[3] for p in res_tlIVabs.u][i][nonvc4[i]] for i in 1:batch_size])
	#
	# pnceld_all_2abs = mean.([[p[3] for p in res_tlIIabs.u][i][1:N] for i in 1:batch_size])
	# pnceld_all_3abs = mean.([[p[3] for p in res_tlIIIabs.u][i][1:N] for i in 1:batch_size])
	# pnceld_all_4abs = mean.([[p[3] for p in res_tlIVabs.u][i][1:N] for i in 1:batch_size])
	#


	# import RDatasets
	# https://github.com/JuliaPlots/StatsPlots.jl
	# singers = RDatasets.dataset("lattice","singer")
	#using StatsPlots

	#y = [pnceld_nonilc_0a pnceld_ilc_1a pnceld_ilc_3 pnceld_ilc_2 pnceld_ilc_4a; pnceld_nonilc_0a pnceld_ilc_1a pnceld_nonilc_3 pnceld_nonilc_2 pnceld_nonilc_4b; pnceld_nonilc_0a pnceld_ilc_1a pnceld_all_3a pnceld_all_2a pnceld_all_4a]'

	#ExpResults_ilc_nodes_0 = DataFrame(Experiment = "no ilc", Result = pnceld_ilc_0)
	ExpResults_ilc_nodes_I = DataFrame(Experiment = "I", Result = pnceld_ilc_1)
	ExpResults_ilc_nodes_III = DataFrame(Experiment = "III", Result = pnceld_ilc_3)
	ExpResults_ilc_nodes_II = DataFrame(Experiment = "II", Result = pnceld_ilc_2)
	ExpResults_ilc_nodes_IV = DataFrame(Experiment = "IV", Result = pnceld_ilc_4)


	begin
	    plt_ex0 = plot(legend=false, dpi=150, # yscale=:log10,
				   ytickfontsize=18,
				   xtickfontsize=14,
				   guidefontsize=20,
				   legendfontsize=10)
		@df ExpResults_ilc_nodes_I boxplot!(plt_ex0, :Experiment,:Result) #,fill=(0,0.5,:orange)) # ,marker=(0.2,:blue,stroke(0))
		@df ExpResults_ilc_nodes_I boxplot!(plt_ex0, :Experiment,:Result) #,fill=(0,0.5,:orange)) # ,marker=(0.2,:blue,stroke(0))
	    @df ExpResults_ilc_nodes_II boxplot!(plt_ex0, :Experiment,:Result) # ,marker=(0.2,:blue,stroke(0))
	    @df ExpResults_ilc_nodes_III boxplot!(plt_ex0, :Experiment,:Result) # ,marker=(0.2,:blue,stroke(0))
	    ylabel!("Control energy [Ws]")

	    #plt_ex1 = plot(legend=false, yscale=:log10)
	    @df ExpResults_ilc_nodes_IV boxplot!(plt_ex0, :Experiment,:Result,fill=(0,0.5,:orange), whisker_width=1) # ,marker=(0.2,:blue,stroke(0))
	    #ylabel!("Exceedance [s]")

	    #l = @layout [a{0.7w} b{0.2w}]
	    plt_ex = plot(plt_ex0)#,plt_ex1, layout = l)
	end
	savefig(plt_ex, "$dir/plots/$(date)/control_ILC_nodes_boxplot_lambda_$lam.pdf")



	#
	#
	#
	# # for el in pnceld_ilc_3
	# #     push!(ExpResults_ilc_nodes_I_III, ("III: local ilc at vc", el))
	# # end
	#
	# plt_ilc_nodes_0 = plot(legend=false)
	#
	#
	# plt_ilc_nodes_I = plot(legend=false,ylims=(14000,20000))
	# @df ExpResults_ilc_nodes_I violin!(plt_ilc_nodes_I, :Experiment,:Result, fill=(0,0.5,:orange)) # ,marker=(0.2,:blue,stroke(0))
	# #@df ExpResults_ilc_nodes_I_III dotplot!(plt_ilc_nodes_I_III, :Experiment,:Result) # ,marker=(0.2,:blue,stroke(0))
	# ylabel!("Control energy [Ws]")
	# title!("ILC nodes")
	# savefig(plt_ilc_nodes_I, "$dir/plots/$(date)/violin_ilc_nodes_I.pdf")
	#
	# plt_ilc_nodes_III = plot(legend=false,ylims=(14000,20000))
	# @df ExpResults_ilc_nodes_III violin!(plt_ilc_nodes_III, :Experiment,:Result) # ,marker=(0.2,:blue,stroke(0))
	# #@df ExpResults_ilc_nodes_I_III dotplot!(plt_ilc_nodes_I_III, :Experiment,:Result) # ,marker=(0.2,:blue,stroke(0))
	# #ylabel!("Control energy [Ws]")
	# #title!("ILC nodes")
	# savefig(plt_ilc_nodes_III, "$dir/plots/$(date)/violin_ilc_nodes_III.pdf")
	#
	#
	# # plt_ilc_nodes_I_III = plot(legend=false,ylims=())
	# # @df ExpResults_ilc_nodes_I_III violin!(plt_ilc_nodes_I_III, :Experiment,:Result) # ,marker=(0.2,:blue,stroke(0))
	# # #@df ExpResults_ilc_nodes_I_III dotplot!(plt_ilc_nodes_I_III, :Experiment,:Result) # ,marker=(0.2,:blue,stroke(0))
	# # ylabel!("Control energy [Ws]")
	# # title!("ILC nodes")
	# # savefig(plt_ilc_nodes_I_III, "$dir/plots/$(date)/violin_ilc_nodes_I_III.pdf")
	#
	# plt_ilc_nodes_II = plot(legend=false)
	# @df ExpResults_ilc_nodes_II violin!(plt_ilc_nodes_II, :Experiment,:Result) # ,marker=(0.2,:blue,stroke(0))
	# #@df ExpResults_ilc_nodes_II dotplot!(plt_ilc_nodes_II, :Experiment,:Result) # ,marker=(0.2,:blue,stroke(0))
	# ylabel!("Control energy [Ws]")
	# #title!("ILC nodes")
	# savefig(plt_ilc_nodes_II, "$dir/plots/$(date)/violin_ilc_nodes_II.pdf")
	#
	# plt_ilc_nodes_IV = plot(legend=false)
	# @df ExpResults_ilc_nodes_IV violin!(plt_ilc_nodes_IV, :Experiment,:Result) # ,marker=(0.2,:blue,stroke(0))
	# #@df ExpResults_ilc_nodes_IV dotplot!(plt_ilc_nodes_IV, :Experiment,:Result) # ,marker=(0.2,:blue,stroke(0))
	# #ylabel!("Control energy [Ws]")
	# #title!("ILC nodes")
	# savefig(plt_ilc_nodes_IV, "$dir/plots/$(date)/violin_ilc_nodes_IV.pdf")
	#
	#
	# ################################################
	#
	# #ExpResults_ilc_nodes_0abs = DataFrame(Experiment = "no ilc", Result = pnceld_ilc_0abs)
	# # ExpResults_ilc_nodes_I_IIIabs = DataFrame(Experiment = "local ilc all nodes", Result = pnceld_ilc_1abs)
	# # ExpResults_ilc_nodes_IIabs = DataFrame(Experiment = "ilc at vc \n+ average neighbors", Result = pnceld_ilc_2abs)
	# # ExpResults_ilc_nodes_IVabs = DataFrame(Experiment = "ilc at 50% random \n + 3 average random", Result = pnceld_ilc_4abs)
	# #
	# # for el in pnceld_ilc_3
	# #     push!(ExpResults_ilc_nodes_I_IIIabs, ("local ilc at vc", el))
	# # end
	# #
	# # plt_ilc_nodes_0abs = plot(legend=false)
	# #
	# # plt_ilc_nodes_I_IIIabs = plot(legend=false)
	# # @df ExpResults_ilc_nodes_I_IIIabs violin!(plt_ilc_nodes_I_IIIabs, :Experiment,:Result) # ,marker=(0.2,:blue,stroke(0))
	# # @df ExpResults_ilc_nodes_I_IIIabs dotplot!(plt_ilc_nodes_I_IIIabs, :Experiment,:Result) # ,marker=(0.2,:blue,stroke(0))
	# # ylabel!("Control energy [Ws]")
	# # title!("ILC nodes")
	# # savefig(plt_ilc_nodes_I_IIIabs, "$dir/plots/$(date)/violin_ilc_nodes_I_IIIabs.pdf")
	# #
	# # plt_ilc_nodes_IIabs = plot(legend=false)
	# # @df ExpResults_ilc_nodes_IIabs violin!(plt_ilc_nodes_IIabs, :Experiment,:Result) # ,marker=(0.2,:blue,stroke(0))
	# # @df ExpResults_ilc_nodes_IIabs dotplot!(plt_ilc_nodes_IIabs, :Experiment,:Result) # ,marker=(0.2,:blue,stroke(0))
	# # ylabel!("Control energy [Ws]")
	# # #title!("ILC nodes")
	# # savefig(plt_ilc_nodes_IIabs, "$dir/plots/$(date)/violin_ilc_nodes_IIabs.pdf")
	# #
	# # plt_ilc_nodes_IVabs = plot(legend=false)
	# # @df ExpResults_ilc_nodes_IVabs violin!(plt_ilc_nodes_IVabs, :Experiment,:Result) # ,marker=(0.2,:blue,stroke(0))
	# # @df ExpResults_ilc_nodes_IVabs dotplot!(plt_ilc_nodes_IVabs, :Experiment,:Result) # ,marker=(0.2,:blue,stroke(0))
	# # ylabel!("Control energy [Ws]")
	# # #title!("ILC nodes")
	# # savefig(plt_ilc_nodes_IVabs, "$dir/plots/$(date)/violin_ilc_nodes_IVabs.pdf")
	#
	#
	# l = @layout [a b; c d]
	# plt_ilc_nodes= plot(plt_ilc_nodes_I,plt_ilc_nodes_III,plt_ilc_nodes_II,plt_ilc_nodes_IV, layout = l)
	# savefig(plt_ilc_nodes, "$dir/plots/$(date)/violin_ilc_nodes.pdf")
	#
	#
	# ####################################
	# ####################################
	# ####################################
	#
	#
	ExpResults_nonilc_nodes_0 = DataFrame(Experiment = "0", Result = pnceld_nonilc_0)
	ExpResults_nonilc_nodes_III = DataFrame(Experiment = "III", Result = pnceld_nonilc_3)
	ExpResults_nonilc_nodes_II = DataFrame(Experiment = "II", Result = pnceld_nonilc_2)
	ExpResults_nonilc_nodes_IV = DataFrame(Experiment = "IV", Result = pnceld_nonilc_4)

	begin
		plt_ex1 = plot(legend=false)#, yscale=:log10)#, ylims=(5e6, 3e7))
	    @df ExpResults_nonilc_nodes_0 boxplot!(plt_ex1, :Experiment,:Result, whisker_width=1) # ,marker=(0.2,:blue,stroke(0))
		ylabel!("Control energy [Ws]")

	    plt_ex0 = plot(legend=false)#, yscale=:log10)#, ylims=(1e4, 2e5))
		#@df ExpResults_all_nodes_0 boxplot!(plt_ex0, :Experiment,:Result) #,fill=(0,0.5,:orange)) # ,marker=(0.2,:blue,stroke(0))
		@df ExpResults_nonilc_nodes_II boxplot!(plt_ex0, :Experiment,:Result) #,fill=(0,0.5,:orange)) # ,marker=(0.2,:blue,stroke(0))
		@df ExpResults_nonilc_nodes_II boxplot!(plt_ex0, :Experiment,:Result) #,fill=(0,0.5,:orange)) # ,marker=(0.2,:blue,stroke(0))
		@df ExpResults_nonilc_nodes_II boxplot!(plt_ex0, :Experiment,:Result) #,fill=(0,0.5,:orange)) # ,marker=(0.2,:blue,stroke(0))
	    @df ExpResults_nonilc_nodes_III boxplot!(plt_ex0, :Experiment,:Result) # ,marker=(0.2,:blue,stroke(0))
	    @df ExpResults_nonilc_nodes_IV boxplot!(plt_ex0, :Experiment,:Result,  fill=(0,0.5,:orange)) # ,marker=(0.2,:blue,stroke(0))

	    l = @layout [a{0.2w} b{0.7w}]
	    plt_ex = plot(plt_ex1, plt_ex0, layout = l)
	end
	savefig(plt_ex, "$dir/plots/$(date)/control_nonILC_nodes_boxplot_lambda_$lam.pdf")

	#
	# plt_nonilc_nodes_0 = plot(legend=false)
	# @df ExpResults_nonilc_nodes_0 violin!(plt_nonilc_nodes_0, :Experiment,:Result, fill=(0,0.5,:orange)) # ,marker=(0.2,:blue,stroke(0))
	# #@df ExpResults_nonilc_nodes_0 dotplot!(plt_nonilc_nodes_0, :Experiment,:Result) # ,marker=(0.2,:blue,stroke(0))
	# ylabel!("Control energy [Ws]")
	# title!("Non-ILC nodes")
	# savefig(plt_nonilc_nodes_0, "$dir/plots/$(date)/violin_nonilc_nodes_0.pdf")
	#
	#
	# plt_nonilc_nodes_III = plot(legend=false,ylims=(3e4,7e4))
	# @df ExpResults_nonilc_nodes_III violin!(plt_nonilc_nodes_III, :Experiment,:Result) # ,marker=(0.2,:blue,stroke(0))
	# #@df ExpResults_nonilc_nodes_III dotplot!(plt_nonilc_nodes_III, :Experiment,:Result) # ,marker=(0.2,:blue,stroke(0))
	# #ylabel!("Control energy [Ws]")
	# #title!("ILC nodes")
	# savefig(plt_nonilc_nodes_III, "$dir/plots/$(date)/violin_nonilc_nodes_III.pdf")
	#
	# plt_nonilc_nodes_II = plot(legend=false,ylims=(3e4,7e4))
	# @df ExpResults_nonilc_nodes_II violin!(plt_nonilc_nodes_II, :Experiment,:Result) # ,marker=(0.2,:blue,stroke(0))
	# #@df ExpResults_nonilc_nodes_II dotplot!(plt_nonilc_nodes_II, :Experiment,:Result) # ,marker=(0.2,:blue,stroke(0))
	# ylabel!("Control energy [Ws]")
	# #title!("ILC nodes")
	# savefig(plt_nonilc_nodes_II, "$dir/plots/$(date)/violin_nonilc_nodes_II.pdf")
	#
	# plt_nonilc_nodes_IV = plot(legend=false,ylims=(3e4,7e4))
	# @df ExpResults_nonilc_nodes_IV violin!(plt_nonilc_nodes_IV, :Experiment,:Result) # ,marker=(0.2,:blue,stroke(0))
	# #@df ExpResults_nonilc_nodes_IV dotplot!(plt_nonilc_nodes_IV, :Experiment,:Result) # ,marker=(0.2,:blue,stroke(0))
	# #ylabel!("Control energy [Ws]")
	# #title!("ILC nodes")
	# savefig(plt_nonilc_nodes_IV, "$dir/plots/$(date)/violin_nonilc_nodes_IV.pdf")
	#
	# ##################################################################
	# #
	# # ExpResults_nonilc_nodes_0abs = DataFrame(Experiment = "no ilc", Result = pnceld_nonilc_0abs)
	# # ExpResults_nonilc_nodes_IIIabs = DataFrame(Experiment = "local ilc at vc", Result = pnceld_nonilc_3abs)
	# # ExpResults_nonilc_nodes_IIabs = DataFrame(Experiment = "ilc at vc \n+ average neighbors", Result = pnceld_nonilc_2abs)
	# # ExpResults_nonilc_nodes_IVabs = DataFrame(Experiment = "ilc at 50% random \n + 3 average random", Result = pnceld_nonilc_4abs)
	# #
	# #
	# #
	# # plt_nonilc_nodes_0abs = plot(legend=false)
	# # @df ExpResults_nonilc_nodes_0abs violin!(plt_nonilc_nodes_0abs, :Experiment,:Result) # ,marker=(0.2,:blue,stroke(0))
	# # @df ExpResults_nonilc_nodes_0abs dotplot!(plt_nonilc_nodes_0abs, :Experiment,:Result) # ,marker=(0.2,:blue,stroke(0))
	# # ylabel!("Control energy [Ws]")
	# # title!("Non-ILC nodes")
	# # savefig(plt_nonilc_nodes_0abs, "$dir/plots/$(date)/violin_nonilc_nodes_0abs.pdf")
	# #
	# #
	# # plt_nonilc_nodes_IIIabs = plot(legend=false)
	# # @df ExpResults_nonilc_nodes_IIIabs violin!(plt_nonilc_nodes_IIIabs, :Experiment,:Result) # ,marker=(0.2,:blue,stroke(0))
	# # @df ExpResults_nonilc_nodes_IIIabs dotplot!(plt_nonilc_nodes_IIIabs, :Experiment,:Result) # ,marker=(0.2,:blue,stroke(0))
	# # ylabel!("Control energy [Ws]")
	# # #title!("ILC nodes")
	# # savefig(plt_nonilc_nodes_IIIabs, "$dir/plots/$(date)/violin_nonilc_nodes_IIIabs.pdf")
	# #
	# # plt_nonilc_nodes_IIabs = plot(legend=false)
	# # @df ExpResults_nonilc_nodes_IIabs violin!(plt_nonilc_nodes_IIabs, :Experiment,:Result) # ,marker=(0.2,:blue,stroke(0))
	# # @df ExpResults_nonilc_nodes_IIabs dotplot!(plt_nonilc_nodes_IIabs, :Experiment,:Result) # ,marker=(0.2,:blue,stroke(0))
	# # ylabel!("Control energy [Ws]")
	# # #title!("ILC nodes")
	# # savefig(plt_nonilc_nodes_IIabs, "$dir/plots/$(date)/violin_nonilc_nodes_IIabs.pdf")
	# #
	# # plt_nonilc_nodes_IVabs = plot(legend=false)
	# # @df ExpResults_nonilc_nodes_IVabs violin!(plt_nonilc_nodes_IVabs, :Experiment,:Result) # ,marker=(0.2,:blue,stroke(0))
	# # @df ExpResults_nonilc_nodes_IVabs dotplot!(plt_nonilc_nodes_IVabs, :Experiment,:Result) # ,marker=(0.2,:blue,stroke(0))
	# # ylabel!("Control energy [Ws]")
	# # #title!("ILC nodes")
	# # savefig(plt_nonilc_nodes_IVabs, "$dir/plots/$(date)/violin_nonilc_nodes_IVabs.pdf")
	# dummy = plot(legend=false)
	#
	# l = @layout [a b; c d]
	# plt_nonilc_nodes= plot(plt_nonilc_nodes_0,plt_nonilc_nodes_III,plt_nonilc_nodes_II,plt_nonilc_nodes_IV, layout = l)
	# savefig(plt_nonilc_nodes, "$dir/plots/$(date)/violin_nonilc_nodes.pdf")
	#
	# ##################################################################
	# #################################################################
	# ##############################################################
	#
	#
	#
	ExpResults_all_nodes_0 = DataFrame(Experiment = "0", Result = pnceld_nonilc_0)
	ExpResults_all_nodes_I = DataFrame(Experiment = "I", Result = pnceld_ilc_1)
	ExpResults_all_nodes_III = DataFrame(Experiment = "III", Result = pnceld_all_3)
	ExpResults_all_nodes_II = DataFrame(Experiment = "II", Result = pnceld_all_2)
	ExpResults_all_nodes_IV = DataFrame(Experiment = "IV", Result = pnceld_all_4)

	begin
		plt_ex1 = plot(legend=false)#, yscale=:log10)#, ylims=(5e6, 3e7))
	    @df ExpResults_all_nodes_0 boxplot!(plt_ex1, :Experiment,:Result, whisker_width=1) # ,marker=(0.2,:blue,stroke(0))
		ylabel!("Control energy [Ws]")

	    plt_ex0 = plot(legend=false)#, yscale=:log10)#, ylims=(1e4, 2e5))
		#@df ExpResults_all_nodes_0 boxplot!(plt_ex0, :Experiment,:Result) #,fill=(0,0.5,:orange)) # ,marker=(0.2,:blue,stroke(0))
		@df ExpResults_all_nodes_I boxplot!(plt_ex0, :Experiment,:Result) #,fill=(0,0.5,:orange)) # ,marker=(0.2,:blue,stroke(0))
		@df ExpResults_all_nodes_I boxplot!(plt_ex0, :Experiment,:Result) #,fill=(0,0.5,:orange)) # ,marker=(0.2,:blue,stroke(0))
	    @df ExpResults_all_nodes_II boxplot!(plt_ex0, :Experiment,:Result) # ,marker=(0.2,:blue,stroke(0))
	    @df ExpResults_all_nodes_III boxplot!(plt_ex0, :Experiment,:Result) # ,marker=(0.2,:blue,stroke(0))
		@df ExpResults_all_nodes_IV boxplot!(plt_ex0, :Experiment,:Result, fill=(0,0.5,:orange)) # ,marker=(0.2,:blue,stroke(0))

	    l = @layout [a{0.2w} b{0.7w}]
	    plt_ex = plot(plt_ex1, plt_ex0, layout = l)
	end
	savefig(plt_ex, "$dir/plots/$(date)/control_all_nodes_boxplot_lambda_$(lam).pdf")
end
	# # for el in pnceld_all_3
	# #     push!(ExpResults_all_nodes_I_III, ("local ilc at vc", el))
	# # end
	#
	# plt_all_nodes_0 = plot(legend=false)
	# @df ExpResults_all_nodes_0 violin!(plt_all_nodes_0, :Experiment,:Result,fill=(0,0.5,:orange)) # ,marker=(0.2,:blue,stroke(0))
	# #@df ExpResults_all_nodes_0 dotplot!(plt_all_nodes_0, :Experiment,:Result) # ,marker=(0.2,:blue,stroke(0))
	# ylabel!("Control energy [Ws]")
	# title!("All nodes")
	# savefig(plt_all_nodes_0, "$dir/plots/$(date)/violin_all_nodes_0.pdf")
	#
	# plt_all_nodes_I = plot(legend=false,ylims=(14000,20000))
	# @df ExpResults_all_nodes_I violin!(plt_all_nodes_I, :Experiment,:Result,fill=(0,0.5,:orange)) # ,marker=(0.2,:blue,stroke(0))
	# #@df ExpResults_all_nodes_I_III dotplot!(plt_all_nodes_I, :Experiment,:Result) # ,marker=(0.2,:blue,stroke(0))
	# #ylabel!("Control energy [Ws]")
	# #title!("ILC nodes")
	# savefig(plt_all_nodes_I, "$dir/plots/$(date)/violin_all_nodes_I.pdf")
	#
	# plt_all_nodes_III = plot(legend=false,ylims=(14000,20000))
	# @df ExpResults_all_nodes_III violin!(plt_all_nodes_III, :Experiment,:Result) # ,marker=(0.2,:blue,stroke(0))
	# #@df ExpResults_all_nodes_I_III dotplot!(plt_all_nodes_I_III, :Experiment,:Result) # ,marker=(0.2,:blue,stroke(0))
	# ylabel!("Control energy [Ws]")
	# #title!("ILC nodes")
	# savefig(plt_all_nodes_III, "$dir/plots/$(date)/violin_all_nodes_III.pdf")
	#
	#
	# # plt_all_nodes_I_III = plot(legend=false)
	# # @df ExpResults_all_nodes_I_III violin!(plt_all_nodes_I_III, :Experiment,:Result) # ,marker=(0.2,:blue,stroke(0))
	# # @df ExpResults_all_nodes_I_III dotplot!(plt_all_nodes_I_III, :Experiment,:Result) # ,marker=(0.2,:blue,stroke(0))
	# # ylabel!("Control energy [Ws]")
	# # #title!("ILC nodes")
	# # savefig(plt_all_nodes_I_III, "$dir/plots/$(date)/violin_all_nodes_I_III.pdf")
	#
	# plt_all_nodes_II = plot(legend=false)
	# @df ExpResults_all_nodes_II violin!(plt_all_nodes_II, :Experiment,:Result) # ,marker=(0.2,:blue,stroke(0))
	# #@df ExpResults_all_nodes_II dotplot!(plt_all_nodes_II, :Experiment,:Result) # ,marker=(0.2,:blue,stroke(0))
	# #ylabel!("Control energy [Ws]")
	# #title!("ILC nodes")
	# savefig(plt_all_nodes_II, "$dir/plots/$(date)/violin_all_nodes_II.pdf")
	#
	# plt_all_nodes_IV = plot(legend=false)
	# @df ExpResults_all_nodes_IV violin!(plt_all_nodes_IV, :Experiment,:Result) # ,marker=(0.2,:blue,stroke(0))
	# #@df ExpResults_all_nodes_IV dotplot!(plt_all_nodes_IV, :Experiment,:Result) # ,marker=(0.2,:blue,stroke(0))
	# #ylabel!("Control energy [Ws]")
	# #title!("ILC nodes")
	# savefig(plt_all_nodes_IV, "$dir/plots/$(date)/violin_all_nodes_IV.pdf")
	#
	# #######################################################################
	#
	# #
	# # ExpResults_all_nodes_0abs = DataFrame(Experiment = "no ilc", Result = pnceld_nonilc_0abs)
	# # ExpResults_all_nodes_I_IIIabs = DataFrame(Experiment = "local ilc all nodes", Result = pnceld_ilc_1abs)
	# # ExpResults_all_nodes_IIabs = DataFrame(Experiment = "ilc at vc \n+ average neighbors", Result = pnceld_all_2abs)
	# # ExpResults_all_nodes_IVabs = DataFrame(Experiment = "ilc at 50% random \n + 3 average random", Result = pnceld_all_4abs)
	# #
	# # for el in pnceld_all_3abs
	# #     push!(ExpResults_all_nodes_I_IIIabs, ("local ilc at vc", el))
	# # end
	# #
	# # plt_all_nodes_0abs = plot(legend=false)
	# # @df ExpResults_all_nodes_0abs violin!(plt_all_nodes_0abs, :Experiment,:Result) # ,marker=(0.2,:blue,stroke(0))
	# # @df ExpResults_all_nodes_0abs dotplot!(plt_all_nodes_0abs, :Experiment,:Result) # ,marker=(0.2,:blue,stroke(0))
	# # ylabel!("Control energy [Ws]")
	# # title!("All nodes")
	# # savefig(plt_all_nodes_0abs, "$dir/plots/$(date)/violin_all_nodes_0abs.pdf")
	# #
	# #
	# # plt_all_nodes_I_IIIabs = plot(legend=false)
	# # @df ExpResults_all_nodes_I_IIIabs violin!(plt_all_nodes_I_IIIabs, :Experiment,:Result) # ,marker=(0.2,:blue,stroke(0))
	# # @df ExpResults_all_nodes_I_IIIabs dotplot!(plt_all_nodes_I_IIIabs, :Experiment,:Result) # ,marker=(0.2,:blue,stroke(0))
	# # ylabel!("Control energy [Ws]")
	# # #title!("ILC nodes")
	# # savefig(plt_all_nodes_I_IIIabs, "$dir/plots/$(date)/violin_all_nodes_I_IIIabs.pdf")
	# #
	# # plt_all_nodes_IIabs = plot(legend=false)
	# # @df ExpResults_all_nodes_IIabs violin!(plt_all_nodes_IIabs, :Experiment,:Result) # ,marker=(0.2,:blue,stroke(0))
	# # @df ExpResults_all_nodes_IIabs dotplot!(plt_all_nodes_IIabs, :Experiment,:Result) # ,marker=(0.2,:blue,stroke(0))
	# # ylabel!("Control energy [Ws]")
	# # #title!("ILC nodes")
	# # savefig(plt_all_nodes_IIabs, "$dir/plots/$(date)/violin_all_nodes_IIabs.pdf")
	# #
	# # plt_all_nodes_IVabs = plot(legend=false)
	# # @df ExpResults_all_nodes_IVabs violin!(plt_all_nodes_IVabs, :Experiment,:Result) # ,marker=(0.2,:blue,stroke(0))
	# # @df ExpResults_all_nodes_IVabs dotplot!(plt_all_nodes_IVabs, :Experiment,:Result) # ,marker=(0.2,:blue,stroke(0))
	# # ylabel!("Control energy [Ws]")
	# # #title!("ILC nodes")
	# # savefig(plt_all_nodes_IVabs, "$dir/plots/$(date)/violin_all_nodes_IVabs.pdf")
	# #
	#
	# l = @layout [a b; c d e]
	# plt_all_nodes= plot(plt_all_nodes_0,plt_all_nodes_I, plt_all_nodes_III,plt_all_nodes_II,plt_all_nodes_IV, layout = l)
	# savefig(plt_all_nodes, "$dir/plots/$(date)/violin_all_nodes.pdf")
