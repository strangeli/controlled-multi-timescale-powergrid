plt_c0 = plot(legend=false, dpi=150,
					ytickfontsize=14,
					xtickfontsize=14,
					guidefontsize=14,
					legendfontsize=14,
					margin = 5Plots.mm,
					xticks = (2:2:48, string.(2:2:24)))
ylabel!("Overall cost (ILC nodes) [a.u.]")
xlabel!("Number of ILC nodes")


plt_c1 = plot(legend=false, dpi=150,
				ytickfontsize=14,
				xtickfontsize=14,
				guidefontsize=14,
				legendfontsize=14,
				margin = 5Plots.mm,
				xticks = (2:2:48, string.(2:2:24)))
ylabel!("Overall cost (non-ILC nodes) [a.u.]")
xlabel!("Number of ILC nodes")


plt_c2 = plot(legend=false, dpi=150,
						   ytickfontsize=14,
						   xtickfontsize=14,
						   guidefontsize=14,
						   legendfontsize=14,
	   					margin = 5Plots.mm,
						 xticks = (2:2:48, string.(2:2:24)))
ylabel!("Overall cost (all nodes) [a.u.]", guidefontsize=14)
xlabel!("Number of ILC nodes", guidefontsize=14)




for m = 1:length(ilc_lst)
        println("m ",m)

	# [(m-1)*batch_size+1:m*batch_size]
	vc3 = [p[7] for p in res_tlIII_ir.u][(m-1)*batch_size+1:m*batch_size]


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


	graph_lst = Graph.([p[5] for p in res_tlIII_ir.u][(m-1)*batch_size+1:m*batch_size])


	##############################################
	######### BOX PLOTS #######################
	##############################################

	# [(m-1)*batch_size+1:m*batch_size]

	cost_ilc_3 = sum.([[p[end] for p in res_tlIII_ir.u][i][vc3[index_calc(i,batch_size)]] for i in (m-1)*batch_size+1:m*batch_size])
	cost_nonilc_3 = sum.([[p[end] for p in res_tlIII_ir.u][i][nonvc3[index_calc(i,batch_size)]] for i in (m-1)*batch_size+1:m*batch_size])
	cost_all_3 = sum.([[p[end] for p in res_tlIII_ir.u][i][1:N] for i in (m-1)*batch_size+1:m*batch_size])

	################################################################
	# ILC NODES
	################################

	ExpResults_ilc_nodes_III = DataFrame(Experiment = m, Result = cost_ilc_3)
	    @df ExpResults_ilc_nodes_III boxplot!(plt_c0, :Experiment,:Result,whisker_width=0.5) # ,marker=(0.2,:blue,stroke(0))
	    plt_ilc = plot(plt_c0)#,plt_c1, layout = l)


#########################
# NON ILC NODES
#########################
	if m != N
		ExpResults_nonilc_nodes_III = DataFrame(Experiment = m, Result = cost_nonilc_3)
	    @df ExpResults_nonilc_nodes_III boxplot!(plt_c1,:Result) # ,marker=(0.2,:blue,stroke(0))
	    plt_nonilc = plot(plt_c1)
	end

	#####################################################
	# ALL NODES
	#######################################################

	ExpResults_all_nodes_III = DataFrame(Experiment = m, Result = cost_all_3)
	@df ExpResults_all_nodes_III boxplot!(plt_c2, :Experiment,:Result) # ,marker=(0.2,:blue,stroke(0))
	plt_all = plot(plt_c2)

end

savefig(plt_c0, "$dir/plots/$(date)/ir_cost_ILC_nodes_boxplot_$(lamda).svg")
savefig(plt_c1, "$dir/plots/$(date)/ir_cost_nonILC_nodes_boxplot_$(lamda).svg")
savefig(plt_c2, "$dir/plots/$(date)/ir_cost_all_nodes_boxplot_$(lamda).svg")
