plt_ex0 = plot(legend=false,  dpi=150, #yscale=:log10,
		   ytickfontsize=18,
		   xtickfontsize=14,
		   guidefontsize=20,
		   legendfontsize=10)
ylabel!("Control energy (ILC) [Ws]")
xlabel!("Number of ILC nodes")


plt_ex1 = plot(legend=false, dpi=150,
				   ytickfontsize=18,
				   xtickfontsize=14,
				   guidefontsize=20,
				   legendfontsize=10)
ylabel!("Control energy (non-ILC) [Ws]")
xlabel!("Number of ILC nodes")


plt_ex2 = plot(legend=false, dpi=150,
						   ytickfontsize=18,
						   xtickfontsize=14,
						   guidefontsize=20,
						   legendfontsize=10)
ylabel!("Control energy (all) [Ws]")
xlabel!("Number of ILC nodes")




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

	pnceld_ilc_3 = mean.([[p[3] for p in res_tlIII_ir.u][i][vc3[index_calc(i,batch_size)]] for i in (m-1)*batch_size+1:m*batch_size])
	pnceld_nonilc_3 = mean.([[p[3] for p in res_tlIII_ir.u][i][nonvc3[index_calc(i,batch_size)]] for i in (m-1)*batch_size+1:m*batch_size])
	pnceld_all_3 = mean.([[p[3] for p in res_tlIII_ir.u][i][1:N] for i in (m-1)*batch_size+1:m*batch_size])

	################################################################
	# ILC NODES
	################################

	ExpResults_ilc_nodes_III = DataFrame(Experiment = "$m", Result = pnceld_ilc_3)
	    @df ExpResults_ilc_nodes_III boxplot!(plt_ex0, :Experiment,:Result) # ,marker=(0.2,:blue,stroke(0))
	    plt_ilc = plot(plt_ex0)#,plt_ex1, layout = l)


#########################
# NON ILC NODES
#########################
	if m != N
		ExpResults_nonilc_nodes_III = DataFrame(Experiment = "$m", Result = pnceld_nonilc_3)
	    @df ExpResults_nonilc_nodes_III boxplot!(plt_ex1, :Experiment,:Result) # ,marker=(0.2,:blue,stroke(0))
	    plt_nonilc = plot(plt_ex1)
	end

	#####################################################
	# ALL NODES
	#######################################################

	ExpResults_all_nodes_III = DataFrame(Experiment = "$m", Result = pnceld_all_3)
	@df ExpResults_all_nodes_III boxplot!(plt_ex2, :Experiment,:Result) # ,marker=(0.2,:blue,stroke(0))
	plt_all = plot(plt_ex2)

end

savefig(plt_ex0, "$dir/plots/$(date)/ir_control_ILC_nodes_boxplot_$(lamda).pdf")
savefig(plt_ex1, "$dir/plots/$(date)/ir_control_nonILC_nodes_boxplot_$(lamda).pdf")
savefig(plt_ex2, "$dir/plots/$(date)/ir_control_all_nodes_boxplot_$(lamda).pdf")
