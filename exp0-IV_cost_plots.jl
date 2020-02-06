for m = 1:length(lambda_lst) # lambda = 0.2:0.1:1
        println("m ",m)
        lam = lambda_lst[m]

    Cost_all_nodes_0 = DataFrame(Experiment = "0", Result = sum.([p[end] for p in res_tl0.u][(m-1)*batch_size+1:m*batch_size])) # no ilc
    Cost_all_nodes_I = DataFrame(Experiment = "I", Result = sum.([p[end] for p in res_tlI.u][(m-1)*batch_size+1:m*batch_size])) #local ilc all nodes
    Cost_all_nodes_III = DataFrame(Experiment = "III", Result = sum.([p[end] for p in res_tlIII.u][(m-1)*batch_size+1:m*batch_size])) # local ilc at vc
    Cost_all_nodes_II = DataFrame(Experiment = "II", Result = sum.([p[end] for p in res_tlII.u][(m-1)*batch_size+1:m*batch_size])) #ilc at vc \n+ neighbor. comm.
    Cost_all_nodes_IV = DataFrame(Experiment = "IV", Result = sum.([p[end] for p in res_tlIV.u][(m-1)*batch_size+1:m*batch_size])) #ilc at 50% random \n + 3 random comm.

    begin


        plt_cost0 = plot(legend=false,dpi=150,#ylims=(-100,2100),
        		   ytickfontsize=14,
                   xtickfontsize=14,
        		   guidefontsize=14,
        		   legendfontsize=3)
        #@df Cost_all_nodes_0 boxplot!(plt_cost0, :Experiment,:Result) #,fill=(0,0.5,:orange)) # ,marker=(0.2,:blue,stroke(0))
        @df Cost_all_nodes_I boxplot!(plt_cost0, :Experiment,:Result, fill=(0,0.5,:red), color = :red) # ,marker=(0.2,:blue,stroke(0))
        @df Cost_all_nodes_II boxplot!(plt_cost0, :Experiment,:Result, fill=(0,0.5,:green), color = :green) # ,marker=(0.2,:blue,stroke(0))
        @df Cost_all_nodes_III boxplot!(plt_cost0, :Experiment,:Result, fill=(0,0.5,:purple), color = :purple) # ,marker=(0.2,:blue,stroke(0))
        @df Cost_all_nodes_IV boxplot!(plt_cost0, :Experiment,:Result, fill=(0,0.5,:orange), color = :orange)#, whisker_width=1) # ,marker=(0.2,:blue,stroke(0))

        plt_cost1 = plot(legend=false,dpi=150,#ylims=(-100,2100),
        		   ytickfontsize=14,
                   xtickfontsize=14,
        		   guidefontsize=14,
        		   legendfontsize=3)#, yscale=:log10)#, ylims=(1e4, 1e6))
        @df Cost_all_nodes_0 boxplot!(plt_cost1, :Experiment,:Result, fill=(0,0.5,:blue), color=:blue, whisker_width=0.6)#,fill=(0,0.5,:orange), whisker_width=1) # ,marker=(0.2,:blue,stroke(0))
            ylabel!("Overall cost [a.u.]")

        l = @layout [a{0.15w} b]
        plt_cost = plot(plt_cost1,plt_cost0,layout = l)
    end

    savefig(plt_cost, "$dir/plots/$(date)/cost_boxplot_lambda_$(lam).pdf")
end

################################################

for m = 1:length(lambda_lst) # lambda = 0.2:0.1:1
        println("m ",m)
        lam = lambda_lst[m]

        vc2 = [p[7] for p in res_tlII.u][(m-1)*batch_size+1:m*batch_size]
        vc3 = [p[7] for p in res_tlIII.u][(m-1)*batch_size+1:m*batch_size]
        vc4 = [p[7] for p in res_tlIV.u][(m-1)*batch_size+1:m*batch_size]

    Cost_ilc_nodes_I = DataFrame(Experiment = "I", Result = sum.([[p[end] for p in res_tlI.u][i][1:N] for i in (m-1)*batch_size+1:m*batch_size])) #local ilc all nodes
    Cost_ilc_nodes_III = DataFrame(Experiment = "III", Result = sum.([[p[end] for p in res_tlIII.u][i][vc3[index_calc(i,batch_size)]] for i in (m-1)*batch_size+1:m*batch_size])) # local ilc at vc
    Cost_ilc_nodes_II = DataFrame(Experiment = "II", Result = sum.([[p[end] for p in res_tlII.u][i][vc2[index_calc(i,batch_size)]] for i in (m-1)*batch_size+1:m*batch_size])) #ilc at vc \n+ neighbor. comm.
    Cost_ilc_nodes_IV = DataFrame(Experiment = "IV", Result = sum.([[p[end] for p in res_tlIV.u][i][vc4[index_calc(i,batch_size)]] for i in (m-1)*batch_size+1:m*batch_size]))#ilc at 50% random \n + 3 random comm.

    begin

        plt_cost_ilc = plot(legend=false,dpi=150,#ylims=(-100,2100),
        		   ytickfontsize=14,
                   xtickfontsize=14,
        		   guidefontsize=14,
        		   legendfontsize=3)
        #@df Cost_all_nodes_0 boxplot!(plt_cost0, :Experiment,:Result) #,fill=(0,0.5,:orange)) # ,marker=(0.2,:blue,stroke(0))
        @df Cost_ilc_nodes_I boxplot!(plt_cost_ilc, :Experiment,:Result, fill=(0,0.5,:red), color = :red) # ,marker=(0.2,:blue,stroke(0))
        @df Cost_ilc_nodes_II boxplot!(plt_cost_ilc, :Experiment,:Result, fill=(0,0.5,:green), color = :green) # ,marker=(0.2,:blue,stroke(0))
        @df Cost_ilc_nodes_III boxplot!(plt_cost_ilc, :Experiment,:Result, fill=(0,0.5,:purple), color = :purple) # ,marker=(0.2,:blue,stroke(0))
        @df Cost_ilc_nodes_IV boxplot!(plt_cost_ilc, :Experiment,:Result, fill=(0,0.5,:orange), color = :orange) # ,marker=(0.2,:blue,stroke(0))

    end
    ylabel!("Overall cost [a.u.]")


    savefig(plt_cost_ilc, "$dir/plots/$(date)/cost_boxplot_ilcnodes_lambda_$(lam).pdf")
end

################################################

for m = 1:length(lambda_lst) # lambda = 0.2:0.1:1
        println("m ",m)
        lam = lambda_lst[m]

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

    #Cost_nonilc_nodes_0 = DataFrame(Experiment = "0", Result = sum.([[p[end] for p in res_tl0.u][i][1:N] for i in (m-1)*batch_size+1:m*batch_size])) #local ilc all nodes
    Cost_nonilc_nodes_I = DataFrame(Experiment = "I", Result = sum.([[p[end] for p in res_tlI.u][i][1:N] for i in (m-1)*batch_size+1:m*batch_size])) #local ilc all nodes
    Cost_nonilc_nodes_III = DataFrame(Experiment = "III", Result = sum.([[p[end] for p in res_tlIII.u][i][nonvc3[index_calc(i,batch_size)]] for i in (m-1)*batch_size+1:m*batch_size])) # local ilc at vc
    Cost_nonilc_nodes_II = DataFrame(Experiment = "II", Result = sum.([[p[end] for p in res_tlII.u][i][nonvc2[index_calc(i,batch_size)]] for i in (m-1)*batch_size+1:m*batch_size])) #ilc at vc \n+ neighbor. comm.
    Cost_nonilc_nodes_IV = DataFrame(Experiment = "IV", Result = sum.([[p[end] for p in res_tlIV.u][i][nonvc4[index_calc(i,batch_size)]] for i in (m-1)*batch_size+1:m*batch_size]))#ilc at 50% random \n + 3 random comm.

    begin

        plt_cost_nonilc = plot(legend=false,dpi=150,#ylims=(-100,2100),
        		   ytickfontsize=14,
                   xtickfontsize=14,
        		   guidefontsize=14,
        		   legendfontsize=3)
        #@df Cost_nonilc_nodes_0 boxplot!(plt_cost_nonilc, :Experiment,:Result) #,fill=(0,0.5,:orange)) # ,marker=(0.2,:blue,stroke(0))
        #@df Cost_nonilc_nodes_I boxplot!(plt_cost_nonilc, :Experiment,:Result, fill=(0,0.5,:red)) # ,marker=(0.2,:blue,stroke(0))
        @df Cost_nonilc_nodes_II boxplot!(plt_cost_nonilc, :Experiment,:Result, fill=(0,0.5,:green), color = :green) # ,marker=(0.2,:blue,stroke(0))
        @df Cost_nonilc_nodes_III boxplot!(plt_cost_nonilc, :Experiment,:Result, fill=(0,0.5,:purple), color = :purple) # ,marker=(0.2,:blue,stroke(0))
        @df Cost_nonilc_nodes_IV boxplot!(plt_cost_nonilc, :Experiment,:Result, fill=(0,0.5,:orange), color = :orange) # ,marker=(0.2,:blue,stroke(0))

    end
    ylabel!("Overall cost [a.u.]")


    savefig(plt_cost_nonilc, "$dir/plots/$(date)/cost_boxplot_nonilcnodes_lambda_$(lam).pdf")
end
