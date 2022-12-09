for m = 1:length(lambda_lst) # lambda = 0.2:0.1:1
        println("m ",m)
        lam = lambda_lst[m]

    Ex_all_nodes_0 = DataFrame(Experiment = "0", Result = [p[2] for p in res_tl0.u][(m-1)*batch_size+1:m*batch_size]* 100. /(3600. * 24. * obs_days)) # no ilc
    #println(Ex_all_nodes_0)
    Ex_all_nodes_I = DataFrame(Experiment = "I", Result = [p[2] for p in res_tlI.u][(m-1)*batch_size+1:m*batch_size] * 100. /(3600. * 24. * obs_days)) #local ilc all nodes
    Ex_all_nodes_III = DataFrame(Experiment = "III", Result = [p[2] for p in res_tlIII.u][(m-1)*batch_size+1:m*batch_size] * 100. /(3600. * 24. * obs_days)) # local ilc at vc
    Ex_all_nodes_II = DataFrame(Experiment = "II", Result = [p[2] for p in res_tlII.u][(m-1)*batch_size+1:m*batch_size] * 100. /(3600. * 24. * obs_days)) #ilc at vc \n+ neighbor. comm.
    Ex_all_nodes_IV = DataFrame(Experiment = "IV", Result = [p[2] for p in res_tlIV.u][(m-1)*batch_size+1:m*batch_size] * 100. /(3600. * 24. * obs_days)) #ilc at 50% random \n + 3 random comm.

    begin
        plt_ex0 = plot(legend=false,dpi=150,#,ylims=(-100,2100),
        		   ytickfontsize=14,
                   xtickfontsize=14,
        		   guidefontsize=14,
        		   legendfontsize=3)
        @df Ex_all_nodes_0 boxplot!(plt_ex0, :Experiment,:Result, fill=(0,0.5,:blue), whisker_width=0.5, color = :blue) #,fill=(0,0.5,:orange)) # ,marker=(0.2,:blue,stroke(0))
        @df Ex_all_nodes_I boxplot!(plt_ex0, :Experiment,:Result, fill=(0,0.5,:red), color = :red) # ,marker=(0.2,:blue,stroke(0))
        @df Ex_all_nodes_II boxplot!(plt_ex0, :Experiment,:Result, fill=(0,0.5,:green), color = :green) # ,marker=(0.2,:blue,stroke(0))
        @df Ex_all_nodes_III boxplot!(plt_ex0, :Experiment,:Result, fill=(0,0.5,:purple),color = :purple) # ,marker=(0.2,:blue,stroke(0))
        @df Ex_all_nodes_IV boxplot!(plt_ex0, :Experiment,:Result,fill=(0,0.5,:orange), whisker_width=1, color = :orange) # ,marker=(0.2,:blue,stroke(0))
        ylabel!("Exceedance [s]",  guidefontsize=14)

        # plt_ex1 = plot(legend=false,dpi=150,#,ylims=(-100,2100),
        # 		   ytickfontsize=14,
        #            xtickfontsize=14,
        # 		   guidefontsize=14,
        # 		   legendfontsize=3, yscale=:log10)#, ylims=(1e4, 1e6))
        # @df Ex_all_nodes_0 boxplot!(plt_ex1, :Experiment,:Result, fill=(0,0.5,:blue), whisker_width=0.5, color = :blue) #,fill=(0,0.5,:orange)) # ,marker=(0.2,:blue,stroke(0))
        # ylabel!("Exceedance [s]",  guidefontsize=14)
        # l = @layout [a{0.15w} b]
        # plt_ex = plot(plt_ex1,plt_ex0,layout = l)
    end

    savefig(plt_ex, "$dir/plots/$(date)/ex_boxplot_lambda_$(lam).pdf")

end
