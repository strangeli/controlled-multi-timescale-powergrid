for m = 1:length(lambda_lst) # lambda = 0.2:0.1:1
        println("m ",m)
        lam = lambda_lst[m]

    Max_freq_all_nodes_0I = DataFrame(Experiment = "0", Result = [p[1] for p in res_tl0.u][(m-1)*batch_size+1:m*batch_size])
    Max_freq_all_nodes_I = DataFrame(Experiment = "I", Result = [p[1] for p in res_tlI.u][(m-1)*batch_size+1:m*batch_size])
    Max_freq_all_nodes_III = DataFrame(Experiment = "III", Result = [p[1] for p in res_tlIII.u][(m-1)*batch_size+1:m*batch_size])
    Max_freq_all_nodes_II = DataFrame(Experiment = "II", Result = [p[1] for p in res_tlII.u][(m-1)*batch_size+1:m*batch_size])
    # Max_freq_all_nodes_II_III = DataFrame(Experiment = "II: ilc at vc \n+ neighbor comm.", Result = [p[1] for p in res_tlII.u])
    Max_freq_all_nodes_IV = DataFrame(Experiment = "IV", Result = [p[1] for p in res_tlIV.u][(m-1)*batch_size+1:m*batch_size])


    begin
        plt_ex0 = plot(legend=false, dpi=150, #ylims = (0.15,0.3),
        		   ytickfontsize=14,
                   xtickfontsize=14,
        		   guidefontsize=14,
        		   legendfontsize=10)
        @df Max_freq_all_nodes_0I boxplot!(plt_ex0, :Experiment,:Result, fill=(0,0.5,:blue)) #,fill=(0,0.5,:orange)) # ,marker=(0.2,:blue,stroke(0))
        @df Max_freq_all_nodes_I boxplot!(plt_ex0, :Experiment,:Result, fill=(0,0.5,:red)) # ,marker=(0.2,:blue,stroke(0))
        @df Max_freq_all_nodes_II boxplot!(plt_ex0, :Experiment,:Result, fill=(0,0.5,:green)) # ,marker=(0.2,:blue,stroke(0))
        @df Max_freq_all_nodes_III boxplot!(plt_ex0, :Experiment,:Result, fill=(0,0.5,:purple)) # ,marker=(0.2,:blue,stroke(0))
        @df Max_freq_all_nodes_IV boxplot!(plt_ex0, :Experiment,:Result,fill=(0,0.5,:orange)) # ,marker=(0.2,:blue,stroke(0))
        ylabel!("Max frequency deviation [rad/s]",guidefontsize=14)

        plt_ex1 = plot(legend=false)#,ylims = (0.27,0.63))
        @df Max_freq_all_nodes_IV boxplot!(plt_ex1, :Experiment,:Result,fill=(0,0.5,:orange), whisker_width=1) # ,marker=(0.2,:blue,stroke(0))
        #ylabel!("Exceedance [s]")

        l = @layout [a b{0.2w}]
        plt_ex = plot(plt_ex0,plt_ex1, layout = l)
    end
    savefig(plt_ex0, "$dir/plots/$(date)/freq_dev_boxplot_lambda_$(lam).pdf")
end



    #
    #
    # for el in [p[1] for p in res_tlI.u]
    # 	push!(Max_freq_all_nodes_0I, ("I: local ilc all nodes", el))
    # end
    #
    # # for el in [p[1] for p in res_tlIII.u]
    # # 	push!(Max_freq_all_nodes_II_III, ("III: local ilc at vc", el))
    # # end
    #
    # plt_max_freq0I = plot(legend=false,ylims = (0.18,0.3))
    # @df Max_freq_all_nodes_0I violin!(plt_max_freq0I, :Experiment,:Result,fill=(0,0.5,:orange)) # ,marker=(0.2,:blue,stroke(0))
    # #@df Max_freq_all_nodes_0I dotplot!(plt_max_freq0I, :Experiment,:Result) # ,marker=(0.2,:blue,stroke(0))
    # ylabel!("Max frequency deviation [pu]")
    # title!("Maximum frequency deviation")
    # savefig(plt_max_freq0I, "$dir/plots/$(date)/violin_max_freq0I.png")
    #
    # #
    # # plt_max_freq0 = plot(legend=false)
    # # @df Max_freq_all_nodes_0 violin!(plt_max_freq0, :Experiment,:Result) # ,marker=(0.2,:blue,stroke(0))
    # # @df Max_freq_all_nodes_0 dotplot!(plt_max_freq0, :Experiment,:Result) # ,marker=(0.2,:blue,stroke(0))
    # # ylabel!("Maximum frequency deviation [pu]")
    # # title!("Max freq deviation")
    # # savefig(plt_max_freq0, "$dir/plots/$(date)/violin_max_freq0.png")
    # #
    # # plt_Max_freq_I= plot(legend=false)
    # # @df Max_freq_all_nodes_I violin!(plt_Max_freq_I, :Experiment,:Result) # ,marker=(0.2,:blue,stroke(0))
    # # @df Max_freq_all_nodes_I dotplot!(plt_Max_freq_I, :Experiment,:Result) # ,marker=(0.2,:blue,stroke(0))
    # # #ylabel!("Maximum frequency deviation [pu]")
    # # savefig(plt_Max_freq_I, "$dir/plots/$(date)/violin_Max_freq_I.png")
    #
    # # plt_Max_freq_II_III= plot(legend=false,ylims = (0.18,0.3))
    # # @df Max_freq_all_nodes_II_III violin!(plt_Max_freq_II_III, :Experiment,:Result) # ,marker=(0.2,:blue,stroke(0))
    # # #@df Max_freq_all_nodes_II_III dotplot!(plt_Max_freq_II_III, :Experiment,:Result) # ,marker=(0.2,:blue,stroke(0))
    # # # ylabel!("Maximum frequency deviation [pu]")
    # # savefig(plt_Max_freq_II_III, "$dir/plots/$(date)/violin_Max_freq_II_III.png")
    #
    # plt_Max_freq_III= plot(legend=false,ylims = (0.18,0.3))
    # @df Max_freq_all_nodes_III violin!(plt_Max_freq_III, :Experiment,:Result) # ,marker=(0.2,:blue,stroke(0))
    # # @df Max_freq_all_nodes_III dotplot!(plt_Max_freq_III, :Experiment,:Result) # ,marker=(0.2,:blue,stroke(0))
    # # ylabel!("Maximum frequency deviation [pu]")
    # savefig(plt_Max_freq_III, "$dir/plots/$(date)/violin_Max_freq_III.png")
    #
    # plt_Max_freq_II = plot(legend=false,ylims = (0.18,0.3))
    # @df Max_freq_all_nodes_II violin!(plt_Max_freq_II, :Experiment,:Result) # ,marker=(0.2,:blue,stroke(0))
    # # @df Max_freq_all_nodes_II dotplot!(plt_Max_freq_II, :Experiment,:Result) # ,marker=(0.2,:blue,stroke(0))
    # #ylabel!("Maximum frequency deviation [pu]")
    # savefig(plt_Max_freq_II, "$dir/plots/$(date)/violin_Max_freq_II.png")
    #
    # plt_Max_freq_IV = plot(legend=false,ylims = (0.1,0.7))
    # @df Max_freq_all_nodes_IV violin!(plt_Max_freq_IV, :Experiment,:Result) # ,marker=(0.2,:blue,stroke(0))
    # #@df Max_freq_all_nodes_IV dotplot!(plt_Max_freq_IV, :Experiment,:Result) # ,marker=(0.2,:blue,stroke(0))
    # #ylabel!("Maximum frequency deviation [pu]")
    # savefig(plt_Max_freq_IV, "$dir/plots/$(date)/violin_Max_freq_IV.png")
    #
    # ###########################
    # #
    # # Max_freq_all_nodes_0abs = DataFrame(Experiment = "0: no ilc", Result = [p[1] for p in res_tl0abs.u])
    # # Max_freq_all_nodes_I_IIIabs = DataFrame(Experiment = "I: local ilc all nodes", Result = [p[1] for p in res_tlIabs.u])
    # # Max_freq_all_nodes_IIabs = DataFrame(Experiment = "II: ilc at vc \n+ average neighbors", Result = [p[1] for p in res_tlIIabs.u])
    # # Max_freq_all_nodes_IVabs = DataFrame(Experiment = "IV: ilc at 50% random \n + 3 averag random", Result = [p[1] for p in res_tlIVabs.u])
    # #
    # #
    # # for el in [p[1] for p in res_tlIIIabs.u]
    # # 	push!(Max_freq_all_nodes_I_IIIabs, ("III: local ilc at vc", el))
    # # end
    # #
    # #
    # #
    # # plt_max_freq0abs = plot(legend=false)
    # # @df Max_freq_all_nodes_0abs violin!(plt_max_freq0abs, :Experiment,:Result) # ,marker=(0.2,:blue,stroke(0))
    # # @df Max_freq_all_nodes_0abs dotplot!(plt_max_freq0abs, :Experiment,:Result) # ,marker=(0.2,:blue,stroke(0))
    # # #ylabel!("Maximum frequency deviation [pu]")
    # # title!("Max freq deviation (abs)")
    # # savefig(plt_max_freq0abs, "$dir/plots/$(date)/violin_max_freq0abs.png")
    # #
    # # plt_Max_freq_I_IIIabs= plot(legend=false)
    # # @df Max_freq_all_nodes_I_IIIabs violin!(plt_Max_freq_I_IIIabs, :Experiment,:Result) # ,marker=(0.2,:blue,stroke(0))
    # # @df Max_freq_all_nodes_I_IIIabs dotplot!(plt_Max_freq_I_IIIabs, :Experiment,:Result) # ,marker=(0.2,:blue,stroke(0))
    # # ylabel!("Maximum frequency deviation [pu]")
    # # savefig(plt_Max_freq_I_IIIabs, "$dir/plots/$(date)/violin_Max_freq_I_IIIabs.png")
    # #
    # #
    # # plt_Max_freq_IIabs = plot(legend=false)
    # # @df Max_freq_all_nodes_IIabs violin!(plt_Max_freq_IIabs, :Experiment,:Result) # ,marker=(0.2,:blue,stroke(0))
    # # @df Max_freq_all_nodes_IIabs dotplot!(plt_Max_freq_IIabs, :Experiment,:Result) # ,marker=(0.2,:blue,stroke(0))
    # # #ylabel!("Maximum frequency deviation [pu]")
    # # savefig(plt_Max_freq_IIabs, "$dir/plots/$(date)/violin_Max_freq_IIabs.png")
    # #
    # #
    # # plt_Max_freq_IVabs = plot(legend=false)
    # # @df Max_freq_all_nodes_IVabs violin!(plt_Max_freq_IVabs, :Experiment,:Result) # ,marker=(0.2,:blue,stroke(0))
    # # @df Max_freq_all_nodes_IVabs dotplot!(plt_Max_freq_IVabs, :Experiment,:Result) # ,marker=(0.2,:blue,stroke(0))
    # # #ylabel!("Maximum frequency deviation [pu]")
    # # savefig(plt_Max_freq_IVabs, "$dir/plots/$(date)/violin_Max_freq_IVabs.png")
    # #
    #
    # l = @layout [a; b c d]
    # plt_max_freq = plot(plt_max_freq0I,plt_Max_freq_III, plt_Max_freq_II, plt_Max_freq_IV, layout = l)
    # savefig(plt_max_freq, "$dir/plots/$(date)/violin_max_freq.png")
