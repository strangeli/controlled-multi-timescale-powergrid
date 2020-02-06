for m = 1:length(lambda_lst) # lambda = 0.2:0.1:1
        println("m ",m)
        lam = lambda_lst[m]

    Ex_all_nodes_0 = DataFrame(Experiment = "0", Result = [p[2] for p in res_tl0.u][(m-1)*batch_size+1:m*batch_size]) # no ilc
    Ex_all_nodes_I = DataFrame(Experiment = "I", Result = [p[2] for p in res_tlI.u][(m-1)*batch_size+1:m*batch_size]) #local ilc all nodes
    Ex_all_nodes_III = DataFrame(Experiment = "III", Result = [p[2] for p in res_tlIII.u][(m-1)*batch_size+1:m*batch_size]) # local ilc at vc
    Ex_all_nodes_II = DataFrame(Experiment = "II", Result = [p[2] for p in res_tlII.u][(m-1)*batch_size+1:m*batch_size]) #ilc at vc \n+ neighbor. comm.
    Ex_all_nodes_IV = DataFrame(Experiment = "IV", Result = [p[2] for p in res_tlIV.u][(m-1)*batch_size+1:m*batch_size]) #ilc at 50% random \n + 3 random comm.

    begin
        plt_ex0 = plot(legend=false,dpi=150,#,ylims=(-100,2100),
        		   ytickfontsize=14,
                   xtickfontsize=14,
        		   guidefontsize=14,
        		   legendfontsize=3)
        @df Ex_all_nodes_I boxplot!(plt_ex0, :Experiment,:Result, fill=(0,0.5,:red), color = :red) # ,marker=(0.2,:blue,stroke(0))
        @df Ex_all_nodes_II boxplot!(plt_ex0, :Experiment,:Result, fill=(0,0.5,:green), color = :green) # ,marker=(0.2,:blue,stroke(0))
        @df Ex_all_nodes_III boxplot!(plt_ex0, :Experiment,:Result, fill=(0,0.5,:purple),color = :purple) # ,marker=(0.2,:blue,stroke(0))
        @df Ex_all_nodes_IV boxplot!(plt_ex0, :Experiment,:Result,fill=(0,0.5,:orange), whisker_width=1, color = :orange) # ,marker=(0.2,:blue,stroke(0))


        plt_ex1 = plot(legend=false,dpi=150,#,ylims=(-100,2100),
        		   ytickfontsize=14,
                   xtickfontsize=14,
        		   guidefontsize=14,
        		   legendfontsize=3, yscale=:log10)#, ylims=(1e4, 1e6))
        @df Ex_all_nodes_0 boxplot!(plt_ex1, :Experiment,:Result, fill=(0,0.5,:blue), whisker_width=0.5, color = :blue) #,fill=(0,0.5,:orange)) # ,marker=(0.2,:blue,stroke(0))
        ylabel!("Exceedance [s]",  guidefontsize=14)
        l = @layout [a{0.15w} b]
        plt_ex = plot(plt_ex1,plt_ex0,layout = l)
    end

    savefig(plt_ex, "$dir/plots/$(date)/ex_boxplot_lambda_$(lam).pdf")
end




# for el in [p[2] for p in res_tlIII.u]
# 	push!(Ex_all_nodes_I_III, ("III: local ilc at vc", el))
# end

# plt_ex0 = plot(legend=false,ylims=(0,2000))
# @df Ex_all_nodes_0 boxplot!(plt_ex0, :Experiment,:Result,fill=(0,0.5,:orange)) # ,marker=(0.2,:blue,stroke(0))
# #@df Ex_all_nodes_0 dotplot!(plt_ex0, :Experiment,:Result) # ,marker=(0.2,:blue,stroke(0))
# ylabel!("Exceedance [s]")
# title!("Exceedance")
# savefig(plt_ex0, "$dir/plots/$(date)/violin_ex0.pdf")
#
# plt_ex_I= plot(legend=false,ylims=(0,2000))
# @df Ex_all_nodes_I violin!(plt_ex_I, :Experiment,:Result,fill=(0,0.5,:orange)) # ,marker=(0.2,:blue,stroke(0))
# #@df Ex_all_nodes_I dotplot!(plt_ex_I, :Experiment,:Result) # ,marker=(0.2,:blue,stroke(0))
# #ylabel!("Exceedance [s]")
# savefig(plt_ex_I, "$dir/plots/$(date)/violin_ex_I.pdf")

#
# plt_ex_III= plot(legend=false,ylims=(0,2000))
# @df Ex_all_nodes_III violin!(plt_ex_III, :Experiment,:Result) # ,marker=(0.2,:blue,stroke(0))
# #@df Ex_all_nodes_III dotplot!(plt_ex_III, :Experiment,:Result) # ,marker=(0.2,:blue,stroke(0))
# ylabel!("Exceedance [s]")
# savefig(plt_ex_III, "$dir/plots/$(date)/violin_ex_III.pdf")
#
# plt_ex_II = plot(legend=false,ylims=(0,2000))
# @df Ex_all_nodes_II violin!(plt_ex_II, :Experiment,:Result) # ,marker=(0.2,:blue,stroke(0))
# #@df Ex_all_nodes_II dotplot!(plt_ex_II, :Experiment,:Result) # ,marker=(0.2,:blue,stroke(0))
# #ylabel!("Exceedance [s]")
# savefig(plt_ex_II, "$dir/plots/$(date)/violin_ex_II.pdf")
#
# plt_ex_IV = plot(legend=false)
# @df Ex_all_nodes_IV violin!(plt_ex_IV, :Experiment,:Result) # ,marker=(0.2,:blue,stroke(0))
# #@df Ex_all_nodes_IV dotplot!(plt_ex_IV, :Experiment,:Result) # ,marker=(0.2,:blue,stroke(0))
# #ylabel!("Exceedance [s]")
# savefig(plt_ex_IV, "$dir/plots/$(date)/violin_ex_IV.pdf")
#
# ###########################
# #
# # Ex_all_nodes_0abs = DataFrame(Experiment = "0: no ilc", Result = [p[2] for p in res_tl0abs.u])
# # Ex_all_nodes_I_IIIabs = DataFrame(Experiment = "I: local ilc all nodes", Result = [p[2] for p in res_tlIabs.u])
# # Ex_all_nodes_IIabs = DataFrame(Experiment = "II: ilc at vc \n+ average neighbors", Result = [p[2] for p in res_tlIIabs.u])
# # Ex_all_nodes_IVabs = DataFrame(Experiment = "IV: ilc at 50% random \n + 3 averag random", Result = [p[2] for p in res_tlIVabs.u])
# #
# #
# # for el in [p[2] for p in res_tlIIIabs.u]
# # 	push!(Ex_all_nodes_I_IIIabs, ("III: local ilc at vc", el))
# # end
# #
# #
# #
# # plt_ex0abs = plot(legend=false)
# # @df Ex_all_nodes_0abs violin!(plt_ex0abs, :Experiment,:Result) # ,marker=(0.2,:blue,stroke(0))
# # @df Ex_all_nodes_0abs dotplot!(plt_ex0abs, :Experiment,:Result) # ,marker=(0.2,:blue,stroke(0))
# # ylabel!("Exceedance [s]")
# # title!("Exceedance (abs)")
# # savefig(plt_ex0abs, "$dir/plots/$(date)/violin_ex0abs.pdf")
# #
# # plt_ex_I_IIIabs= plot(legend=false)
# # @df Ex_all_nodes_I_IIIabs violin!(plt_ex_I_IIIabs, :Experiment,:Result) # ,marker=(0.2,:blue,stroke(0))
# # @df Ex_all_nodes_I_IIIabs dotplot!(plt_ex_I_IIIabs, :Experiment,:Result) # ,marker=(0.2,:blue,stroke(0))
# # ylabel!("Exceedance [s]")
# # savefig(plt_ex_I_IIIabs, "$dir/plots/$(date)/violin_ex_I_IIIabs.pdf")
# #
# #
# # plt_ex_IIabs = plot(legend=false)
# # @df Ex_all_nodes_IIabs violin!(plt_ex_IIabs, :Experiment,:Result) # ,marker=(0.2,:blue,stroke(0))
# # @df Ex_all_nodes_IIabs dotplot!(plt_ex_IIabs, :Experiment,:Result) # ,marker=(0.2,:blue,stroke(0))
# # ylabel!("Exceedance [s]")
# # savefig(plt_ex_IIabs, "$dir/plots/$(date)/violin_ex_IIabs.pdf")
# #
# #
# # plt_ex_IVabs = plot(legend=false)
# # @df Ex_all_nodes_IVabs violin!(plt_ex_IVabs, :Experiment,:Result) # ,marker=(0.2,:blue,stroke(0))
# # @df Ex_all_nodes_IVabs dotplot!(plt_ex_IVabs, :Experiment,:Result) # ,marker=(0.2,:blue,stroke(0))
# # ylabel!("Exceedance [s]")
# # savefig(plt_ex_IVabs, "$dir/plots/$(date)/violin_ex_IVabs.pdf")
#
#
# l = @layout [a b; c d e]
# plt_ex = plot(plt_ex0,plt_ex_I,plt_ex_III,plt_ex_II, plt_ex_IV, layout = l)
# savefig(plt_ex, "$dir/plots/$(date)/violin_ex.pdf")
