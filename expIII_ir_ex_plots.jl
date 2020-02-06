plt_ex0 = plot(legend=false,dpi=150,ylims=(-100,2500),
                ytickfontsize=12,
                xtickfontsize=12,
                guidefontsize=12,
                legendfontsize=12,
                margin = 5Plots.mm,
                xticks = (2:2:48, string.(2:2:24)))
xlabel!("Number of ILC nodes")
ylabel!("Exceedance [s]")


for m = 1:length(ilc_lst) # ilc_lst = 1:N
        println("ilc_nodes ",m)

    Ex_all_nodes_III = DataFrame(Experiment = m, Result = [p[2] for p in res_tlIII_ir.u][(m-1)*batch_size+1:m*batch_size]) # local ilc at m random nodes
    @df Ex_all_nodes_III boxplot!(plt_ex0, :Experiment,:Result) # ,marker=(0.2,:blue,stroke(0))

end

savefig(plt_ex0, "$dir/plots/$(date)/ir_ex_boxplot_ilc_$(lamda).pdf")
