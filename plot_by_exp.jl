
plt_cost0 = plot(legend=false, margin=8Plots.mm)#, yscale=:log10)#, ylims=(1e4, 1e6))

for m = 1:length(lambda_lst) # lambda = 0.2:0.1:1
        println("m ",m)
        lam = lambda_lst[m]

    Cost_all_nodes_0 = DataFrame(Experiment = "$(lam)", Result = mean.([p[end] for p in res_tl0.u][(m-1)*batch_size+1:m*batch_size])) # no ilc

    @df Cost_all_nodes_0 boxplot!(plt_cost0, :Experiment,:Result,fill=(0,0.5,:orange), whisker_width=1) # ,marker=(0.2,:blue,stroke(0))

end

ylabel!("Cost (mean)")
xlabel!("Lambda")
title!("0: no ILC")

savefig("$dir/plots/$(date)/box_exp0_cost.png")

#######################################################

plt_cost1 = plot(legend=false, margin=8Plots.mm)#, yscale=:log10)#, ylims=(1e4, 1e6))

for m = 1:length(lambda_lst) # lambda = 0.2:0.1:1
        println("m ",m)
        lam = lambda_lst[m]

    Cost_all_nodes_1 = DataFrame(Experiment = "$(lam)", Result = mean.([p[end] for p in res_tlI.u][(m-1)*batch_size+1:m*batch_size])) # no ilc

    @df Cost_all_nodes_1 boxplot!(plt_cost1, :Experiment,:Result,fill=(0,0.5,:blue), whisker_width=1) # ,marker=(0.2,:blue,stroke(0))

end

ylabel!("Cost (mean)")
xlabel!("Lambda")
title!("I: local ILC all nodes")

savefig("$dir/plots/$(date)/box_expI_cost.png")

#######################################################


plt_cost2 = plot(legend=false, margin=8Plots.mm)#, yscale=:log10)#, ylims=(1e4, 1e6))

for m = 1:length(lambda_lst) # lambda = 0.2:0.1:1
        println("m ",m)
        lam = lambda_lst[m]

    Cost_all_nodes_2 = DataFrame(Experiment = "$(lam)", Result = mean.([p[end] for p in res_tlII.u][(m-1)*batch_size+1:m*batch_size])) # no ilc

    @df Cost_all_nodes_2 boxplot!(plt_cost2, :Experiment,:Result,fill=(0,0.5,:green), whisker_width=1) # ,marker=(0.2,:blue,stroke(0))

end

ylabel!("Cost (mean)")
xlabel!("Lambda")
title!("II: local ILC at vc + neighbor. com")

savefig("$dir/plots/$(date)/box_expII_cost.png")

########################################################

plt_cost3 = plot(legend=false, margin=8Plots.mm)#, yscale=:log10)#, ylims=(1e4, 1e6))

for m = 1:length(lambda_lst) # lambda = 0.2:0.1:1
        println("m ",m)
        lam = lambda_lst[m]

    Cost_all_nodes_3 = DataFrame(Experiment = "$(lam)", Result = mean.([p[end] for p in res_tlIII.u][(m-1)*batch_size+1:m*batch_size])) # no ilc

    @df Cost_all_nodes_3 boxplot!(plt_cost3, :Experiment,:Result,fill=(0,0.5,:red), whisker_width=1) # ,marker=(0.2,:blue,stroke(0))

end

ylabel!("Cost (mean)")
xlabel!("Lambda")
title!("III: local ILC at vc (no com)")

savefig("$dir/plots/$(date)/box_expIII_cost.png")

##########################################################

plt_cost4 = plot(legend=false, margin=8Plots.mm)#, yscale=:log10)#, ylims=(1e4, 1e6))

for m = 1:length(lambda_lst) # lambda = 0.2:0.1:1
        println("m ",m)
        lam = lambda_lst[m]

    Cost_all_nodes_4 = DataFrame(Experiment = "$(lam)", Result = mean.([p[end] for p in res_tlIV.u][(m-1)*batch_size+1:m*batch_size])) # no ilc

    @df Cost_all_nodes_4 boxplot!(plt_cost4, :Experiment,:Result,fill=(0,0.5,:purple), whisker_width=1) # ,marker=(0.2,:blue,stroke(0))

end

ylabel!("Cost (mean)")
xlabel!("Lambda")
title!("IV: local ILC random 50% + rand com")

savefig("$dir/plots/$(date)/box_expIV_cost.png")

########################################################################
y0 = zeros(length(lambda_lst))
y1 = zeros(length(lambda_lst))
y2 = zeros(length(lambda_lst))
y3 = zeros(length(lambda_lst))
y4 = zeros(length(lambda_lst))

for m = 1:length(lambda_lst) # lambda = 0.2:0.1:1
        println("m ",m)
        lam = lambda_lst[m]
    y0[m] = mean(sum.([p[end] for p in res_tl0.u][(m-1)*batch_size+1:m*batch_size]))
    y1[m] = mean(sum.([p[end] for p in res_tlI.u][(m-1)*batch_size+1:m*batch_size]))
    y2[m] = mean(sum.([p[end] for p in res_tlII.u][(m-1)*batch_size+1:m*batch_size]))
    y3[m] = mean(sum.([p[end] for p in res_tlIII.u][(m-1)*batch_size+1:m*batch_size]))
    y4[m] = mean(sum.([p[end] for p in res_tlIV.u][(m-1)*batch_size+1:m*batch_size]))
end

plta = plot(collect(0.2:0.1:1), y0, label = "0: no ILC", margin=8Plots.mm, dpi=600,
            ytickfontsize=12,
            xtickfontsize=12,
            guidefontsize=12,
            legendfontsize=10,lw=3,color=:blue, lt=:scatter,alpha = 0.4)
pltb = plot(collect(0.2:0.1:1), y1, lw = 3, color = :red, ls = :dot, alpha = 0.4, label = "", ytickfontsize=12,
xtickfontsize=12,
guidefontsize=12,
legendfontsize=10)
plot!(collect(0.2:0.1:1), y1, label = "I: local ILC all nodes",
           margin=8Plots.mm, dpi=600,
           ytickfontsize=12, yscale=:log10,
           xtickfontsize=12,
           guidefontsize=12,
           legendfontsize=10,lw=3, color = :red, alpha = 0.4,  lt=:scatter)
           plot!(collect(0.2:0.1:1), y2, lw = 3, color = :green, alpha = 0.4, label = "")
plot!(collect(0.2:0.1:1), y2, label = "II: local ILC at vc + neighbor. com", lw = 3, lt=:scatter, m = :cross, color = :green, alpha = 0.4)
plot!(collect(0.2:0.1:1), y3, label = "", ls = :dash, lw=3,  color = :purple, alpha = 0.4)
plot!(collect(0.2:0.1:1), y3, label = "III: local ILC at vc", ls = :dash, lw=3,  lt=:scatter, m = :hex, color = :purple, alpha = 0.4)
plot!(collect(0.2:0.1:1), y4, label = "", ls=:dashdot, lw = 3, color = :orange, alpha = 0.4)
plot!(collect(0.2:0.1:1), y4, label = "IV: local ILC rand 50% + rand com", ls=:dashdot, lw = 3, m = :star7, color = :orange,lt=:scatter, alpha = 0.4)
xlabel!(L"\lambda [-]")
ylabel!("Average total cost [a.u.]")


#l = @layout [a; b]
#plt = plot(plta,pltb, layout = l)


savefig("$dir/plots/$(date)/box_cost_all_exp_total.pdf")
