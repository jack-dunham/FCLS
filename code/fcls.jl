using Makie
using CairoMakie
using MakieExtra
using DelimitedFiles
using Statistics

function plot_magn_length(D = 9)
    read = f -> collect(eachrow(readdlm(joinpath("data/dlm", f))))

    Zs = read("Zs$D.txt")
    Ls = read("Ls$D.txt")

    Ts = read("Ts.txt")
    Xis = [i for i in 4:0.2:5]

    ZsInf = read("ZsInf.txt")
    LsInf = read("LsInf.txt")

    set_theme!(
        theme_latexfonts(),
        ScatterLines = (
            markersize = 4,
            strokewidth=1/2,
            cycle = Cycle([:color, :marker]; covary=true),
        ),
        Axis = (
            xtickalign=1,
            ytickalign=1,
            xgridvisible = false,
            ygridvisible = false,
        ),
        # palette = (color = reverse(cgrad(:viridis, 6; categorical=true)),)
   )
    pt = 4/3


    fig = Figure(size = (450,400), fontsize=8pt)
    ax2 = Axis(fig[1,1], ylabel="magnetisation", xticklabelsvisible=false)
    ax4 = Axis(fig[1,2], yticklabelsvisible=false, xticklabelsvisible=false,palette = (color = reverse(cgrad(:viridis, 6; categorical=true)),))
    ax1 = Axis(fig[2,1], xlabel=L"T", ylabel="correlation length")
    ax3 = Axis(fig[2,2], xlabel=L"T", yticklabelsvisible=false,palette = (color = reverse(cgrad(:viridis, 6; categorical=true)),))

    vlines!(ax1, 1.2737, color=(:gray,0.5), linewidth=1)
    vlines!(ax2, 1.2737, color=(:gray,0.5), linewidth=1)
    vlines!(ax3, 1.2737, color=(:gray,0.5), linewidth=1)
    vlines!(ax4, 1.2737, color=(:gray,0.5), linewidth=1)

    # cp = last.(findmax.(Ls9))
    # Lc = first.(findmax.(Ls9))
    # Tc = Ts[1][cp]
    # println(Tc)


    for (Xi,Z,L) in zip(Xis,Zs,Ls)

        ZL = map((z,l) -> z .* (getindex.(l,1)).^(1/8), Z, L)
        scatterlines!(ax3,Ts[1],L; label=L"\zeta=%$Xi")
        scatterlines!(ax4,Ts[1],Z; label=L"\zeta=10^{%$Xi}")
        # scatterlines!(ax3,Ts[1],ZL)
    end
    for (Z,L,D) in zip(ZsInf,LsInf,(6,7,8,9))
        ZL = map((z,l) -> z .* (getindex.(l,1)).^(1/8), Z, L)
        scatterlines!(ax1,Ts[1],L; label=L"D=%$D")
        scatterlines!(ax2,Ts[1],Z; label=L"D=%$D")
        # scatterlines!(ax3,Ts[1],ZL)
    end

    linkxaxes!(ax1,ax2,ax3,ax4)
    linkyaxes!(ax2,ax4)
    linkyaxes!(ax1,ax3)


    axislegend(ax2, L"\zeta = \infty";position=:lb, rowgap=-0, framevisible=false, patchsize=fill(15,4))
    axislegend(ax4, L"D = %$D";position=:lb,rowgap=-0, framevisible=false, patchsize=fill(15,4))

    current_figure()

    
    return fig
end

function plotcollapse(D = 9)
    read = f -> collect(eachrow(readdlm(joinpath("data/dlm", f))))

    Zs = read("Zs$D.txt")
    Ls = read("Ls$D.txt")

    Ts = read("Ts.txt")
    Xis = [i for i in 4:0.2:5]

    ZsInf = read("ZsInf.txt")
    LsInf = read("LsInf.txt")

    set_theme!(
        theme_latexfonts(),
        ScatterLines = (
            markersize = 4,
            strokewidth=1/2,
            cycle = Cycle([:color, :marker]; covary=true),
        ),
        Axis = (
            xtickalign=1,
            ytickalign=1,
            xgridvisible = false,
            ygridvisible = false,
        ),
        # palette = (color = reverse(cgrad(:viridis, 6; categorical=true)),)
   )
    pt = 4/3


    fig = Figure(size = (500,200), fontsize=8pt)
    ax1 = Axis(fig[1,1], xlabel=L"T", ylabel=L"m \xi_{D}^{\tilde{\beta}/\nu}", )
    ax2 = Axis(fig[1,2], xlabel=L"T", ylabel=L"m \xi_{\zeta}^{\tilde{\beta}/\nu}", palette = (color = reverse(cgrad(:viridis, 6; categorical=true)),))

    vlines!(ax1, 1.2737, color=(:gray,0.5), linewidth=1)
    vlines!(ax2, 1.2737, color=(:gray,0.5), linewidth=1)

    # cp = last.(findmax.(Ls9))
    # Lc = first.(findmax.(Ls9))
    # Tc = Ts[1][cp]
    # println(Tc)


    for (Z,L,D) in zip(ZsInf,LsInf,(6,7,8,9))
        ZL = map((z,l) -> z .* (getindex.(l,1)).^(1/8), Z, L)
        scatterlines!(ax1,Ts[1],ZL; label=L"D=%$D")
    end
    for (Xi,Z,L) in zip(Xis,Zs,Ls)
        ZL = map((z,l) -> z .* (getindex.(l,1)).^(1/8), Z, L)
        scatterlines!(ax2,Ts[1],ZL; label=L"\zeta=10^{%$Xi}")
        # scatterlines!(ax3,Ts[1],ZL)
    end

    linkyaxes!(ax1,ax2)

    axislegend(ax1, L"\zeta = \infty";position=:lb, rowgap=-0, framevisible=false, patchsize=fill(15,4))
    axislegend(ax2, L"D = %$D";position=:lb,rowgap=-0, framevisible=false, patchsize=fill(15,4))

    current_figure()

    return fig
end

function plotvars()
    read = f -> collect(eachrow(readdlm(joinpath("data/dlm", f))))

    Ts = read("Ts.txt")
    Xis = [i for i in 4:0.2:5]

    ZsInf = read("ZsInf.txt")
    LsInf = read("LsInf.txt")

    set_theme!(
        theme_latexfonts(),
        ScatterLines = (
            markersize = 5,
            # strokewidth=1/2,
            strokewidth=0,
            cycle = Cycle([:color, :marker]; covary=true),
        ),
        Axis = (
            xtickalign=1,
            ytickalign=1,
            xgridvisible = false,
            ygridvisible = false,
            yscale=log10,
        ),
        # palette = (color = reverse(cgrad(:viridis, 6; categorical=true)),)
   )
    pt = 4/3


    fig = Figure(size = (400,200), fontsize=8pt)
    ax = Axis(fig[1:2,1], xlabel=L"T", ylabel=L"\mathrm{Var.}")

    # cp = last.(findmax.(Ls9))
    # Lc = first.(findmax.(Ls9))
    # Tc = Ts[1][cp]
    # println(Tc)

    # linkyaxes!(ax...)

    rowgap!(fig.layout, 0)

    Tcs = []
    Err = []
    scs = []
    for (i,D) in enumerate((6,7,8,9))
        ZLXi = []

        Zs = read("Zs$D.txt")
        Ls = read("Ls$D.txt")

        for (Xi,Z,L) in zip(Xis,Zs,Ls)
            ZL = map((z,l) -> z .* (getindex.(l,1)).^(1/8), Z, L)
            push!(ZLXi, ZL)
            # scatterlines!(ax3,Ts[1],ZL)
        end

        vars = [var(getindex.(ZLXi, i)) for i in 1:length(Ts[1])]

        sc = scatterlines!(ax,Ts[1],vars; label=L"\zeta \text{finite}, D=%$D", color=(Makie.wong_colors()[i],0.5), markercolor=(Makie.wong_colors()[i],1), marker=Cycled(i))
        push!(scs,sc)

        minv = last(findmin(vars))

        next_smallest = Ts[1][minv + 1] < Ts[1][minv - 1] ? Ts[1][minv + 1] : Ts[1][minv - 1]

        push!(Tcs, mean([Ts[1][minv], next_smallest]))
        push!(Err, abs(Ts[1][minv] - next_smallest)/2)
    end

    ZLD = []

    for (Z,L,D) in zip(ZsInf,LsInf,(6,7,8,9))
        ZL = map((z,l) -> z .* (getindex.(l,1)).^(1/8), Z, L)
        push!(ZLD, ZL)
    end

    vlines!(ax, 1.2737, color=(:gray,0.5), linewidth=1)

    vars = [var(getindex.(ZLD, i)) for i in 1:length(Ts[1])]
    scinf = scatterlines!(ax,Ts[1],vars; label=L"\zeta=\infty", color=(:black,1), markercolor=(:black,0), marker=Cycled(5), strokewidth=1/2)


        minv = last(findmin(vars))

        next_smallest = Ts[1][minv + 1] < Ts[1][minv - 1] ? Ts[1][minv + 1] : Ts[1][minv - 1]

        push!(Tcs, mean([Ts[1][minv], next_smallest]))
        push!(Err, abs(Ts[1][minv] - next_smallest)/2)
    # linkyaxes!(ax1,ax2)
    Legend(fig[1:2,2],
           [scs, [scinf]],
           [[L"D = %$D" for D in (6,7,8,9)], [L"\zeta = \infty"]], 
           [L"\zeta \in [10^{4.0},10^{4.2},\dots,10^{5.0}]",L"D \in [6,\dots,9]"], 
           rowgap=-5, titlegap=1, groupgap=10, framevisible=false, gridshalign=:left,titlehalign=:left)

    # axislegend(ax1, L"\zeta = \infty";position=:lb, rowgap=-0, framevisible=false, patchsize=fill(15,4))
    # axislegend(ax2, L"D = %$D";position=:lb,rowgap=-0, framevisible=false, patchsize=fill(15,4))

    current_figure()

    return fig, Tcs, Err
end
