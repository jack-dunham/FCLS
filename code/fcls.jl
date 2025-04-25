using TimeEvolutionPEPO
using TensorKit
using Makie
using CairoMakie
using MakieExtra
using DelimitedFiles
using Statistics
using DrWatson
using Tar
using DataFrames
using LsqFit
using EasyFit
using LaTeXStrings
using Measurements

function extract()
    for (root, _, files) in walkdir(datadir("raw"))
        for file in files
            name, ext = splitext(file)
            if ext == ".gz"
                @info "Extracting" filename = file
                tarball = joinpath(root, file)

                outpath = datadir("raw", "extract", splitext(name) |> first)

                if isdir(outpath)
                    @info "Skipping"
                    continue
                end

                Tar.list(`gzcat $tarball`)
                Tar.extract(`gzcat $tarball`, outpath)
            end
        end
    end


    df = collect_results!(
        datadir("pro", "results.jld2"),
        datadir("raw", "extract"),
        rinclude=[r".*obs"],
        black_list=["dm", "psi", "time"],
        subfolders=true,
    )
    return df
end

function plot_magn_length(D=9)
    read = f -> collect(eachrow(readdlm(joinpath("data/dlm", f))))

    Zs = read("Zs$D.txt")
    Ls = read("Ls$D.txt")

    Ts = read("Ts.txt")
    Xis = [i for i in 4:0.2:5]

    ZsInf = read("ZsInf.txt")
    LsInf = read("LsInf.txt")

    set_theme!(
        theme_latexfonts(),
        ScatterLines=(
            markersize=4,
            strokewidth=1 / 2,
            cycle=Cycle([:color, :marker]; covary=true),
        ),
        Axis=(
            xtickalign=1,
            ytickalign=1,
            xgridvisible=false,
            ygridvisible=false,
        ),
        # palette = (color = reverse(cgrad(:viridis, 6; categorical=true)),)
    )
    pt = 4 / 3


    fig = Figure(size=(450, 400), fontsize=8pt)
    ax2 = Axis(fig[1, 1], ylabel="magnetisation", xticklabelsvisible=false)
    ax4 = Axis(fig[1, 2], yticklabelsvisible=false, xticklabelsvisible=false, palette=(color=reverse(cgrad(:viridis, 6; categorical=true)),))
    ax1 = Axis(fig[2, 1], xlabel=L"T", ylabel="correlation length")
    ax3 = Axis(fig[2, 2], xlabel=L"T", yticklabelsvisible=false, palette=(color=reverse(cgrad(:viridis, 6; categorical=true)),))

    vlines!(ax1, 1.2737, color=(:gray, 0.5), linewidth=1)
    vlines!(ax2, 1.2737, color=(:gray, 0.5), linewidth=1)
    vlines!(ax3, 1.2737, color=(:gray, 0.5), linewidth=1)
    vlines!(ax4, 1.2737, color=(:gray, 0.5), linewidth=1)

    # cp = last.(findmax.(Ls9))
    # Lc = first.(findmax.(Ls9))
    # Tc = Ts[1][cp]
    # println(Tc)


    for (Xi, Z, L) in zip(Xis, Zs, Ls)

        ZL = map((z, l) -> z .* (getindex.(l, 1)) .^ (1 / 8), Z, L)
        scatterlines!(ax3, Ts[1], L; label=L"\zeta=%$Xi")
        scatterlines!(ax4, Ts[1], Z; label=L"\zeta=10^{%$Xi}")
        # scatterlines!(ax3,Ts[1],ZL)
    end
    for (Z, L, D) in zip(ZsInf, LsInf, (6, 7, 8, 9))
        ZL = map((z, l) -> z .* (getindex.(l, 1)) .^ (1 / 8), Z, L)
        scatterlines!(ax1, Ts[1], L; label=L"D=%$D")
        scatterlines!(ax2, Ts[1], Z; label=L"D=%$D")
        # scatterlines!(ax3,Ts[1],ZL)
    end

    linkxaxes!(ax1, ax2, ax3, ax4)
    linkyaxes!(ax2, ax4)
    linkyaxes!(ax1, ax3)


    axislegend(ax2, L"\zeta = \infty"; position=:lb, rowgap=-0, framevisible=false, patchsize=fill(15, 4))
    axislegend(ax4, L"D = %$D"; position=:lb, rowgap=-0, framevisible=false, patchsize=fill(15, 4))

    current_figure()


    return fig
end

function plotcollapse(D=9)
    read = f -> collect(eachrow(readdlm(joinpath("data/dlm", f))))

    Zs = read("Zs$D.txt")
    Ls = read("Ls$D.txt")

    Ts = read("Ts.txt")
    Xis = [i for i in 4:0.2:5]

    ZsInf = read("ZsInf.txt")
    LsInf = read("LsInf.txt")

    set_theme!(
        theme_latexfonts(),
        ScatterLines=(
            markersize=4,
            strokewidth=1 / 2,
            cycle=Cycle([:color, :marker]; covary=true),
        ),
        Axis=(
            xtickalign=1,
            ytickalign=1,
            xgridvisible=false,
            ygridvisible=false,
        ),
        # palette = (color = reverse(cgrad(:viridis, 6; categorical=true)),)
    )
    pt = 4 / 3


    fig = Figure(size=(500, 200), fontsize=8pt)
    ax1 = Axis(fig[1, 1], xlabel=L"T", ylabel=L"m \xi_{D}^{\tilde{\beta}/\nu}",)
    ax2 = Axis(fig[1, 2], xlabel=L"T", ylabel=L"m \xi_{\zeta}^{\tilde{\beta}/\nu}", palette=(color=reverse(cgrad(:viridis, 6; categorical=true)),))

    vlines!(ax1, 1.2737, color=(:gray, 0.5), linewidth=1)
    vlines!(ax2, 1.2737, color=(:gray, 0.5), linewidth=1)

    # cp = last.(findmax.(Ls9))
    # Lc = first.(findmax.(Ls9))
    # Tc = Ts[1][cp]
    # println(Tc)


    for (Z, L, D) in zip(ZsInf, LsInf, (6, 7, 8, 9))
        ZL = map((z, l) -> z .* (getindex.(l, 1)) .^ (1 / 8), Z, L)
        scatterlines!(ax1, Ts[1], ZL; label=L"D=%$D")
    end
    for (Xi, Z, L) in zip(Xis, Zs, Ls)
        ZL = map((z, l) -> z .* (getindex.(l, 1)) .^ (1 / 8), Z, L)
        scatterlines!(ax2, Ts[1], ZL; label=L"\zeta=10^{%$Xi}")
        # scatterlines!(ax3,Ts[1],ZL)
    end

    linkyaxes!(ax1, ax2)

    axislegend(ax1, L"\zeta = \infty"; position=:lb, rowgap=-0, framevisible=false, patchsize=fill(15, 4))
    axislegend(ax2, L"D = %$D"; position=:lb, rowgap=-0, framevisible=false, patchsize=fill(15, 4))

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
        ScatterLines=(
            markersize=5,
            # strokewidth=1/2,
            strokewidth=0,
            cycle=Cycle([:color, :marker]; covary=true),
        ),
        Axis=(
            xtickalign=1,
            ytickalign=1,
            xgridvisible=false,
            ygridvisible=false,
            yscale=log10,
        ),
        # palette = (color = reverse(cgrad(:viridis, 6; categorical=true)),)
    )
    pt = 4 / 3


    fig = Figure(size=(400, 200), fontsize=8pt)
    ax = Axis(fig[1:2, 1], xlabel=L"T", ylabel=L"\mathrm{Var.}")

    # cp = last.(findmax.(Ls9))
    # Lc = first.(findmax.(Ls9))
    # Tc = Ts[1][cp]
    # println(Tc)

    # linkyaxes!(ax...)

    rowgap!(fig.layout, 0)

    Tcs = []
    Err = []
    scs = []
    for (i, D) in enumerate((6, 7, 8, 9))
        ZLXi = []

        Zs = read("Zs$D.txt")
        Ls = read("Ls$D.txt")

        for (Xi, Z, L) in zip(Xis, Zs, Ls)
            ZL = map((z, l) -> z .* (getindex.(l, 1)) .^ (1 / 8), Z, L)
            push!(ZLXi, ZL)
            # scatterlines!(ax3,Ts[1],ZL)
        end

        vars = [var(getindex.(ZLXi, i)) for i in 1:length(Ts[1])]

        sc = scatterlines!(ax, Ts[1], vars; label=L"\zeta \text{finite}, D=%$D", color=(Makie.wong_colors()[i], 0.5), markercolor=(Makie.wong_colors()[i], 1), marker=Cycled(i))
        push!(scs, sc)

        minv = last(findmin(vars))

        next_smallest = Ts[1][minv+1] < Ts[1][minv-1] ? Ts[1][minv+1] : Ts[1][minv-1]

        push!(Tcs, mean([Ts[1][minv], next_smallest]))
        push!(Err, abs(Ts[1][minv] - next_smallest) / 2)
    end

    ZLD = []

    for (Z, L, D) in zip(ZsInf, LsInf, (6, 7, 8, 9))
        ZL = map((z, l) -> z .* (getindex.(l, 1)) .^ (1 / 8), Z, L)
        push!(ZLD, ZL)
    end

    vlines!(ax, 1.2737, color=(:gray, 0.5), linewidth=1)

    vars = [var(getindex.(ZLD, i)) for i in 1:length(Ts[1])]
    scinf = scatterlines!(ax, Ts[1], vars; label=L"\zeta=\infty", color=(:black, 1), markercolor=(:black, 0), marker=Cycled(5), strokewidth=1 / 2)


    minv = last(findmin(vars))

    next_smallest = Ts[1][minv+1] < Ts[1][minv-1] ? Ts[1][minv+1] : Ts[1][minv-1]

    push!(Tcs, mean([Ts[1][minv], next_smallest]))
    push!(Err, abs(Ts[1][minv] - next_smallest) / 2)
    # linkyaxes!(ax1,ax2)
    Legend(fig[1:2, 2],
        [scs, [scinf]],
        [[L"D = %$D" for D in (6, 7, 8, 9)], [L"\zeta = \infty"]],
        [L"\zeta \in [10^{4.0},10^{4.2},\dots,10^{5.0}]", L"D \in [6,\dots,9]"],
        rowgap=-5, titlegap=1, groupgap=10, framevisible=false, gridshalign=:left, titlehalign=:left)

    # axislegend(ax1, L"\zeta = \infty";position=:lb, rowgap=-0, framevisible=false, patchsize=fill(15,4))
    # axislegend(ax2, L"D = %$D";position=:lb,rowgap=-0, framevisible=false, patchsize=fill(15,4))

    current_figure()

    return fig, Tcs, Err
end


function plotquantum(df=extract(); D=9)

    ax = Axis(Figure()[1, 1])

    df = filter(row -> row["D"] == D, df)

    for gdf in groupby(df, "ζ")

        sort!(gdf, "χ")

        xi = Measurement{Float64}[]
        xi32 = Float64[]
        yv32 = Float64[]

        yv = Measurement{Float64}[]

        for i in 1:50

            delt = Float64[]
            lens = Float64[]


            local xi32i
            local m32i

            for gdfi in groupby(gdf, "χ")
                push!(delt, gdfi[1, "δ"][i])
                push!(lens, gdfi[1, "ξ"][i])
                xi32i = gdfi[1, "ξ"][i]

                m32i = real(gdfi[1, "m"][i])
            end

            p, _, s = myfitnonlinear(delt, lens)

            xi_extrap = p[3] ± s[3]

            push!(xi, xi_extrap)
            push!(xi32, xi32i)

            push!(yv, m32i * xi_extrap^(1 / 8))
            push!(yv32, m32i * xi32i^(1 / 8))
        end

        # lines!(ax, xi)
        # lines!(ax, xi32; color=:black)
        ζ = gdf[1, "ζ"]
        time = load(gdf[1, "path"], "time")
        # lines!(ax, 1 ./ time, yv; label=L"\zeta = %$ζ")
        # lines!(ax, 1 ./ time, yv32; label=L"\zeta = %$ζ")

        if ζ == Inf
            label = L"\zeta = \infty"
        else
            label = L"\zeta = %$(log10(ζ))"
        end

        lines!(ax, 1 ./ time, xi; label=label)
        lines!(ax, 1 ./ time, xi32; linestyle=:dash)

        vlines!(ax, 1.2737, color=(:gray, 0.5), linewidth=1)
    end
    axislegend(ax)
    return
end

function myfitnonlinear(x, y, wt=nothing; verbose=true)

    init = fitlinear(x, y)

    verbose && println(init)

    fitfunc(t, p) = @. p[2] * t^p[1] + p[3]

    p0 = [1.0, 1.0, 1.0]

    if isnothing(wt)
        fit = curve_fit(fitfunc, x, y, [1.0, init.a, init.b])
        # fit = curve_fit(fitfunc, x, y, p0)
    else
        fit = curve_fit(fitfunc, x, y, wt, [1.0, init.a, init.b])
    end

    println(fit.param)
    println(stderror(fit))

    return fit.param, t -> fitfunc(t, fit.param), stderror(fit)
    # return fit.param, t -> fitfunc(t, fit.param), (0, 0, 0)
end
