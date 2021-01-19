# This file contains functions used for the production of the plots in: 
# Figure 6 of the main manuscript, as well as Supplementary Figures 4, 5, 6, and 7.
# For plotting we relied on the PGFPlotsX julia package, that essentially creates
# pgfplots code that is then compliled into the final figure. 

# The parameters used for the specific plots in our manuscript are explained within the text. 
# For simplicity, we have added a Jupyter Notebook that recreates the exact plots seen in the manuscript


using PGFPlotsX, Colors, LaTeXStrings

# Formatting instructions for the production of the final vector graphic using PGFPlotsX 
push!(PGFPlotsX.CUSTOM_PREAMBLE, raw"\usepgfplotslibrary{fillbetween}")
push!(PGFPlotsX.CUSTOM_PREAMBLE, raw"\usepackage{cmbright}")
push!(PGFPlotsX.CUSTOM_PREAMBLE, raw"\usepackage[T1]{fontenc}")
push!(PGFPlotsX.CUSTOM_PREAMBLE, raw"\DeclareFontFamily{U}{futm}{}
    \DeclareFontShape{U}{futm}{m}{n}{
        <-> s * [.92] fourier-bb
    }{}
    \DeclareSymbolFont{Ufutm}{U}{futm}{m}{n}
    \DeclareSymbolFontAlphabet{\mathbb}{Ufutm}")
push!(PGFPlotsX.CUSTOM_PREAMBLE, raw"\definecolor{colorBm}{RGB}{178,223,138}") 
push!(PGFPlotsX.CUSTOM_PREAMBLE, raw"\definecolor{colorBp}{RGB}{31,120,180}") 



function plotBmPop(xs, ys, error, title)
    ggplot2_plot_theme = @pgf {
        "no marks",
        "error bars/y dir=both",
        "error bars/y explicit",
        line_width = "1.64pt",
    }
    mycolor = raw"colorBm"

    @pgf SemiLogXAxis(
        {
            xmin = minimum(xs),
            xmax = 1.05 * maximum(xs),
            ymin = 0,
            ymax = 1.05 * maximum(maximum(ys[1] .+ error)),
            title = title,
            xlabel = "Dissemination Time", #"Cell Population",
            ylabel = L"\large$\mathbb{E}[B_{md}]$",
        },
        Plot(
            {"name path=f", color = mycolor, ggplot2_plot_theme...},
            Coordinates(xs, ys[1]),
        ),
        Plot(
            {"name path=g", color = mycolor, opacity = 0.35, ggplot2_plot_theme...},
            Coordinates(xs, ys[1] .+ error),
        ),
        Plot(
            {"name path=q", color = mycolor, opacity = 0.35, ggplot2_plot_theme...},
            Coordinates(xs, ys[1] .- error),
        ),
        Plot(
            {thick, color = mycolor, fill = mycolor, opacity = 0.35},
            raw"fill between [of=f and g]",
        ),
        Plot(
            {thick, color = mycolor, fill = mycolor, opacity = 0.35},
            raw"fill between [of=f and q]",
        ),
    )
end


function plotBp(xs, ys, error, title)

    ggplot2_plot_theme = @pgf {
        #         mark="square",
        #         mark_size = 1,
        "no marks",
        "error bars/y dir=both",
        "error bars/y explicit",
        line_width = "1.64pt",
    }

    mycolor = raw"colorBp"
    
    @pgf SemiLogXAxis(
        {
            xmin = minimum(xs),
            xmax = 1.05 * maximum(xs),
            ymin = 0,
            ymax = 1.05 * maximum(maximum(ys[1] .+ error)),
            title = title,
            xlabel = "Dissemination Time", #"Cell Population",
            ylabel = L"\large$\mathbb{E}[B_{p}]$",
        },
        # [
        #     PlotInc(
        #         {ggplot2_plot_theme..., color = colors[i]},
        #         Coordinates(xs, _y; yerror = error),
        #     ) for (i, _y) in enumerate(ys)
        # ]...,
        Plot(
            {"name path=f", color = mycolor, ggplot2_plot_theme...},
            Coordinates(xs, ys[1]),
        ),
        Plot(
            {"name path=g", color = mycolor, opacity = 0.35, ggplot2_plot_theme...},
            Coordinates(xs, ys[1] .+ error),
        ),
        Plot(
            {"name path=q", color = mycolor, opacity = 0.35, ggplot2_plot_theme...},
            Coordinates(xs, ys[1] .- error),
        ),
        Plot(
            {thick, color = mycolor, fill = mycolor, opacity = 0.35},
            raw"fill between [of=f and g]",
        ),
        Plot(
            {thick, color = mycolor, fill = mycolor, opacity = 0.35},
            raw"fill between [of=f and q]",
        ),
    )
end





function BmMutIndexPlotPrepare(detectionFrequency, MAX_MUT, mutationProbability, ρ)
    k_max = MAX_MUT
    y = [Bm(detectionFrequency, mutationProbability, ρ, k_seed) for k_seed = 1:k_max]
    ys = [y]
    xs = 1:k_max
    error = sqrt.(BmVar.(detectionFrequency, mutationProbability, ρ, 1:k_max))
    options = Dict("u" => mutationProbability, "rho" => ρ, "f" => detectionFrequency)
    return xs, ys, error, options
end


function plotBmVar(xs, ys, options::Dict{String,Float64})
  
    ggplot2_plot_theme =
        @pgf {mark = "square", mark_size = 0, mark_options = "solid", line_width = "1.64pt"}

    n = length(ys)


    # Evenly spread out colors
     colors = range(HSL(colorant"red"), stop=HSL(colorant"green"), length=n)
    @pgf SemiLogXAxis(
        {
            xmin = minimum(xs),
            xmax = 1.05 * maximum(xs),
            ymin = minimum(minimum(ys)),
            ymax = 20,
            title = L"Neutral Growth $(u = $"*
                   "$(options["u"]),"*
                    L"$f = $"*
                    "$(options["f"]))",
            xlabel = "Dissemination Time",
            ylabel = L"$\sqrt{\mathrm{Var}(B_{md}^{neutral})}$",
            legend_pos = "north west",
        },
        [
            PlotInc({ggplot2_plot_theme..., color = colors[i]}, Coordinates(xs, _y))
            for (i, _y) in enumerate(ys)
        ]...,
        Legend([L"$\rho = 0.72$", L"$\rho = 0.81$", L"$\rho = 0.90$", L"$\rho = 0.99$"]),
    )
end

# This is the old one with the error bars
function plotBmPop_(xs, ys, error, title)
    ggplot2_axis_theme = @pgf {
        tick_align = "outside",
        tick_pos = "left",
        xmajorgrids,
        x_grid_style = "white",
        ymajorgrids,
        y_grid_style = "white",
        axis_line_style = "white",
        "axis_background/.style" = {fill = "white!89.803921568627459!black"},
    }

    ggplot2_plot_theme = @pgf {
        #         mark="square",
        #         mark_size = 1,
        "no marks",
        "error bars/y dir=both",
        "error bars/y explicit",
        line_width = "1.64pt",
    }

    n = length(ys)

    # Evenly spread out colors
    colors = [LCHuv(65, 100, h) for h in range(15; stop = 360 + 15, length = n + 1)][1:n]
    @pgf SemiLogXAxis(
        #@pgf Axis(
        {
            ggplot2_axis_theme...,
            xmin = minimum(xs),
            xmax = 1.05 * maximum(xs),
            ymin = 0,
            ymax = 1.05 * maximum(maximum(ys)),
            title = title,
            xlabel = "Cell Population",
            ylabel = L"$B_m$",
        },
        [
            PlotInc(
                {ggplot2_plot_theme..., color = colors[i]},
                Coordinates(xs, _y; yerror = error),
            ) for (i, _y) in enumerate(ys)
        ]...,
    )
end


function plotBmPop(xs, ys, error, title)
    ggplot2_axis_theme = @pgf {
        tick_align = "outside",
        tick_pos = "left",
        # xmajorgrids,
        # x_grid_style = "white",
        # ymajorgrids,
        # y_grid_style = "white",
        axis_line_style = "white",
        "axis_background/.style" = {fill = "white!89.803921568627459!black"},
    }

    ggplot2_plot_theme = @pgf {
        #         mark="square",
        #         mark_size = 1,
        "no marks",
        "error bars/y dir=both",
        "error bars/y explicit",
        line_width = "1.64pt",
    }

    n = length(ys)
    n = 3
    # Evenly spread out colors
    colors = [LCHuv(65, 100, h) for h in range(15; stop = 360 + 15, length = n + 1)][1:n]
    mycolor = colors[2]
    push!(PGFPlotsX.CUSTOM_PREAMBLE, raw"\usepgfplotslibrary{fillbetween}")
    push!(PGFPlotsX.CUSTOM_PREAMBLE, raw"\usepackage{cmbright}")
    push!(PGFPlotsX.CUSTOM_PREAMBLE, raw"\usepackage[T1]{fontenc}")

    @pgf SemiLogXAxis(
        #@pgf Axis(
        {
            ggplot2_axis_theme...,
            xmin = minimum(xs),
            xmax = 1.05 * maximum(xs),
            ymin = 0,
            ymax = 1.05 * maximum(maximum(ys[1] .+ error)),
            title = title,
            xlabel = "Seeding Time", #"Cell Population",
            ylabel = L"$B_{md}$",
        },
        # [
        #     PlotInc(
        #         {ggplot2_plot_theme..., color = colors[i]},
        #         Coordinates(xs, _y; yerror = error),
        #     ) for (i, _y) in enumerate(ys)
        # ]...,
        Plot(
            {"name path=f", color = mycolor, ggplot2_plot_theme...},
            Coordinates(xs, ys[1]),
        ),
        Plot(
            {"name path=g", color = mycolor, opacity = 0.25, ggplot2_plot_theme...},
            Coordinates(xs, ys[1] .+ error),
        ),
        Plot(
            {"name path=q", color = mycolor, opacity = 0.25, ggplot2_plot_theme...},
            Coordinates(xs, ys[1] .- error),
        ),
        Plot(
            {thick, color = mycolor, fill = mycolor, opacity = 0.25},
            raw"fill between [of=f and g]",
        ),
        Plot(
            {thick, color = mycolor, fill = mycolor, opacity = 0.25},
            raw"fill between [of=f and q]",
        ),
    )
end


function plotBm(xs, ys, error, title)
    ggplot2_axis_theme = @pgf {
        tick_align = "outside",
        tick_pos = "left",
        xmajorgrids,
        x_grid_style = "white",
        ymajorgrids,
        y_grid_style = "white",
        axis_line_style = "white",
        "axis_background/.style" = {fill = "white!89.803921568627459!black"},
    }

    ggplot2_plot_theme = @pgf {
        #         mark="square",
        #         mark_size = 1,
        "no marks",
        "error bars/y dir=both",
        "error bars/y explicit",
        line_width = "1.64pt",
    }

    n = length(ys)

    # Evenly spread out colors
    colors = [LCHuv(65, 100, h) for h in range(15; stop = 360 + 15, length = n + 1)][1:n]
    @pgf Axis(
        {
            ggplot2_axis_theme...,
            xmin = minimum(xs),
            xmax = 1.05 * maximum(xs),
            ymin = 0,
            ymax = 1.05 * maximum(maximum(ys)),
            title = title,
            xlabel = "Time",
            ylabel = L"$B_{md}$",
        },
        [
            PlotInc(
                {ggplot2_plot_theme..., color = colors[i]},
                Coordinates(xs, _y; yerror = error),
            ) for (i, _y) in enumerate(ys)
        ]...,
    )
end

function plotBmNoVar(xs, ys, options::Dict{String,Float64})
    ggplot2_axis_theme = @pgf {
        tick_align = "outside",
        tick_pos = "left",
        xmajorgrids,
        x_grid_style = "white",
        ymajorgrids,
        y_grid_style = "white",
        axis_line_style = "white",
        "axis_background/.style" = {fill = "white!89.803921568627459!black"},
    }

    ggplot2_plot_theme = @pgf {mark = "square", mark_size = 1, line_width = "1.64pt"}

    n = length(ys)

    # Evenly spread out colors
    colors = [LCHuv(65, 100, h) for h in range(15; stop = 360 + 15, length = n + 1)][1:n]
    @pgf Axis(
        {
            ggplot2_axis_theme...,
            xmin = minimum(xs),
            xmax = 1.05 * maximum(xs),
            ymin = 0,
            ymax = 1.05 * maximum(maximum(ys)),
            title = options["title"],
            xlabel = "Time",
            ylabel = L"$B_m$",
        },
        [
            PlotInc({ggplot2_plot_theme..., color = colors[i]}, Coordinates(xs, _y))
            for (i, _y) in enumerate(ys)
        ]...,
    )
end

function plotPdetFreq(xs, ys, options::Dict{String,Float64})
    ggplot2_axis_theme = @pgf {
        tick_align = "outside",
        tick_pos = "left",
        xmajorgrids,
        x_grid_style = "white",
        ymajorgrids,
        y_grid_style = "white",
        axis_line_style = "white",
        "axis_background/.style" = {fill = "white!89.803921568627459!black"},
    }

    ggplot2_plot_theme =
        @pgf {mark = "square", mark_size = 1, mark_options = "solid", line_width = "1.64pt"}

    n = length(ys)

    # Evenly spread out colors
    colors = [LCHuv(65, 100, h) for h in range(15; stop = 360 + 15, length = n + 1)][1:n]
    @pgf SemiLogXAxis(
        {
            ggplot2_axis_theme...,
            xmin = minimum(xs),
            xmax = 1.05 * maximum(xs),
            ymin = 0,
            ymax = 1.05 * maximum(maximum(ys)),
            title = L"Neutral Growth $(u = $"*
                   "$(options["u"]),"*
                    L"$\rho = $"*
                    "$(options["rho"]))",
            xlabel = "Frequency",
            ylabel = L"$\phi(c_j,f,\rho)$",
        },
        [
            PlotInc({ggplot2_plot_theme..., color = colors[i]}, Coordinates(xs, _y))
            for (i, _y) in enumerate(ys)
        ]...,
    )
end

function plotQuick(xs, ys, title)
    ggplot2_axis_theme = @pgf {
        tick_align = "outside",
        tick_pos = "left",
        xmajorgrids,
        x_grid_style = "white",
        ymajorgrids,
        y_grid_style = "white",
        axis_line_style = "white",
        "axis_background/.style" = {fill = "white!89.803921568627459!black"},
    }

    ggplot2_plot_theme =
        @pgf {mark = "square", mark_size = 0.5, mark_options = "solid", line_width = "1pt"}

    n = length(ys)

    # Evenly spread out colors
    colors = [LCHuv(65, 100, h) for h in range(15; stop = 360 + 15, length = n + 1)][1:n]
    @pgf Axis(
        {
            ggplot2_axis_theme...,
            xmin = minimum(xs),
            xmax = 1.05 * maximum(xs),
            ymin = minimum(minimum(ys)),
            ymax = 1,
            title = title,
            xlabel = "Seeding time",
            ylabel = L"$\Pr[\mathrm{sampling\ from\ type}_i](t)$",
            legend_pos = "north west",
        },
        [
            PlotInc({ggplot2_plot_theme..., color = colors[i]}, Coordinates(xs, _y))
            for (i, _y) in enumerate(ys)
        ]...,
        Legend(["type-0", "type-1", "type-2"]),
    )
end


function plotDecays(xs, ys, title)
    ggplot2_axis_theme = @pgf {
        # tick_align = "outside",
        # tick_pos = "left",
        # xmajorgrids,
        # x_grid_style = "white",
        # ymajorgrids,
        # y_grid_style = "white",
        # axis_line_style = "white",
        # "axis_background/.style" = {fill = "white!89.803921568627459!black"},
    }

    ggplot2_plot_theme =
        @pgf {mark = "square", mark_size = 0.5, mark_options = "solid", line_width = "1pt"}

    n = length(ys)

    # Evenly spread out colors
    colors = range(HSL(colorant"red"), stop=HSL(colorant"green"), length=n)
    @pgf Axis(
        {
            ggplot2_axis_theme...,
            xmin = minimum(xs),
            xmax = 1.05 * maximum(xs),
            ymin = minimum(minimum(ys)),
            ymax = 1,
            title = title,
            xlabel = L"$j$",
            ylabel = L"$d_{c_j}$",
            # legend_pos = " east",
        },
        [
            PlotInc({ggplot2_plot_theme..., color = colors[i]}, Coordinates(xs, _y))
            for (i, _y) in enumerate(ys)
        ]...,
        Legend(["Exponential Decay", "Logarithmic Decay", "Polynomial Decay"]),
    )
end

function plotBmPMFs(xs, ys, title)
    ggplot2_axis_theme = @pgf {
        # tick_align = "outside",
        # tick_pos = "left",
        # xmajorgrids,
        # x_grid_style = "white",
        # ymajorgrids,
        # y_grid_style = "white",
        # axis_line_style = "white",
        # "axis_background/.style" = {fill = "white!89.803921568627459!black"},
    }

    ggplot2_plot_theme =
        @pgf {mark = "square", mark_size = 0.5, mark_options = "solid", line_width = "1pt"}

    n = length(ys)

    push!(PGFPlotsX.CUSTOM_PREAMBLE, raw"\usepgfplotslibrary{fillbetween}")
    push!(PGFPlotsX.CUSTOM_PREAMBLE, raw"\usepackage{cmbright}")
    push!(PGFPlotsX.CUSTOM_PREAMBLE, raw"\usepackage[T1]{fontenc}")

    # Evenly spread out colors
    colors = range(HSL(colorant"red"), stop=HSL(colorant"green"), length=n)
    @pgf Axis(
        {
            ggplot2_axis_theme...,
            xmin = minimum(xs),
            xmax = 1.05 * maximum(xs),
            ymin = minimum(minimum(ys)),
            ymax = 1,
            title = title,
            xlabel = L"$j$",
            ylabel = L"$B_{md}$",
            legend_pos = "north east",
        },
        [
            PlotInc({ggplot2_plot_theme..., color = colors[i]}, Coordinates(xs, _y))
            for (i, _y) in enumerate(ys)
        ]...,
        Legend(["Exponential Decay", "Logarithmic Decay", "Polynomial Decay"]),
    )
end



function BmPopPlotPrepare(detectionFrequency, detectionPopulation, mutationProbability, ρ)
    k_max = kmaxEstimate(detectionPopulation, mutationProbability, ρ)
    y = [Bm(detectionFrequency, mutationProbability, ρ, k_seed) for k_seed = 1:k_max]
    ys = [y]
    xs = [popEstimate(s, mutationProbability, ρ) for s = 1:k_max]
    error = sqrt.(BmVar.(detectionFrequency, mutationProbability, ρ, 1:k_max))
    options = Dict("u" => mutationProbability, "rho" => ρ, "f" => detectionFrequency)
    return xs, ys, error, options
end
function BmMutIndexPlotPrepare(detectionFrequency, MAX_MUT, mutationProbability, ρ)
    k_max = MAX_MUT
    y = [Bm(detectionFrequency, mutationProbability, ρ, k_seed) for k_seed = 1:k_max]
    ys = [y]
    xs = 1:k_max
    error = sqrt.(BmVar.(detectionFrequency, mutationProbability, ρ, 1:k_max))
    options = Dict("u" => mutationProbability, "rho" => ρ, "f" => detectionFrequency)
    return xs, ys, error, options
end
function BpMutIndexPlotPrepare(MAX_MUT, Bm_s, McMs_sum)
    k_max = MAX_MUT
#     y = [max(McMs_sum - k_seed + Bm_s[k_seed], 0) for k_seed = 1:k_max]
    y = [McMs_sum - k_seed + Bm_s[k_seed] for k_seed = 1:k_max]
    ys = [y]
    xs = 1:k_max
    return xs, ys
end


function plotBpUpperLower(xs, ys, title)
    ggplot2_axis_theme = @pgf {
        tick_align = "outside",
        tick_pos = "left",
        xmajorgrids,
        x_grid_style = "white",
        ymajorgrids,
        y_grid_style = "white",
        axis_line_style = "white",
        "axis_background/.style" = {fill = "white!89.803921568627459!black"},
    }

    ggplot2_plot_theme = @pgf {
        #         mark="square",
        #         mark_size = 1,
        "no marks",
        "error bars/y dir=both",
        "error bars/y explicit",
        line_width = "1.64pt",
    }

    n = length(ys)

    # Evenly spread out colors
    colors = [LCHuv(65, 100, h) for h in range(15; stop = 360 + 15, length = n + 1)][1:n]
    @pgf SemiLogXAxis(
        {
            ggplot2_axis_theme...,
            xmin = minimum(xs),
            xmax = 1.05 * maximum(xs),
            #              ymin = minimum(ys),   ymax = 1.05*maximum(maximum(ys)),
            title = title,
            xlabel = "Seeding Time",
            ylabel = L"$B_p$",
        },
        [
            PlotInc({ggplot2_plot_theme..., color = colors[i]}, Coordinates(xs, _y))
            for (i, _y) in enumerate(ys)
        ]...,
    )
end



function plotBp(xs, ys, error, title)
    ggplot2_axis_theme = @pgf {
        tick_align = "outside",
        tick_pos = "left",
        # xmajorgrids,
        # x_grid_style = "white",
        # ymajorgrids,
        # y_grid_style = "white",
        axis_line_style = "white",
        "axis_background/.style" = {fill = "white!89.803921568627459!black"},
    }

    ggplot2_plot_theme = @pgf {
        #         mark="square",
        #         mark_size = 1,
        "no marks",
        "error bars/y dir=both",
        "error bars/y explicit",
        line_width = "1.64pt",
    }

    n = length(ys)
    n = 3

    # Evenly spread out colors
    colors = [LCHuv(65, 100, h) for h in range(15; stop = 360 + 15, length = n + 1)][1:n]
    push!(PGFPlotsX.CUSTOM_PREAMBLE, raw"\usepgfplotslibrary{fillbetween}")
    push!(PGFPlotsX.CUSTOM_PREAMBLE, raw"\usepackage{cmbright}")
    push!(PGFPlotsX.CUSTOM_PREAMBLE, raw"\usepackage[T1]{fontenc}")
    mycolor = colors[3]
    @pgf SemiLogXAxis(
        #@pgf Axis(
        {
            ggplot2_axis_theme...,
            xmin = minimum(xs),
            xmax = 1.05 * maximum(xs),
            ymin = 0,
            ymax = 1.05 * maximum(maximum(ys[1] .+ error)),
            title = title,
            xlabel = "Seeding Time", #"Cell Population",
            ylabel = L"$B_{p}$",
        },
        # [
        #     PlotInc(
        #         {ggplot2_plot_theme..., color = colors[i]},
        #         Coordinates(xs, _y; yerror = error),
        #     ) for (i, _y) in enumerate(ys)
        # ]...,
        Plot(
            {"name path=f", color = mycolor, ggplot2_plot_theme...},
            Coordinates(xs, ys[1]),
        ),
        Plot(
            {"name path=g", color = mycolor, opacity = 0.25, ggplot2_plot_theme...},
            Coordinates(xs, ys[1] .+ error),
        ),
        Plot(
            {"name path=q", color = mycolor, opacity = 0.25, ggplot2_plot_theme...},
            Coordinates(xs, ys[1] .- error),
        ),
        Plot(
            {thick, color = mycolor, fill = mycolor, opacity = 0.25},
            raw"fill between [of=f and g]",
        ),
        Plot(
            {thick, color = mycolor, fill = mycolor, opacity = 0.25},
            raw"fill between [of=f and q]",
        ),
    )
end
