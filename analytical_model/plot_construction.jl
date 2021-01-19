include("plot_utils.jl")

function BmFig(ρ, u, f, ff_1, ff_2)
    DetectedPopulation = 10^9
    ft0 = 1 - ff_1
    ft1 = ff_1 - ff_2
    ft2 = ff_2
    if ft0 <= 0
        println("Wrong selected frequencies. Exiting....")
        return
    end

    N(s, ρ, u) = exp((ρ / u - 1 / u) * s) + 2 * sinh((1 / u - ρ / u) * s) / ρ

    td0 = findtd(u, ρ, DetectedPopulation * ft0) # time till the number of exclusively type-0 cells are at least DetectedPopulation*ft0
    ts1 = Int.(floor(td0 / 3)) # time of appearance of type-1 cells
    ts2 = Int.(ceil(2 * td0 / 3)) # time of appearance of type-2 cells
    s_1 = finds_i(u, ρ, td0 - ts1, DetectedPopulation * ft1) # calculates the selective growth advantage necessary for the new type to reach
    # at least DetectedPopulation*ft1 number of cells after td0 - ts1 time.
    # this is a helper function to allow us to select the values that make the pattern visually apparent for the example.
    ρ1 = ρ / (s_1 + 1)
    s_2 = finds_i(u, ρ1, td0 - ts2, DetectedPopulation * ft2)
    ρ2 = ρ1 / (s_2 + 1)
#     @show td0, ts1, s_1, ts2, s_2, ρ1, ρ2
    overallPopulation = (N(td0, ρ, u) + N(td0 - ts1, ρ1, u) + N(td0 - ts2, ρ2, u))
    type1Ratio = N(td0 - ts1, ρ1, u) / overallPopulation
    type2Ratio = N(td0 - ts2, ρ2, u) / overallPopulation

    println("========================================SUMMARY=========================================== ")
    println(
        "The growth parameters of the example are:\n\t ρ = ",
        ρ,
        "\n\t mutation probability: u = ",
        u,
        "\n\t detectability frequency: f = ",
        f,
        "\n\t rho of type 1: ρ_1 = ",
        ρ1,
        "\n\t rho of type 2: ρ_2 = ",
        ρ2,
        "",
    )


    println("The population at detection time is ", overallPopulation)
    println(
        "The type-1 mutation appears in ",
        100 * (type1Ratio + type2Ratio),
        "% of the population",
    )
    println("The type-2 mutation appears in ", 100 * type2Ratio, "% of the population")
    println("========================================================================================== \n\n")

    # Create the 0-->ts1 part of Bm
    xs, ys, error, options = BmMutIndexPlotPrepare(f / ft0, td0, u, ρ)

    # Create the ts1-->ts2 part of Bm
    xs1, ys1, error1, options1 = BmMutIndexPlotPrepare(f / ft1, td0 - ts1, u, ρ1)

    # Create the ts1-->detection part of the plot
    xs2, ys2, error2, options2 = BmMutIndexPlotPrepare(f / ft2, td0 - ts2, u, ρ2)

    # @show length(ys[1]), length(ys1[1]), length(ys2[1]), td0-ts2

    xs_pop = [popEstimate(s, u, ρ) for s = 1:td0]
    [xs_pop[s] += popEstimate(s - ts1, u, ρ1) for s = ts1:ts2]
    [xs_pop[s] += popEstimate(s - ts2, u, ρ2) for s = ts2:td0]


    # Probabilities are a function of the population

    # Construct final ys
    _y0 = deepcopy(ys[1])
    _y1 = deepcopy(ys1[1])
    _y2 = deepcopy(ys2[1])
    # after the first type is introduced
    #p0, p1, p2 = calculate_sampling_probabilities(u, ts1, ts2, td0, ρ, ρ1, ρ2 )
    p0, p1, p2 = calculate_sampling_probabilities_populations(u, ts1, ts2, td0, ρ, ρ1, ρ2)
    h_timeseries = plotQuick(1:length(p0), [p0, p1, p2], "Sampling Probabilities over Time")
    [_y0[s+ts1] = p1[s+ts1] * _y1[s] + p0[s+ts1] * _y0[s+ts1] for s = 1:(td0-ts1)]
    # after the second type is introduced
    [
        _y0[ts2+s] =
            p0[ts2+s] * _y0[ts2+s] + p1[ts2+s] * _y1[(ts2-ts1)+s] + p2[ts2+s] * _y2[s]
        for s = 1:td0-ts2
    ]
    y_out = [_y0]
    _m = deepcopy(_y0)
    
    ## Seeding randomly from the living population  (3rd strategy. See paper for details)
    title = "Seeding randomly from the living population"
    error = zeros(length(_y0))
    [
        error[s+ts1] = sqrt(
            (_m[s+ts1] - _y0[s+ts1])^2 * p0[s+ts1] + (_m[s+ts1] - _y1[s])^2 * p1[s+ts1],
        ) for s = 1:(td0-ts1)
    ]
    [
        error[s+ts2] = sqrt(
            (_m[s+ts2] - _y0[s+ts2])^2 * p0[s+ts2] +
            (_m[s+ts2] - _y1[s+ts2-ts1])^2 * p1[s+ts2] +
            (_m[s+ts2] - _y2[s])^2 * p2[s+ts2],
        ) for s = 1:td0-ts2
    ]
    axs1 = plotBmPop(xs_pop, y_out, error, title)

    ## Only New type can Seed
    title = "Seeding from the most advanced clone"
    # Construct final ys
    _y0 = deepcopy(ys[1])
    _y1 = deepcopy(ys1[1])
    _y2 = deepcopy(ys2[1])
    # after the first type is introduced
    p1 = 1
    p0 = 0
    [_y0[s+ts1] = p1 * _y1[s] + p0 * _y0[s+ts1] for s = 1:td0-ts1]
    # after the second type is introduced
    p0 = 0
    p1 = 0
    p2 = 1
    [_y0[ts2+s] = p0 * _y0[ts2+s] + p1 * _y1[(ts2-ts1)+s] + p2 * _y2[s] for s = 1:(td0-ts2)]
    y_out = [_y0]
    _m = deepcopy(_y0)
    error = zeros(length(_y0))
    [
        error[s+ts1] = sqrt((_m[s+ts1] - _y0[s+ts1])^2 * p0 + (_m[s+ts1] - _y1[s])^2 * p1)
        for s = 1:ts2-ts1
    ]
    [
        error[s+ts2] = sqrt(
            (_m[s+ts2] - _y0[s+ts2])^2 * p0 +
            (_m[s+ts2] - _y1[s+ts2-ts1])^2 * p1 +
            (_m[s+ts2] - _y2[s])^2 * p2,
        ) for s = 1:td0-ts2
    ]
    axs2 = plotBmPop(xs_pop, y_out, error, title)
    # @show _m[ts2+1:end] - _y2

    ## Probabilities Proportional to fitness
    title = "Seeding proportional to subclone fitness"
    # Construct final ys
    _y0 = deepcopy(ys[1])
    _y1 = deepcopy(ys1[1])
    _y2 = deepcopy(ys2[1])
    # after the first type is introduced
    p1 = (1 + s_1) / (2 + s_1)
    p0 = (1) / (2 + s_1)
    [_y0[s+ts1] = p1 * _y1[s] + p0 * _y0[s+ts1] for s = 1:(td0-ts1)]
    # after the second type is introduced
    p0 = 1 / (2 + s_1 + (1 + s_1) * (1 + s_2))
    p1 = (1 + s_1) / (2 + s_1 + (1 + s_1) * (1 + s_2))
    p2 = (1 + s_1) * (1 + s_2) / (2 + s_1 + (1 + s_1) * (1 + s_2))
    [
        _y0[ts2+s] = p0 * _y0[ts2+s] + p1 * _y1[(ts2-ts1)+s] + p2 * _y2[s]
        for s = 1:length(_y2)
    ]
    y_out = [_y0]
    _m = deepcopy(_y0)
    error = zeros(length(_y0))
    [
        error[s+ts1] = sqrt((_m[s+ts1] - _y0[s+ts1])^2 * p0 + (_m[s+ts1] - _y1[s])^2 * p1)
        for s = 1:(td0-ts1)
    ]
    [
        error[s+ts2] = sqrt(
            (_m[s+ts2] - _y0[s+ts2])^2 * p0 +
            (_m[s+ts2] - _y1[s+ts2-ts1])^2 * p1 +
            (_m[s+ts2] - _y2[s])^2 * p2,
        ) for s = 1:length(_y2)
    ]
    axs3 = plotBmPop(xs_pop, y_out, error, title)

    h = @pgf GroupPlot(
        {
            group_style =
                {group_size = "3 by 1", horizontal_sep = "2cm", vertical_sep = "3cm"},
        },
        axs2,
        axs3,
        axs1,
    )

    return h, h_timeseries, axs1, axs2, axs3
end

# Z is the total number of mutations comprising Bp at detection time. It is assumed to be a measurable quantity 
function BpFig(Z, ρ, u, f, ff1, ff2)
    DetectedPopulation = 10^9
    ft0 = 1 - ff_1
    ft1 = ff_1 - ff_2
    ft2 = ff_2
    if ft0 <= 0
        println("Wrong selected frequencies. Exiting....")
        return
    end
    f_bp = f
    f_bm = f
    N(s, ρ, u) = exp((ρ / u - 1 / u) * s) + 2 * sinh((1 / u - ρ / u) * s) / ρ

    td0 = findtd(u, ρ, DetectedPopulation * ft0) # time till the number of exclusively type-0 cells 
                                                 # are at least DetectedPopulation*ft0
    ts1 = Int.(floor(td0 / 3)) # time of appearance of type-1 cells
    ts2 = Int.(ceil(2 * td0 / 3)) # time of appearance of type-2 cells
    s_1 = finds_i(u, ρ, td0 - ts1, DetectedPopulation * ft1) # calculates the selective growth advantage necessary for the new type to reach
    # at least DetectedPopulation*ft1 number of cells after td0 - ts1 time.
    # this is a helper function to allow us to select the values that make the pattern visually apparent for the example.
    ρ1 = ρ / (s_1 + 1)
    s_2 = finds_i(u, ρ1, td0 - ts2, DetectedPopulation * ft2)
    ρ2 = ρ1 / (s_2 + 1)
    # @show ft0, ft1, ft2
    overallPopulation = (N(td0, ρ, u) + N(td0 - ts1, ρ1, u) + N(td0 - ts2, ρ2, u))
    type1Ratio = N(td0 - ts1, ρ1, u) / overallPopulation
    type2Ratio = N(td0 - ts2, ρ2, u) / overallPopulation


    # Create the 0-->ts1 part of Bm
#     @show f, ft0 ,ft1, ft2
    xs, ys, error, options = BmMutIndexPlotPrepare(f / ft0, td0, u, ρ)

    # Create the ts1-->ts2 part of Bm
    xs1, ys1, error1, options1 = BmMutIndexPlotPrepare(f / ft1, td0 - ts1, u, ρ1)

    # Create the ts1-->detection part of the plot
    xs2, ys2, error2, options2 = BmMutIndexPlotPrepare(f / ft2, td0 - ts2, u, ρ2)

    # @show length(ys[1]), length(ys1[1]), length(ys2[1]), td0-ts2

    xs_pop = [popEstimate(s, u, ρ) for s = 1:td0]
    [xs_pop[s] += popEstimate(s - ts1, u, ρ1) for s = ts1:ts2]
    [xs_pop[s] += popEstimate(s - ts2, u, ρ2) for s = ts2:td0]


    ## Seeding randomly from the population
    title = "Seeding randomly from the population"
    # Construct final ys
    _y0 = deepcopy(ys[1])
    _y1 = deepcopy(ys1[1])
    _y2 = deepcopy(ys2[1])
    # after the first type is introduced
    p0, p1, p2 = calculate_sampling_probabilities_populations(u, ts1, ts2, td0, ρ, ρ1, ρ2)


    h_timeseries = plotQuick(1:length(p0), [p0, p1, p2], "Sampling Probabilities over Time")
    [_y0[s+ts1] = p1[s+ts1] * _y1[s] + p0[s+ts1] * _y0[s+ts1] for s = 1:(td0-ts1)]
    # after the second type is introduced
    [
        _y0[ts2+s] =
            p0[ts2+s] * _y0[ts2+s] + p1[ts2+s] * _y1[(ts2-ts1)+s] + p2[ts2+s] * _y2[s]
        for s = 1:td0-ts2
    ]
    y_out = [_y0]
    _m = deepcopy(_y0)
    error = zeros(length(_y0))
    [
        error[s+ts1] = sqrt(
            (_m[s+ts1] - _y0[s+ts1])^2 * p0[s+ts1] + (_m[s+ts1] - _y1[s])^2 * p1[s+ts1],
        ) for s = 1:(td0-ts1)
    ]
    [
        error[s+ts2] = sqrt(
            (_m[s+ts2] - _y0[s+ts2])^2 * p0[s+ts2] +
            (_m[s+ts2] - _y1[s+ts2-ts1])^2 * p1[s+ts2] +
            (_m[s+ts2] - _y2[s])^2 * p2[s+ts2],
        ) for s = 1:td0-ts2
    ]
    axs1 = plotBmPop(xs_pop, y_out, error, title)

    # Calculate Bp
    bpys = [
        [_m[k_seed] - k_seed + Z - 1 for k_seed = 1:td0]
    ]
    Bpaxs1 = plotBp(xs_pop, bpys, error, title)

    ## Only New type can Seed
    title = "Seeding from the most advanced clone"
    # Construct final ys
    _y0 = deepcopy(ys[1])
    _y1 = deepcopy(ys1[1])
    _y2 = deepcopy(ys2[1])
    # after the first type is introduced
    p1 = 1
    p0 = 0
    [_y0[s+ts1] = p1 * _y1[s] + p0 * _y0[s+ts1] for s = 1:td0-ts1]
    # after the second type is introduced
    p0 = 0
    p1 = 0
    p2 = 1
    [_y0[ts2+s] = p0 * _y0[ts2+s] + p1 * _y1[(ts2-ts1)+s] + p2 * _y2[s] for s = 1:(td0-ts2)]
    y_out = [_y0]
    _m = deepcopy(_y0)
    error = zeros(length(_y0))
    [
        error[s+ts1] = sqrt((_m[s+ts1] - _y0[s+ts1])^2 * p0 + (_m[s+ts1] - _y1[s])^2 * p1)
        for s = 1:ts2-ts1
    ]
    [
        error[s+ts2] = sqrt(
            (_m[s+ts2] - _y0[s+ts2])^2 * p0 +
            (_m[s+ts2] - _y1[s+ts2-ts1])^2 * p1 +
            (_m[s+ts2] - _y2[s])^2 * p2,
        ) for s = 1:td0-ts2
    ]
    axs2 = plotBmPop(xs_pop, y_out, error, title)


    # Calculate Bp
    bpys = [
        [_m[k_seed] - k_seed + Z - 1 for k_seed = 1:td0]
    ]
    Bpaxs2 = plotBp(xs_pop, bpys, error, title)


    # Probabilities Proportional to fitness
    title = "Seeding proportional to subclone fitness"
    # Construct final ys
    _y0 = deepcopy(ys[1])
    _y1 = deepcopy(ys1[1])
    _y2 = deepcopy(ys2[1])
    # after the first type is introduced
    p1 = (1 + s_1) / (2 + s_1)
    p0 = (1) / (2 + s_1)
    [_y0[s+ts1] = p1 * _y1[s] + p0 * _y0[s+ts1] for s = 1:(td0-ts1)]
    # after the second type is introduced
    p0 = 1 / (2 + s_1 + (1 + s_1) * (1 + s_2))
    p1 = (1 + s_1) / (2 + s_1 + (1 + s_1) * (1 + s_2))
    p2 = (1 + s_1) * (1 + s_2) / (2 + s_1 + (1 + s_1) * (1 + s_2))
    [
        _y0[ts2+s] = p0 * _y0[ts2+s] + p1 * _y1[(ts2-ts1)+s] + p2 * _y2[s]
        for s = 1:length(_y2)
    ]
    y_out = [_y0]
    _m = deepcopy(_y0)
    
    error = zeros(length(_y0))
    [
        error[s+ts1] = sqrt((_m[s+ts1] - _y0[s+ts1])^2 * p0 + (_m[s+ts1] - _y1[s])^2 * p1)
        for s = 1:(td0-ts1)
    ]
    [
        error[s+ts2] = sqrt(
            (_m[s+ts2] - _y0[s+ts2])^2 * p0 +
            (_m[s+ts2] - _y1[s+ts2-ts1])^2 * p1 +
            (_m[s+ts2] - _y2[s])^2 * p2,
        ) for s = 1:length(_y2)
    ]
    axs3 = plotBmPop(xs_pop, y_out, error, title)


    # Calculate Bp
    bpys = [
        [_m[k_seed] - k_seed + Z - 1 for k_seed = 1:td0]
    ]
    Bpaxs3 = plotBp(xs_pop, bpys, error, title)

    h = @pgf GroupPlot(
        {
            group_style =
                {group_size = "3 by 2", horizontal_sep = "2cm", vertical_sep = "3cm"},
        },
        axs1,
        axs3,
        axs2,
        Bpaxs1,
        Bpaxs3,
        Bpaxs2,
    )

    return h, h_timeseries, Bpaxs1, Bpaxs2, Bpaxs3
end

function plotBmNeutral(xs, ys, title)
    ggplot2_axis_theme = @pgf {
    }

    ggplot2_plot_theme =
        @pgf {mark = "square", mark_size = 0, mark_options = "solid", line_width = "1.64pt"}

    n = length(ys)

    # Evenly spread out colors
    colors = range(HSL(colorant"red"), stop=HSL(colorant"green"), length=n)
    @pgf SemiLogXAxis(
        {
            ggplot2_axis_theme...,
            xmin = minimum(xs),
            xmax = 1.05 * maximum(xs),
            ymin = minimum(minimum(ys)),
            ymax = maximum(maximum(ys)),
            title = title,
            xlabel = "Seeding time",
            ylabel = L"$\mathbb{E}[B^{\mathit{neutral}}_{md}]$",
            legend_pos = "north west",
        },
        [
            PlotInc({ggplot2_plot_theme..., color = colors[i]}, Coordinates(xs, _y))
            for (i, _y) in enumerate(ys)
        ]...,
        Legend([L"$f = 0.005$", L"$f = 0.01$", L"$f = 0.05$", L"$f = 0.1$"]),
    )
end