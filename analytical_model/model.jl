using LinearAlgebra, SparseArrays, Statistics, StatsBase


# This expresses the probability the frequency of cells that carry variant a certain variant c_j is at least f. 
φ(f,u,ρ,c) = (u / (u-log(1+(ρ-1)*f)))^c

# This calculates the expected number of cells in the population at a given time (measured here based on the discussion 
# of the paper (see Equations 9 and 10 for details)
N(s,ρ,u) = exp((ρ/u -1/u)*s) + 2*sinh((1/u - ρ/u)*s)/ρ


function gmap(seed,u,ρ, debug = false)
    # Calculation of the approximate mapping g(c_j) 
    # For details see Equation (10) of the manuscript
    if debug
        println("exp part = ",exp((ρ/u - 1/u)*seed))
        println("sinh part = ",2*sinh((1/u - ρ/u)*seed)/ρ)
    end
    return seed == 0 ? 0 : ceil(u*( exp((ρ/u - 1/u)*seed) + 2*sinh((1/u - ρ/u)*seed)/ρ))
end

function popEstimate(seed,u,ρ, debug = false)
    # This function returns an estimate of the expected number of cells in the population, 
    # at the time  each of the seeding-cell-specific variants appeared. For details see the manuscript (Equations 9 and 10)
    if debug
        println("exp part = ",exp((ρ/u - 1/u)*seed))
        println("sinh part = ",2*sinh((1/u - ρ/u)*seed)/ρ)
    end
    return exp((ρ/u - 1/u)*seed) + 2*sinh((1/u - ρ/u)*seed)/ρ
end


function calculate_sampling_probabilities_populations(u, ts1, ts2, t_det, ρ, ρ1, ρ2 )
    # This refers to the 3rd population based seeding strategy
    # For context see the related discussion in Section 4.1.3 of the paper.
    p0 = zeros(t_det)
    p1 = zeros(t_det)
    p2 = zeros(t_det)
    for k=1:t_det
        if k < ts1
            N0 = N(k,ρ,u)
            p0[k] = 1
            continue
        end
        if k < ts2
            N0 = N(k,ρ,u)
            N1 = N(k-ts1,ρ1,u)
            p0[k] = N0 / (N0 + N1)
            p1[k] = N1 / (N0 + N1)
            continue
        end
        N0 = N(k,ρ,u)
        N1 = N(k-ts1,ρ1,u)
        N2 = N(k-ts2,ρ2,u)
        p0[k] = (N0) / (N0+N1+N2)
        p1[k] = (N1) / (N0+N1+N2)
        p2[k] = (N2) / (N0+N1+N2)
    end
    p0, p1, p2
end


function Bm(f,u,ρ,seed, debug = false)
    Bm = 0
    for j=0:seed
        Ptemp = 1.0
        if j != seed 
            p = (1-φ(f,u,ρ,gmap(j+1,u,ρ))) - (1-φ(f,u,ρ,gmap(j,u,ρ)))
        else
            p = φ(f,u,ρ,gmap(j,u,ρ))
        end
        Bm += (seed-j)*p
        if debug
            println("P(U_j+1) = ", (1-φ(f,u,ρ,gmap(j+1,u,ρ))))
            println("P(U_j) = ", (1-φ(f,u,ρ,gmap(j,u,ρ))))
            println("Bm = ", Bm)
        end
    end
    return Bm
end

function BmGeneric(d, k_max, debug = false)
    Bm = 0
    BmPMF = zeros(k_max+1)
    for j=0:k_max
        if j != k_max
            p = (1-d[j+2]) - (1-d[j+1])
        else
            p = d[k_max]
        end
        Bm += (k_max-j)*p
        BmPMF[k_max-j+1] = p
    end
    return Bm, BmPMF
end


function BmDist(f,u,ρ,seed, debug = false)
    # Calculates the B_md distribution (See equation 6, and the related discussion in the paper)
    Bm = []
    for j=0:seed
        Ptemp = 1.0
        if j != seed 
            p = (1-φ(f,u,ρ,gmap(j+1,u,ρ))) - (1-φ(f,u,ρ,gmap(j,u,ρ)))
        else
            p = φ(f,u,ρ,gmap(j,u,ρ))
        end
        append!(Bm,p)
        if debug
            println("P(U_j+1) = ", (1-φ(f,u,ρ,gmap(j+1,u,ρ))))
            println("P(U_j) = ", (1-φ(f,u,ρ,gmap(j,u,ρ))))
            println("Bm = ", Bm)
        end
    end
    return Bm
end

function BmVar(f,u,ρ,seed)
    # To measure the innate variance of the random variable Bmd that follows the above distribution BmDist
    Bm = 0
    Bm2 = 0
    for j=0:seed
        Ptemp = 1.0
        if j != seed 
            p = (1-φ(f,u,ρ,gmap(j+1,u,ρ))) - (1-φ(f,u,ρ,gmap(j,u,ρ)))
        else
            p = φ(f,u,ρ,gmap(j,u,ρ))
        end
        Bm += (seed-j)*p
        Bm2 += (seed-j)^2*p
    end
    return Bm2 - Bm^2
end



# The following are simple Helper functions necessary to facilitate the production of the plots

function kmaxEstimate(population,u,ρ)
    k = 1
    while exp((ρ/u - 1/u)*k) + 2*sinh((1/u - ρ/u)*k)/ρ < population
        k += 1
    end
    return k
end

function findtd(u, ρ, population)
    s = 1
    while exp((ρ / u - 1 / u) * s) + 2 * sinh((1 / u - ρ / u) * s) / ρ < population
        s += 1
    end
    return s
end


function finds_i(u, ρ, time, population)
    x = 0.0
    ρ1 = ρ / (x + 1)
    while exp((ρ1 / u - 1 / u) * time) + 2 * sinh((1 / u - ρ1 / u) * time) / ρ1 < population
        x += 0.0005
        ρ1 = ρ / (x + 1)
    end
    return x
end

function findts_i(u, ρ, td0, population)
    t = 0
    while exp((ρ / u - 1 / u) * t) + 2 * sinh((1 / u - ρ / u) * t) / ρ < population
        t += 1
    end
    return td0 - t # returns the type origination time measured in [0..td]  (integer value, i.e., the number of mutations...)
end





