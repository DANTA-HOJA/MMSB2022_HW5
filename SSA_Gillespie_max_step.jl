using StatsBase # Weights() and sample()
using Random    # randexp()
using Plots
using Interpolations
using Statistics # mean()
Random.seed!(2022)


#=
Stochastic chemical reaction: Gillespie Algorithm (direct method)
Adapted from: Chemical and Biomedical Enginnering Calculations Using Python Ch.4-3

# alpha 越大時前兩個 propensities 越大 => dt會越小 => 所以同樣的 t 要跑更多步驟，因此改成限制最大反應次數（step）
=#
function ssa_direct(model, u0::AbstractArray, max_step, p, stoich; tstart=zero(max_step))
    t = tstart   # Current time
    ts = [t]     # Time points
    u = copy(u0) # Current state
    us = copy(u) # States over time
    while length(ts) < max_step + 1
        a = model(u, p, t)               # propensities
        dt = randexp() / sum(a)          # Time step for the direct method
        du = sample(stoich, Weights(a))  # Choose the stoichiometry for the next reaction
        u .+= du  # Update state
        t += dt   # Update time

        us = [us u]  # Append state
        push!(ts, t) # Append time point
    end
    # Trasnpose to make columns as variables, rows as observations
    us = collect(us')
    return (t = ts, u = us)
end

#=
Stochastic chemical reaction: Gillespie Algorithm (first reaction method)
Adapted from: Chemical and Biomedical Enginnering Calculations Using Python Ch.4-3

# alpha 越大時前兩個 propensities 越大 => dt會越小 => 所以同樣的 t 要跑更多步驟，因此改成限制最大反應次數（step）
=#
function ssa_first(model, u0, max_step, p, stoich; tstart=zero(max_step))
    t = tstart   # Current time
    ts = [t]     # Time points
    u = copy(u0) # Current state
    us = copy(u) # States over time
    while length(ts) < max_step + 1
        a = model(u, p, t)               # propensities
        dts = randexp(length(a)) ./ a    # dts from all reactions
        # Choose the reaction 
        i = argmin(dts)
        # print(i)
        dt = dts[i]
        du = stoich[i]
        # Update state and time
        u .+= du
        t += dt
        us = [us u]  # Append state variable to record
        push!(ts, t) # Append time point to record
    end
    # Make column as variables, rows as observations
    us = collect(us')
    return (t = ts, u = us)
end

"""
Propensity model for this reaction.
Reaction of A <-> B with rate constants k1 & k2
"""

#= parameters：

Stoichiometry：
        for Reaction_1(R1), s1 = [1, 0]
        for Reaction_2(R2), s2 = [0, 1]
        for Reaction_3(R3), s3 = [-1, 0]
        for Reaction_4(R4), s4 = [0, -1]


3 sets of initial conditions (u0)：
    (p1, p2) = (100, 50),
               (50, 100),
               (75, 75)

alpha(α) = 5, 50, 500, and 5000

beta(β) = 4, delta(δ) = 1

N1 and N2 stand for the amount of P1 and P2, respectively.

=#

# Reaction_1, Reaction_2 實際是 alpha * (1^b / (1^b + Nx^b)，alpha是最大生成速率，若沒有beta (即beta =1)應該就是效果比較弱
model(u, p, t) = [p.alpha / (1+u[2]^p.beta),
                  p.alpha / (1+u[1]^p.beta),
                  p.delta * u[1],
                  p.delta * u[2]]

# model(u, p, t) = [p.alpha / (1+u[2]),
#                   p.alpha / (1+u[1]),
#                   p.delta * u[1],
#                   p.delta * u[2]]

# change alpha(α) => 5, 50, 500, and 5000
# change initial conditions => (p1, p2) = (100, 50), (50, 100), and (75, 75)
parameters = (alpha=5.0, beta=4.0, delta=1.0, stoich=[[1, 0], [0, 1], [-1, 0], [0, -1]])
max_step = 110000.0

# alpha 越大時前兩個 propensities 越大 => dt會越小 => 所以同樣的 t 要跑更多步驟，因此改成限制最大反應次數（step）

# create figure
fig_ssa_direct = plot(xlabel="time", ylabel="# of molecules", legend=:topleft,
                      title = "SSA (direct method), α = $(parameters.alpha)", titlefontsize=12)

fig_ssa_first = plot(xlabel="time", ylabel="# of molecules", legend=:topleft,
                     title = "SSA (1st reaction method), α = $(parameters.alpha)", titlefontsize=12)



"""
initial condition 1
"""

u0 = [100, 50] # [N1(P1), N2(P2)]

soldirect_cond1 = ssa_direct(model, u0, max_step, parameters, parameters.stoich)
solfirst_cond1 = ssa_first(model, u0, max_step, parameters, parameters.stoich)

plot!(fig_ssa_direct, soldirect_cond1.t, soldirect_cond1.u, label=["P1_init($(u0[1]))" "P2_init($(u0[2]))"])
plot!(fig_ssa_first, solfirst_cond1.t, solfirst_cond1.u, label=["P1_init($(u0[1]))" "P2_init($(u0[2]))"])

# initial condition 1 獨立畫圖
plot(xlabel="time", ylabel="# of molecules", legend=:topleft,
     title = "SSA (direct method), initial condition 1, α = $(parameters.alpha)", titlefontsize=12,
     soldirect_cond1.t, soldirect_cond1.u, label=["P1_init($(u0[1]))" "P2_init($(u0[2]))"])

plot(xlabel="time", ylabel="# of molecules", legend=:topleft,
     title = "SSA (1st reaction method), initial condition 1, α = $(parameters.alpha)", titlefontsize=12,
     solfirst_cond1.t, solfirst_cond1.u, label=["P1_init($(u0[1]))" "P2_init($(u0[2]))"])


"""
initial condition 2
"""

u0 = [50, 100] # [N1(P1), N2(P2)]

soldirect_cond2 = ssa_direct(model, u0, max_step, parameters, parameters.stoich)
solfirst_cond2 = ssa_first(model, u0, max_step, parameters, parameters.stoich)

plot!(fig_ssa_direct, soldirect_cond2.t, soldirect_cond2.u, label=["P1_init($(u0[1]))" "P2_init($(u0[2]))"])
plot!(fig_ssa_first, solfirst_cond2.t, solfirst_cond2.u, label=["P1_init($(u0[1]))" "P2_init($(u0[2]))"])

# initial condition 2 獨立畫圖
plot(xlabel="time", ylabel="# of molecules", legend=:topleft,
     title = "SSA (direct method), initial condition 2, α = $(parameters.alpha)", titlefontsize=12,
     soldirect_cond2.t, soldirect_cond2.u, label=["P1_init($(u0[1]))" "P2_init($(u0[2]))"])

plot(xlabel="time", ylabel="# of molecules", legend=:topleft,
     title = "SSA (1st reaction method), initial condition 2, α = $(parameters.alpha)", titlefontsize=12,
     solfirst_cond2.t, solfirst_cond2.u, label=["P1_init($(u0[1]))" "P2_init($(u0[2]))"])


"""
initial condition 3
"""

u0 = [75, 75] # [N1(P1), N2(P2)]

soldirect_cond3 = ssa_direct(model, u0, max_step, parameters, parameters.stoich)
solfirst_cond3 = ssa_first(model, u0, max_step, parameters, parameters.stoich)

plot!(fig_ssa_direct, soldirect_cond3.t, soldirect_cond3.u, label=["P1_init($(u0[1]))" "P2_init($(u0[2]))"])
plot!(fig_ssa_first, solfirst_cond3.t, solfirst_cond3.u, label=["P1_init($(u0[1]))" "P2_init($(u0[2]))"])

# initial condition 3 獨立畫圖
plot(xlabel="time", ylabel="# of molecules", legend=:topleft,
     title = "SSA (direct method), initial condition 3, α = $(parameters.alpha)", titlefontsize=12,
     soldirect_cond3.t, soldirect_cond3.u, label=["P1_init($(u0[1]))" "P2_init($(u0[2]))"])

plot(xlabel="time", ylabel="# of molecules", legend=:topleft,
     title = "SSA (1st reaction method), initial condition 3, α = $(parameters.alpha)", titlefontsize=12,
     solfirst_cond3.t, solfirst_cond3.u, label=["P1_init($(u0[1]))" "P2_init($(u0[2]))"])