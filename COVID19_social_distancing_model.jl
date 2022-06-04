using Agents
using Random
using InteractiveDynamics
using CairoMakie

mutable struct Agent <: AbstractAgent
    id::Int                 # Mandatory Agent identifier
    pos::NTuple{2,Float64}  # Position, required for agents in the ContinuousSpace
    vel::NTuple{2,Float64}  # Moving speeds
    mass::Float64           # Can move or not
end

function ball_model(; speed = 0.002)
    space2d = ContinuousSpace((1, 1), 0.02)
    model = ABM(Agent, space2d, properties = Dict(:dt => 1.0), rng = MersenneTwister(42))

    # Add agents to the model
    for ind in 1:500
        pos = Tuple(rand(model.rng, 2))
        vel = sincos(2π * rand(model.rng)) .* speed
        mass = 1.0
        add_agent!(pos, model, vel, mass)
    end
    return model
end

model = ball_model()

mutable struct PoorSoul <: AbstractAgent
    id::Int
    pos::NTuple{2,Float64}
    vel::NTuple{2,Float64}
    mass::Float64
    days_infected::Int  # number of days since is infected
    status::Symbol  # :S, :I or :R
    β::Float64
end

const steps_per_day = 24 # One tick per hour

function sir_initiation(;
    infection_period = 30 * steps_per_day,
    detection_time = 14 * steps_per_day,
    reinfection_probability = 0.02,
    isolated = 0.0, # *** in percentage ***
    interaction_radius = 0.010,
    dt = 1.0,
    speed = 0.001,
    death_rate = 0.044,
    N = 1000,
    initial_infected = 2,
    seed = 42,
    βmin = 0.01,
    βmax = 0.04,
)

    properties = (;
        infection_period,
        reinfection_probability,
        detection_time,
        death_rate,
        interaction_radius,
        dt,
    )
    space = ContinuousSpace((1,1), 0.02)
    model = ABM(PoorSoul, space, properties = Dict(pairs(properties)), rng = MersenneTwister(seed))

    # Add initial individual agents
    for ind in 1:N
        pos = Tuple(rand(model.rng, 2))
        status = ind ≤ N - initial_infected ? :S : :I
        isisolated = ind ≤ isolated * N
        mass = isisolated ? Inf : 1.0
        vel = isisolated ? (0.0, 0.0) : sincos(2π * rand(model.rng)) .* speed

        β = (βmax - βmin) * rand(model.rng) + βmin
        add_agent!(pos, model, vel, mass, 0, status, β)
    end

    return model
end

# sir_model = sir_initiation()

# sir_colors(a) = a.status == :S ? "#2b2b33" : a.status == :I ? "#bf2642" : "#338c54"

# fig, abmstepper = abmplot(sir_model; ac = sir_colors)
# fig # display figure

function transmit!(a1, a2, reinfection_probability)

    infected, other = a1.status == :I ? (a1, a2) : (a2, a1)

    # No one infected
    if infected.status != :I
        return
    end

    # if the other agent is already infected, do nothing
    if other.status == :I
        return
    end

    # Lucky and not infected
    if rand(model.rng) <= infected.β
        return
    end

    # Risk of reinfection
    if other.status == :R && rand(model.rng) > reinfection_probability
        return
    end

    # You got virus
    other.status = :I
end

function sir_model_step!(model)
    r = model.interaction_radius
    for (a1, a2) in interacting_pairs(model, r, :nearest)
        transmit!(a1, a2, model.reinfection_probability)
        elastic_collision!(a1, a2, :mass)
    end
end

# Agent-specific functions
function update!(agent) 
    if agent.status == :I
        agent.days_infected += 1
    end
end

function recover_or_die!(agent, model)
    if agent.days_infected ≥ model.infection_period
        if rand(model.rng) ≤ model.death_rate
            kill_agent!(agent, model)
        else
            agent.status = :R
            agent.days_infected = 0
        end
    end
end

function sir_agent_step!(agent, model)
    move_agent!(agent, model, model.dt)
    update!(agent)
    recover_or_die!(agent, model)
end



#=

SIR model Agent Color：

Susceptible  => "#2b2b33"（black）
Infectious   => "#bf2642"（red）
Recovered    => "#338c54"（green）

check color here：https://htmlcolorcodes.com/

=#
sir_model_test = sir_initiation()

sir_colors(a) = a.status == :S ? "#2b2b33" : a.status == :I ? "#bf2642" : "#338c54"

# fig, abmstepper = abmplot(sir_model_test; ac = sir_colors)

abmvideo(
    "_static/sir_model_test.mp4",
    sir_model_test,
    sir_agent_step!,
    sir_model_step!;
    title = "SIR model Video test", # title shown on output_video（it is not a file_name）.
    frames = 200, # How many frames to record in total, including the starting frame.
    ac = sir_colors, # color of agents
    as = 10, # size of agents
    am = '♠', # marker of agents
    spf = 1, # Steps-per-frame, i.e. how many times to step the model before recording a new frame.
    framerate = 20, # The frame rate of the exported video.
)

#=

"Analyzing exponential spread"

infected(x) = count(i == :I for i in x)
recovered(x) = count(i == :R for i in x)
# Aggregated data for number of infected and recovered indivisuals
adata = [(:status, infected), (:status, recovered)]

# Try different parameters
r1, r2 = 0.02, 0.05
β1, β2 = 0.05, 0.02
sir_model1 = sir_initiation(reinfection_probability = r1, βmax = β1)
sir_model2 = sir_initiation(reinfection_probability = r2, βmax = β1)
sir_model3 = sir_initiation(reinfection_probability = r1, βmax = β2)

data1, _ = run!(sir_model1, sir_agent_step!, sir_model_step!, 2000; adata)
data2, _ = run!(sir_model2, sir_agent_step!, sir_model_step!, 2000; adata)
data3, _ = run!(sir_model3, sir_agent_step!, sir_model_step!, 2000; adata)

# data1[(end-10):end, :] # shown partial of data

# using CairoMakie to create figure

figure = Figure()
ax = figure[1, 1] = Axis(figure; ylabel = "Infected", xlabel="Steps")
l1 = lines!(ax, data1[:, dataname((:status, infected))], color = :orange)
l2 = lines!(ax, data2[:, dataname((:status, infected))], color = :blue)
l3 = lines!(ax, data3[:, dataname((:status, infected))], color = :green)

figure[1, 2] = Legend(figure, [l1, l2, l3],
    ["r=$r1, beta=$β1", "r=$r2, beta=$β1", "r=$r1, beta=$β2"]
)
figure


"""
Adding Social distancing：

The best way to model social distancing is to make some agents simply not move (which feels like it approximates reality better).
"""

# sir_model = sir_initiation(isolated = 0.85)
# abmvideo(
#     "_static/socialdist5.mp4",
#     sir_model,
#     sir_agent_step!,
#     sir_model_step!;
#     title = "Social Distancing",
#     frames = 100,
#     spf = 2,
#     ac = sir_colors,
#     framerate = 20,
# )

r4 = 0.02
sir_model4 = sir_initiation(reinfection_probability = r4, βmax = β1, isolated = 0.85)

abmvideo(
    "_static/socialdist6.mp4",
    sir_model4,
    sir_agent_step!,
    sir_model_step!;
    title = "Social Distancing",
    frames = 2000,
    spf = 2,
    ac = sir_colors,
    framerate = 20,
)

data4, _ = run!(sir_model4, sir_agent_step!, sir_model_step!, 2000; adata)

l4 = lines!(ax, data4[:, dataname((:status, infected))], color = :red)

figure[1, 2] = Legend(figure, [l1, l2, l3, l4],
    ["r=$r1, beta=$β1", "r=$r2, beta=$β1", "r=$r1, beta=$β2", "r=$r4, social distancing"],
)
figure

=#

"Home Work task description"
#=
Question 1：

Please simulate a combination of 
- beta_max = 0.05, 0.04, 0.03, 0.02, 0.01 (representing personal hygiene: e.g. Wearing a mask) and
- isolated = 0.0, 0.5, 0.7, 0.8, 0.9. (Representing social distancing and/or lockdown measures) 
And thus there will be 25 simulations in total. 

Plot the number of infected individuals over 5000 steps with 5 beta_max parametegrs, each isolated parameter on a different plot. 
(That is, plot #1 is isolated = 0.0 with 5 time series: beta_max = 0.05, 0.04, 0.03, 0.02, 0.01, plot #2 is isolated = 0.5 with 5 time series: beta_max = 0.05, 0.04, 0.03, 0.02, 0.01, and so on.)
=#

#=
Question 2：

- List the maxima (peaks) of infected individuals in the 25 simulations. Which parameter set is the most effective in "flattening the curve" (having the lowest peak infected individuals)? 
- Compared to decreasing the beta_max parameter, is increasing the "isolated" parameter more effective in flattening the curve?
=#

infected(x) = count(i == :I for i in x)
recovered(x) = count(i == :R for i in x)
# Aggregated data for number of infected and recovered indivisuals
adata = [(:status, infected), (:status, recovered)]

"Try different parameters => isolated = 0.0, 0.5, 0.7, 0.8, 0.9"
hw5_isolated = 0.9

sir_model1 = sir_initiation(isolated = hw5_isolated, βmax = 0.05)
sir_model2 = sir_initiation(isolated = hw5_isolated, βmax = 0.04)
sir_model3 = sir_initiation(isolated = hw5_isolated, βmax = 0.03)
sir_model4 = sir_initiation(isolated = hw5_isolated, βmax = 0.02)
sir_model5 = sir_initiation(isolated = hw5_isolated, βmax = 0.01)

data1, _ = run!(sir_model1, sir_agent_step!, sir_model_step!, 5000; adata)
data2, _ = run!(sir_model2, sir_agent_step!, sir_model_step!, 5000; adata)
data3, _ = run!(sir_model3, sir_agent_step!, sir_model_step!, 5000; adata)
data4, _ = run!(sir_model4, sir_agent_step!, sir_model_step!, 5000; adata)
data5, _ = run!(sir_model5, sir_agent_step!, sir_model_step!, 5000; adata)

# data1[(end-10):end, :] # shown partial of data

# using CairoMakie to create figure

figure = Figure()
ax = figure[1, 1] = Axis(figure; ylabel = "Infected", xlabel="Steps", title = "SIR Model, isolated = $(hw5_isolated)", titlefontsize = 12)
# l1 = lines!(ax, data1[:, dataname((:status, infected))], color = :orange)
# l2 = lines!(ax, data2[:, dataname((:status, infected))], color = :blue)
# l3 = lines!(ax, data3[:, dataname((:status, infected))], color = :green)
# l4 = lines!(ax, data4[:, dataname((:status, infected))], color = :red)
# l5 = lines!(ax, data5[:, dataname((:status, infected))], color = :purple)
l1 = lines!(ax, data1[:, dataname((:status, infected))])
l2 = lines!(ax, data2[:, dataname((:status, infected))])
l3 = lines!(ax, data3[:, dataname((:status, infected))])
l4 = lines!(ax, data4[:, dataname((:status, infected))])
l5 = lines!(ax, data5[:, dataname((:status, infected))])

figure[1, 2] = Legend(figure, [l1, l2, l3, l4, l5],
    ["βmax = 0.05", "βmax = 0.04", "βmax = 0.03", "βmax = 0.02", "βmax = 0.01"]
)
figure


# find the maxima (peaks) of infected individuals in the 5 different βmax conditions.
peak1 = maximum(data1[:, dataname((:status, infected))])
peak2 = maximum(data2[:, dataname((:status, infected))])
peak3 = maximum(data3[:, dataname((:status, infected))])
peak4 = maximum(data4[:, dataname((:status, infected))])
peak5 = maximum(data5[:, dataname((:status, infected))])

print(
"Q2, maxima (peaks), isolated = $(hw5_isolated)：

- βmax = 0.05  ==> $(peak1)
- βmax = 0.04  ==> $(peak2)
- βmax = 0.03  ==> $(peak3)
- βmax = 0.02  ==> $(peak4)
- βmax = 0.01  ==> $(peak5)

"
)