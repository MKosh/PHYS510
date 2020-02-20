using Plots

const N = 100
const L = 30
const g = 9.8
const dt = 0.001
const time_steps = 5000
const m = 1.0
const σ = 1
const ϵ = 1



function Initialize(N)
    # Positions
    r = zeros(N, 2)
    x = repeat(collect(1.:10), 10)
    y = repeat(collect(1.:10), inner=10)
    r[:,1] = x .+ 10
    r[:,2] = y .+ 10

    # Forces

    #Velocities
    
    return r
end

r = Initialize(N)

display(r)
display(scatter(r[:, 1], r[:, 2], xlims = (0, 30), ylims = (0, 30)))
println("Done \n")
#function Force()
