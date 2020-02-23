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
    R = zeros(N, 2)
    x = repeat(collect(1.:10), 10)
    y = repeat(collect(1.:10), inner=10)
    R[:,1] = x .+ 10
    R[:,2] = y .+ 10

    # Forces
    F = zeros(N, 2)
    
    #Velocities
    V = zeros(N, 2)

    return R, F, V
end

function Force()

end

function Verlet()

end


R, F, V = Initialize(N)

display(scatter(R[:, 1], R[:, 2], xlims = (0, 30), ylims = (0, 30)))
println("Done \n")
#function Force()
