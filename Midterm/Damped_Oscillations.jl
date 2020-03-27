#=----------------------------------------------#
Damped Oscillations

Plots the angular position, angular velocity, & potential as a function of time
Plots the angular velocity as a function of the angle
Plots the shape of the potential well as a function of displacement

Using:
    - Runge Kutta 3/8 Method

Equations of Motion:
    ddot{θ} = -ω^2 * f(θ) - α * dot{θ}
    => dot{θ} = ν, dot{ν} = -ω^2*f(θ) - αν

    f(θ) = θ*H(|θ| - 1) where H is the Heaviside step function

-----------------------------------------------=#

using Plots

ω = 1
α = 0.05
dt = 0.01
N = 15000
θo = 3.0
vo = 0.0

#------------------- Start Function Declarations --------------------#
function f4(θ) # f_k given
    if (abs(θ) - 1) > 0
        return θ
    else
        return 0.0
    end
end

function Potential(θ) # Potential corresponding to f_4
    if (abs(θ) - 1) > 0
        return θ^2
    else
        return 0
    end
end

function F(θ, v)
    global α, ω
    dv = -ω^2*f4(θ) - α*v
    dθ = v
    return dθ, dv
end

function RK38!(θ, v) # RK 3/8 step
    global dt
    Y1 = [θ, v]
    F1 = F(Y1[1], Y1[2])
    Y2 = @. Y1 + dt*(F1)/3
    F2 = F(Y2[1], Y2[2])
    Y3 = @. Y1 + dt*(F2 - F1)/3
    F3 = F(Y3[1], Y3[2])
    Y4 = @. Y1 + dt*(F1 - F2 + F3)
    F4 = F(Y4[1], Y4[2])
    Y = @. Y1 + dt/8 * (F1 + 3F2 + 3F3 + F4)
    return Y
end

function Plot(Θ, V)
    x = -3:0.01:3
    line = zeros(length(Θ))

    lay = @layout [a ; [b c]]

    θ_plot = plot(line, color=:black, label="")
    plot!(test.(Θ), linestyle=:dot, color=:red, label="Potential")
    plot!(Θ, linewidth=1.5, c=:steelblue, gridalpha=0.5, ylims=(-3,4), yticks=(-3:4), framestyle=:box, title="Damped Oscillations", label="theta", xlabel="Time", ylabel="Angle, Angular Velocity, Potential")
    plot!(V, linewidth=1.5, c=:purple, label="v")

    phase = plot(Θ, V, linewidth=1.5, framestyle=:box, c=:maroon, xlabel="theta", ylabel="v", title="Phase Plot", legend=false)
    potential = plot(collect(x), Potential.(x), linewidth=1.5, framestyle=:box,c=:green, legend=false, xlabel="Displacement", ylabel="Potential", title="Potential vs Displacement")

    plot(θ_plot, phase, potential,  layout=lay, size=(1280,720))
    gui()
end
#------------------- End Function Declarations --------------------#

#------------------- Start Main --------------------#
let
    pos_vel = [θo, vo]
    Θ = [θo]
    V = [vo]
    for i in 1:N
        pos_vel = RK38!(pos_vel[1], pos_vel[2])
        push!(Θ, pos_vel[1])
        push!(V, pos_vel[2])
    end
    Plot(Θ, V)
end
#------------------- End Main --------------------#
println("\nDone!\n")
