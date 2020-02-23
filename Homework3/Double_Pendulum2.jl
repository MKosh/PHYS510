using Plots, DelimitedFiles

# Initializing values
const m = 1.0
const g = 9.8
const l = 1.0
const dt = 0.01

N = 10000
ϕ1o = 0.0
ϕ2o = 0.0
p1o = 0.0
p2o = 1.5

Time = collect(0.0:dt:(N*dt)-dt)

function F(ϕ1, ϕ2, p1, p2)
    global m
    global l
    global g
    c = cos(ϕ1 - ϕ2)
    s = sin(ϕ1 - ϕ2)
    d = 1 / (m * l^2)
    a = 1 / (1 + s^2)

    f1 = d * a * (p1 - p2 * c)
    f2 = d * a * (2p2 - p1 * c)
    f3 = d * a * (-p1 * p2 * s + a * (p1^2 + 2 * p2^2 - 2p1 * p2 * c) * c * s) -
         2 * m * g * l * sin(ϕ1)
    f4 = d * a * (p1 * p2 * s - a * (p1^2 + 2 * p2^2 - 2p1 * p2 * c) * s * c) -
         m * g * l * sin(ϕ2)

    return f1, f2, f3, f4
end


function RK4(ϕ1, ϕ2, p1, p2)
    global m
    global l
    global g

    Y1 = [ϕ1, ϕ2, p1, p2]
    F1 = F(ϕ1, ϕ2, p1, p2)
    Y2 = @. Y1 + 0.5 * dt * F1
    F2 = F(Y2[1], Y2[2], Y2[3], Y2[4])
    Y3 = @. Y1 + 0.5 * dt * F2
    F3 = F(Y3[1], Y3[2], Y3[3], Y3[4])
    Y4 = @. Y1 + dt * F3
    F4 = F(Y4[1], Y4[2], Y4[3], Y4[4])
    Ynew = @. (Y1 + dt * (F1 + 2 * (F2 + F3) + F4) / 6)
    return Ynew
end

function draw_plot(dist, time)
    display(plot(time, dist, xlabel="Time",ylabel="Distance",title="Something",legend=false))
end

Dist(ϕ1, ϕ2, ϕ1p, ϕ2p, p1, p2, p1p, p2p) =
    sqrt((ϕ1 - ϕ1p)^2 + (ϕ2 - ϕ2p)^2 + (p1 - p1p)^2 + (p2 - p2p)^2)



    pp1 = [ϕ1o, ϕ2o, p1o, p2o]
    pp2 = [ϕ1o, ϕ2o, p1o+0.1, p2o]

    function pend(pp, dt)
        Φ = [[] for i in 1:2]
        P = [[] for i in 1:2]
        for i = 1:N
            pp = RK4(pp[1], pp[2], pp[3], pp[4])
            push!(Φ[1], pp[1])
            push!(Φ[2], pp[2])
            push!(P[1], pp[3])
            push!(P[2], pp[4])
        end
        return Φ, P
    end

Phi1, Mom1 = pend(pp1, dt)
Phi2, Mom2 = pend(pp2, dt)
dist = Dist.(Phi1[1], Phi1[2], Phi2[1], Phi2[2], Mom1[1], Mom1[2], Mom2[1],
    Mom2[2])
draw_plot(dist, Time)


    println("Doneish! \n")
