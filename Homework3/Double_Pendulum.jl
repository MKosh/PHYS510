using Plots

const m = 1.0
const g = 9.8
const l = 1.0
const dt = 0.01

N = 10000
ϕ1o = 1.0
ϕ2o = 0.0
p1o = 0.0
p2o = 3.0

Φ = [[] for i in 1:2]
P = [[] for i in 1:2]
Y_values = [[] for i in 1:4]

push!(Φ[1], ϕ1o)
push!(Φ[2], ϕ2o)
push!(P[1], p1o)
push!(P[2], p2o)

push!(Y_values[1], ϕ1o)
push!(Y_values[2], ϕ2o)
push!(Y_values[3], p1o)
push!(Y_values[4], p2o)

function F(ϕ1, ϕ2, p1, p2)
    global m
    global l
    global g
    c = cos(ϕ1 - ϕ2)
    s = sin(ϕ1 - ϕ2)
    d = 1/(m*l^2)
    a = 1/(1 + s^2)

    f1 = d * a * (p1 - p2 * c)
    f2 = d * a * (2p2 - p1 * c)
    f3 = d * a * (-p1*p2*s + a*(p1^2 + 2p2^2 - 2p1*p2*c)*c*s) - 2*m*g*l*sin(ϕ1)
    f4 = d * a * (p1*p2*s - a *(p1^2 + 2p2^2 - 2p1*p2*c)*s*c) - m*g*l*sin(ϕ2)

    return f1, f2, f3, f4
end


function RK4(ϕ1, ϕ2, p1, p2)
    global m
    global l
    global g

    Y1 = [ϕ1, ϕ2, p1, p2]
    F1 = F(ϕ1, ϕ2, p1, p2)
    Y2 = @. Y1 + 0.5*dt*F1
    F2 = F(Y2[1], Y2[2], Y2[3], Y2[4])
    Y3 = @. Y1 + 0.5*dt*F2
    F3 = F(Y3[1], Y3[2], Y3[3], Y3[4])
    Y4 = @. Y1 + dt * F3
    F4 = F(Y4[1], Y4[2], Y4[3], Y4[4])
    Ynew = @. (Y1 + dt * (F1 + 2*(F2+F3) + F4)/6)
    return Ynew
end


let
pp = [ϕ1o, ϕ2o, p1o, p2o]

for i in 1:N
    pp = RK4(pp[1], pp[2], pp[3], pp[4])
    push!(Φ[1], pp[1])
    push!(Φ[2], pp[2])
    push!(P[1], pp[3])
    push!(P[2], pp[4])
end

display(plot(Φ[1], Φ[2], xlabel="phi1",ylabel="phi2",legend=false))
display(plot(P[1], P[2], xlabel="p1",ylabel="p2",legend=false))

println("Done!")

end
