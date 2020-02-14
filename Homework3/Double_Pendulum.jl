using Plots

const m = 1.0
const g = 9.8
const l = 1.0
const dt = 0.01

N = 10000
ϕ1o = 0.0
ϕ2o = 0.0
p1o = 0.0
p2o = 6.5

Φ = [[] for i in 1:2]
P = [[] for i in 1:2]
Y_values = [[] for i in 1:4]

push!(Φ[1], ϕ1o)
push!(Φ[2], ϕ2o)
push!(P[1], p1o)
push!(P[2], p2o)

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

x = [l*sin(ϕ1o), l*sin(ϕ2o)]
y = [-l*cos(ϕ1o), -l*cos(ϕ1o) - l*cos(ϕ2o)]
xstring = [[] for i in 1:2]
ystring = [[] for i in 1:2]
xTrace = [[] for i in 1:2]
yTrace = [[] for i in 1:2]

push!(xstring[1], 0,x[1])
push!(xstring[2], x[1],x[2])
push!(ystring[1], 0,y[1])
push!(ystring[2], y[1],y[2])

push!(xTrace[1], x[1])
push!(xTrace[2], x[2])
push!(yTrace[1], y[1])
push!(yTrace[2], y[2])

lay = @layout [[a{0.5h} ; b{0.5h}] c]

for i in 1:N
    pp = RK4(pp[1], pp[2], pp[3], pp[4])
    push!(Φ[1], pp[1])
    push!(Φ[2], pp[2])
    push!(P[1], pp[3])
    push!(P[2], pp[4])

    x1 = l*sin(Φ[1][i])
    x2 = x1 + l*sin(Φ[2][i])
    y1 = -l*cos(Φ[1][i])
    y2 = y1 - l*cos(Φ[2][i])

    x[1] = x1
    x[2] = x2
    y[1] = y1
    y[2] = y2

    xstring[1][2] = x1
    xstring[2][1] = x1
    xstring[2][2] = x2
    ystring[1][2] = y1
    ystring[2][2] = y2
    ystring[2][1] = y1

    push!(xTrace[1], x1)
    push!(xTrace[2], x2)
    push!(yTrace[1], y1)
    push!(yTrace[2], y2)

    if i > 5000
        deleteat!(xTrace[1], 1)
        deleteat!(xTrace[2], 1)
        deleteat!(yTrace[1], 1)
        deleteat!(yTrace[2], 1)
    end
    #pend = scatter(x,y, xlims=(-2l,2l),ylims=(-3l,0.5))
    #plot!(pend, xstring, ystring, seriescolor=:black)
    #plot!(pend, xTrace, yTrace)

    #phi = plot(Φ[1], Φ[2], xlabel="phi1",ylabel="phi2",legend=false)
    #plt = plot(P[1], P[2], xlabel="p1",ylabel="p2",legend=false)

    #plot(phi, plt, pend, layout=lay)
    #gui()
end

pend = scatter(x,y, xlims=(-2l,2l),ylims=(-2.5l,0.5), legend=false)
plot!(pend, xstring, ystring, seriescolor=:black)
plot!(pend, xTrace, yTrace)

phi = plot(Φ[1], Φ[2], xlabel="phi1",ylabel="phi2",legend=false)
plt = plot(P[1], P[2], xlabel="p1",ylabel="p2",legend=false)

plot(phi, plt, pend, layout=lay)
gui()

println("Doneish!")

end
