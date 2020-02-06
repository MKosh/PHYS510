using Plots
let

function H(q1, q2, p1, p2)
    (p1^2 + p2^2)/2 - 1/(sqrt(q1^2 + q2^2))
end

e = 0.6
p1 = 0.0
q1 = 1-e

p2 = sqrt((1+e)/(1-e))
q2 = 0.0

t = 0.0
dt = 0.01
N = 1000000

Time = [t]
p1_values = [p1]
p2_values = [p2]
q1_values = [q1]
q2_values = [q2]

for i in 1:N
    p1new = p1 - (q1*dt)/((q1)^2 + (q2)^2)^(3/2)
    p2new = p2 - (q2*dt)/((q1)^2 + (q2)^2)^(3/2)
    q1new = q1 + p1*dt - (q1*dt^2)/(q1^2 +q2^2)^(3/2)
    q2new = q2 + p2*dt - (q2*dt^2)/(q1^2 +q2^2)^(3/2)
    p1 = p1new
    p2 = p2new
    q1 = q1new
    q2 = q2new

    push!(p1_values, p1)
    push!(p2_values, p2)
    push!(q1_values, q1)
    push!(q2_values, q2)
    if i%5 == 0
        plot(q1_values,q2_values,title="Symplectic Euler, e = 0.6",xlabel="q1",ylabel="q2",xlims=(-4,4),ylims=(-4,4))
        gui()
    end
    t += dt
    push!(Time, t)
end

#H_values = H.(q1_values, q2_values, p1_values, p2_values)
#display(plot(q1_values, q2_values, title = "Symplectic Explicit Euler, e = 0.6",xlabel="q1", ylabel="q2"))
#display(plot(p1_values, p2_values, title = "Symplectic Explicit Euler, e = 0.6",xlabel="p1", ylabel="p2"))
#display(plot(Time, H_values, title="Energy of the system",xlabel="Time",ylabel="Energy",ylims=(-1,0)))
end
