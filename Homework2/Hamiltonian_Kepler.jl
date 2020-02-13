using Plots
let

function H(q1, q2, p1, p2)
    (p1^2 + p2^2)/2 - 1/(sqrt(q1^2 + q2^2))
end

e = 0.9
Title = "Explicit Euler, e = "*string(e)
p1 = 0.0
q1 = 1-e

p2 = sqrt((1+e)/(1-e))
q2 = 0.0

t = 0.0
dt = 0.01
N = 10000

Time = [t]
p1_values = [p1]
p2_values = [p2]
q1_values = [q1]
q2_values = [q2]

for i in 1:N
    p1new = p1 - (q1*dt)/((q1)^2 + (q2)^2)^(3/2)
    p2new = p2 - (q2*dt)/((q1)^2 + (q2)^2)^(3/2)
    q1new = q1 + p1*dt
    q2new = q2 + p2*dt
    p1 = p1new
    p2 = p2new
    q1 = q1new
    q2 = q2new

    push!(p1_values, p1)
    push!(p2_values, p2)
    push!(q1_values, q1)
    push!(q2_values, q2)

    t += dt
    push!(Time, t)
end

H_values = H.(q1_values, q2_values, p1_values, p2_values)

q_plt = plot(q1_values, q2_values,title=Title, xlabel="q1", ylabel="q2", aspect_ratio=1,legend=false,size=(1280,720))
p_plt = plot(p1_values, p2_values,title=Title, xlabel="p1", ylabel="p2", aspect_ratio=1,legend=false,size=(1280,720))
H_plt = plot(Time, H_values,title = "Energy of the system, e = 0.9", xlabel="Time", ylabel="Energy",legend=false,ylims=(-1,0),size=(1280,720))

savefig(q_plt, "HEqe0_9")
savefig(p_plt, "HEpe0_9")
savefig(H_plt, "HEHe0_9")
end
