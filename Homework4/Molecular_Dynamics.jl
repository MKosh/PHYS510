using Plots

const N = 100
const L = 30
const g = 9.8
const dt = 0.001
const time_steps = 5000
const m = 1.0
const σ = 1
const ϵ = 1
const rows = 10

function Initialize(N, rows)
    global g
    global m
    R = zeros(N, 2)
    x = repeat(collect(1.:rows), rows)
    y = repeat(collect(1.:rows), inner=rows)
    R[:,1] = x .+ 10
    R[:,2] = y .+ 10

    #F = zeros(N, 2)
    #F[:,2] .= -m*g
    F = [0.0 -9.8]

    V = zeros(N, 2)

    return R, F, V
end

function dist(R, particle)
    distance = zeros(N)
    temp = [R[particle,1] R[particle,2]]
    Z = repeat(temp, N)

    C = Z - R
    dir = zeros(N, 2)
    for i in 1:N
        if i == particle
            distance[i] = 0.0
            dir[i,1] = 0.0
            dir[i,2] = 0.0
        else
            distance[i] = sqrt(C[i, 1]^2 + C[i, 2]^2)
            dir[i,1] = C[i,1]/distance[i]
            dir[i,2] = C[i,2]/distance[i]
        end
    end

    #temp = findall(x->x<1e-8, distance)
    #for i in 1:length(temp)
    #    distance[temp[i]] *=1000
    #end

    return distance, dir
end



function LJ(r, rhat)
    global σ, ϵ
    ρ = @. 1/abs(r)^2
    f = @. 24*σ*ρ * (2*(ϵ^2*ρ)^6 - (ϵ*ρ)^3) .* rhat
    temp = findall(isnan, f)
    f[temp[1]] = 0.0
    f[temp[2]] = 0.0
    return f
end


function find_force!(R, particle)
    global m, g
    Fext = zeros(N, 2)
    Fext[:,2] .= -m*g

    r, rhat = dist(R, particle)
    Fij = LJ(r, rhat)
    Force = [0.5*(sum(Fij[:,1])) (0.5*sum(Fij[:,2])-m*g*10)]
    return Force
end

function update_position!(R, V, F, part)
    global dt
    R[part,:] = @. R[part,:] + V[part,:]*dt + 0.5*F[:]*dt^2
    return R
end

function update_velocity!(R, V, F, newF, part)
    global dt
    V[part,:] = @. V[part,:] + 0.5*(F[:] + newF[:])*dt
    return V
end

function boundaries!(R, V, part)
    global L
        if R[part, 1] < 0
            R[part, 1] = -R[part,1]
            V[part, 1] = -V[part,1]
        elseif R[part, 1] > L
            R[part, 1] = L-(R[part,1] -L)
            V[part, 1] = -V[part, 1]
        end
        if R[part, 2] < 0
            R[part, 2] = -R[part, 2]
            V[part, 2] = -V[part, 2]
        end
    #temp = findall(x->x>L, R[:,1])
    #R[temp] = @. L - (R[temp] - L)
    #V[temp] = -V[temp]
    #temp = findall(x->x<0, R[:,1])
    #R[temp] = -R[temp]
    #V[temp] = -V[temp]
    #temp = findall(x->x<0, R[:,2])
    #R[temp] = -R[temp]
    #V[temp] = -V[temp]
    return R, V
end

let
pos, forces, vel = Initialize(N, rows)

for j in 1:time_steps

    for i in 1:N
        Fk = forces
        forces = find_force!(pos, i)
        vel = update_velocity!(pos, vel, Fk, forces, i)
            pos = update_position!(pos, vel, forces, i)
            pos, vel = boundaries!(pos, vel, i)
    end
scatter(pos[:,1],pos[:,2], xlims=(0,L),ylims=(0,2L),legend=false)
gui()
end
#display(scatter(pos[:,1],pos[:,2],xlims=(0,L),ylims=(0,L)))

end
println("Done!\n")
