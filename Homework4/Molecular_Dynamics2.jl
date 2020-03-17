using Plots

const N = 100
const L = 30
const g = 9.8
const dt = 0.001
const dt2 = dt^2
const time_steps = 5000
const m = 1.0
const σ = 1
const ϵ = 1
const eps = 1e-8


function Initialize(N)
    global m, g
    x = repeat(collect(10.5:19.5), 10)
    y = repeat(collect(10.5:19.5), inner=10)
    vx = zeros(N)
    vy = zeros(N)
    ax = zeros(N)
    ay = zeros(N)
    ax1 = zero(N)
    ay1 = zeros(N)
    return x, y, vx, vy, ax, ay, ax1, ay1
end

function Force(dx, dy)
    global σ, ϵ, g, eps
    r2 = (dx^2 + dy^2)
    if r2 < eps
        r2 = eps
    end
    r2inv = 1.0/r2
    er2 = ϵ^2 * r2inv
    a = er2^3
    f = 24.0 * σ * r2inv * a * (2.0*a - 1.0)
    fx = f*dx
    fy = f*dy
    return fx, fy
end

function Verlet!(x, y, vx, vy, ax, ay, ax1, ay1)
    global L, dt, dt2
    for i in 1:N
        xnew = x[i] + vx[i]*dt + 0.5*ax[i]*dt2
        ynew = y[i] + vy[i]*dt + 0.5*ay[i]*dt2

        x[i], vx[i] = xBoundary!(xnew, vx[i], L)
        y[i], vy[i] = yBoundary!(ynew, vy[i], L)

    #    vx[i] = vx[i] + 0.5*ax[i]*dt
    #    vy[i] = vy[i] + 0.5*ay[i]*dt

    end

    ax1, ay1 = Accel!(x, y, ax, ay)
    for i in 1:N
        vx[i] += 0.5*(ax[i] + ax1[i])*dt
        vy[i] += 0.5*(ay[i] + ay1[i])*dt
    end
    return x, y
end

function xBoundary!(pos, v, L)
    if pos > L
        pos = L-(pos-L)
        v = -v
    elseif pos < 0
        pos = -pos
        v = -v
    end
    return pos, v
end

function yBoundary!(pos, v, L)
    if pos < 0
        pos = -pos
        v = -v
    end
    return pos, v
end

function Plot(x, y, L)
    scatter(x, y, xlims=(0,L), ylims=(0,L), legend=false)
    gui()
end

function Accel!(x, y, ax, ay)
    global g, m
    ax .= 0.0
    ay .= -m*g

    for i in 1:N
        for j in i+1:N
            dx = x[i] - x[j]
            dy = y[i] - y[j]
            fx, fy = Force(dx, dy)
            ax[i] += fx
            ay[i] += fy

            ax[j] -= fx
            ay[j] -= fy
        end
    end
        return ax, ay
end

let
    x, y, vx, vy, ax, ay, ax1, ay1 = Initialize(N)
    for i in 1:time_steps
        ax, ay = Accel!(x, y, ax, ay)
        x, y = Verlet!(x, y, vx, vy, ax, ay, ax1, ay1)
        Plot(x, y, L)
        #sleep(0.05)
    end
end

println("Doneish!\n")
