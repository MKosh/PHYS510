using Plots

const N = 100
const L = 30
const g = 9.8
const dt = 0.0001
const dt2 = dt^2
const time_steps = 5000
const m = 1.0
const σ = 1
const ϵ = 1


function Initialize(N)
    x = repeat(collect(10.0:2:28), 10)
    y = repeat(collect(10.0:2:28), inner=10)
    vx = zeros(N)
    vy = zeros(N)
    ax = zeros(N)
    ay = zeros(N)
    return x, y, vx, vy, ax, ay
end

function Force(dx, dy)
    global σ, ϵ, g
    a = σ/ϵ^2
    r2 = (dx^2 + dy^2)
    r2inv = 1/r2
    r2inv3 = r2inv^3
    f = 24 * a * r2inv3 * (2*r2inv3 - 1) * r2inv
    fx = f*dx
    fy = f*dy
    return fx, fy
end

function Verlet!(x, y, vx, vy, ax, ay)
    global L, dt, dt2
    for i in 1:N
        xnew = x[i] + vx[i]*dt + 0.5*ax[i]*dt2
        ynew = y[i] + vy[i]*dt + 0.5*ay[i]*dt2

        vx[i] = vx[i] + 0.5*ax[i]*dt
        vy[i] = vy[i] + 0.5*ay[i]*dt

        x[i], vx[i] = xBoundary!(xnew, vx[i], L)
        y[i], vy[i] = yBoundary!(ynew, vy[i], L)

    end

    ax, ay = Accel!(x, y, ax, ay)
    for i in 1:N
        vx[i] += 0.5*ax[i]*dt
        vy[i] += 0.5*ay[i]*dt
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
    scatter(x, y, xlims=(-L,2L), ylims=(0,2L), legend=false)
    gui()
end

function Accel!(x, y, ax, ay)
    global g
    for i in 1:N-1
        for j in i+1:N
            dx = x[i] - x[j]
            dy = y[i] - y[j]
            fx, fy = Force(dx, dy)
            ax[i] += fx
            ay[i] += fy - g

            ax[j] -= fx
            ay[j] -= fy + g
        end
    end
        return ax, ay
end

let
    x, y, vx, vy, ax, ay = Initialize(N)
    for i in 1:time_steps
        x, y = Verlet!(x, y, vx, vy, ax, ay)
        Plot(x, y, L)
        sleep(0.05)
    end
end

println("Doneish!\n")
