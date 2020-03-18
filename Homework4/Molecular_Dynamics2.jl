#=--------------------------------------------------------#
Molecular Dynamics Simulation

Using:
    - Lennard-Jones Potential
    - Velocity Verlet algorithm

Options:
    -Initial shape of the particles
    -Saves a gif of the animation to the appdata/local/temp folder
    -Run a slow animation without waiting for a gif to be made
        -Particle postions
        -Temperature
        -Average Temperature
    -Display a static plot of the Temp/AveTemp and the final
        state of the particle simulation
#--------------------------------------------------------=#

using Plots

const L = 30
const g = 9.8
const dt = 0.001
const dt2 = dt^2
const time_steps = 5000
const m = 1.0
const σ = 1
const ϵ = 1
const eps = 1e-8
const dims = 2
const Kb = 1.38e-23
const C = 1/(dims*(N-1))

#------------------- Start Function Declarations --------------------#
function Initialize(N, shape)  # Initialize the variables for position, velocity, and force
    global m, g, ϵ
    if shape == "box"
        x = repeat(collect(20.5:29.5), 10)    # Set up an array for the initial x components
        y = repeat(collect(10.5:19.5), inner=10)  # Set up an array for the inital y components
    elseif shape == "triangle"
        N = 91 # Change N to 91 for symmetry
        x = zeros(N)
        y = zeros(N)
        for i in 1:13
            for j in 1:i
                pos = collect((L/2-(i+1)):2:(L/2+(i-1)))
                x[(sum(1:i-1)+j)] = pos[j]
                y[(sum(1:i-1)+j)] = 14-i
            end
        end
        x = x .+ 2
        y = y .+ 10
    end
    vx = zeros(N)   # Arrays of length N, all elements are zero
    vy = zeros(N)
    ax = zeros(N)   # F = ma, but for our case m = 1
    ay = zeros(N)
    ax1 = zero(N)
    ay1 = zeros(N)
    T = []
    Tave = []
    return x, y, vx, vy, ax, ay, ax1, ay1, T, Tave, N
end

function LJ(dx, dy) # Calculate the force on the particles from the Lennard-Jones potential
    global σ, ϵ, g, eps
    r2 = (dx^2 + dy^2)

    if r2 < eps # If the distance between two particles is smaller than eps, set the distance to eps
        r2 = eps
    end

    r2inv = 1.0/r2
    er2 = ϵ^2 * r2inv
    a = er2^3
    f = 24.0 * σ * r2inv * a * (2.0*a - 1.0)
    fx = f*dx
    fy = f*dy
    return fx, fy   # Return the x and y components of the force
end

function Verlet!(x, y, vx, vy, ax, ay, ax1, ay1, N)    # Implement the velocity Verlet method
    global L, dt, dt2
    for i in 1:N    # Update the positions of the particles, check boundary conditions and correct as needed
        xnew = x[i] + vx[i]*dt + 0.5*ax[i]*dt2
        ynew = y[i] + vy[i]*dt + 0.5*ay[i]*dt2

        x[i], vx[i] = xBoundary!(xnew, vx[i], L)
        y[i], vy[i] = yBoundary!(ynew, vy[i], L)

    end

    ax1, ay1 = Accel!(x, y, ax, ay, N) # Calculate the updated force acting on the particles based on their new positions
    for i in 1:N    # Update the velocities of the particles
        vx[i] += 0.5*(ax[i] + ax1[i])*dt
        vy[i] += 0.5*(ay[i] + ay1[i])*dt
    end
    return x, y # Return the position vectors
end

function xBoundary!(pos, v, L)  # Check the boundary conditions in the x direction (Nonperiodic boundaries)
    if pos > L
        pos = L-(pos-L)
        v = -v
    elseif pos < 0
        pos = -pos
        v = -v
    end
    return pos, v
end

function yBoundary!(pos, v, L)  # Check the boundary condition in the y direction (No top on the box)
    if pos < 0
        pos = -pos
        v = -v
    end
    return pos, v
end

function Plot(x, y, L, T, Tave, type)  # Plot the positions of the particles. Used for real time visualization
    #=This function plots the particle simulation, temperature at each time step,
    or the average temperature over 100 time step intervals.
    It's mostly just for testing, to check if things are running properly.
    It runs very slowly when run in real time.
    The main plotting is handled by the @gif macro so the animation can be saved, and
    it looks much smoother.
    =#
    global time_steps
    if type == "particle"
        particles = scatter(x, y, xlims=(0,L), ylims=(0,L), legend=false, size=(1280,720))
        plot(particles)
        gui()
    elseif type == "Temp"
        plot(T, xlims=(0,time_steps), ylims=(0,9000),legend=false, xlabel=("Time"),ylabel=("Temperature"), size=(1280,720))
        gui()
    elseif type == "Tave"
        plot(Tave, xlabel=("Time"), ylabel=("Average Temperature"),title="Average Temp (Intervals of 100 Time Steps)", legend=false, size=(1280,720))
        gui()
    elseif type == "all"
        lay = @layout([[a ; b] c])
        p1 = plot(T, xlims=(0,time_steps), ylims=(0,9000),legend=false, xlabel=("Time"),ylabel=("Temperature"))
        p2 = plot(Tave, xlabel=("Time"), ylabel=("Average Temperature"),title="Average Temp (Intervals of 100 Time Steps)", legend=false)
        p3 = scatter(x, y, xlims=(0,L), ylims=(0,L), legend=false)
        plot(p1, p2, p3, layout=lay, size=(1280,720))
        gui()
    end
end

function Accel!(x, y, ax, ay, N)   # Calculate the force/acceleration on each particle
    global g, m
    ax .= 0.0   # Set x componet of force to 0
    ay .= -m*g  # Add gravitational force in the y direction

    for i in 1:N-1    # Sum over all particles
        for j in i+1:N
            dx = x[i] - x[j]    #Calculate the x component of the distance between particles i and j
            dy = y[i] - y[j]    # y component
            fx, fy = LJ(dx, dy) #Calculate the force from the Lennard-Jones potential
            ax[i] += fx     # Add the forces on the ith particle from the jth particle
            ay[i] += fy

            ax[j] -= fx     # Force from the ith on the jth particle is opposite to the force above
            ay[j] -= fy
        end
    end
        return ax, ay
end

function Temperature!(T, vx, vy)
    global C
    push!(T, C*(sum(vx)^2 + sum(vy)^2))
end

function AveTemp!(Tave, T, i)
    push!(Tave, 0.01*sum(T[i-99:i]))
end
#------------------- End Function Declarations --------------------#

#------------------- Start Main Loop --------------------#
let
    N = 100 # Number of particles, can't be declared a const because it changes for triangles
    shape = "box" # Initial shape of the particles box or triangle
    Gif = false # Only set to true if you want to save a gif of the simulation
    Anim = false # Set to true for a painfully slow realtime animation, anything else gives a static plot at the end of the time steps.
    if (Gif == true)
        lay = @layout [a b]
        x, y, vx, vy, ax, ay, ax1, ay1, T, Tave, N = Initialize(N, shape)
        #p = plot(plot([], xlims=(0,time_steps),ylims=(0,9000),xlabel="Time",
            #ylabel="Temperature"), scatter(x, y,xlims=(0,L),ylims=(0,L),legend=false),
            #size=(1280,720),layout=lay) #$
            p = scatter(x, y, xlims=(0,L), ylims=(0,L), legend=false, size=(1280,720)) #*
        anim = Animation() # Set up the animaion.
        @gif for i in 1:time_steps
            ax, ay = Accel!(x, y, ax, ay, N)
            x, y = Verlet!(x, y, vx, vy, ax, ay, ax1, ay1, N)
            Temperature!(T, vx, vy)
            #push!(p, T[i]) #$
            #p[2][1][:x] = x #$
            #p[2][1][:y] = y #$
            p[1][1][:x] = x #*
            p[1][1][:y] = y #*
        end every 75
        #* Uncomment these lines and comment out the lines with $ after them for the animation of just the particles
    else
        if Anim == true
            x, y, vx, vy, ax, ay, ax1, ay1, T, Tave, N = Initialize(N, shape)
            for i in 1:time_steps
                ax, ay = Accel!(x, y, ax, ay, N)
                x, y = Verlet!(x, y, vx, vy, ax, ay, ax1, ay1, N)
                #Temperature!(T, vx, vy)
                #if i % 100 == 0
                #    AveTemp!(Tave, T, i)
                #end
                Plot(x, y, L, T, Tave, "particle")
            end
        else
            x, y, vx, vy, ax, ay, ax1, ay1, T, Tave, N = Initialize(N, shape)
            for i in 1:time_steps
                ax, ay = Accel!(x, y, ax, ay, N)
                x, y = Verlet!(x, y, vx, vy, ax, ay, ax1, ay1, N)
                Temperature!(T, vx, vy)
                if i % 100 == 0
                    AveTemp!(Tave, T, i)
                end
            end
            Plot(x, y, L, T, Tave, "all")
        end
    end
end
#------------------- End Main Loop --------------------#

println("\nDone!\n")
