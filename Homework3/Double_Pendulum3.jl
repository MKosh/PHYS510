using Plots

####################################################
# Double Pendulum simulation using RK4
# Can be run as an animation or just give static plots of the results

####################################################

# Initializing constants
const m1 = 1.0
const m2 = 3.0
const g = 9.8
const l1 = 1.0
const l2 = 2.0
const dt = 0.01

# Initializing variables
N = 10000
ϕ1o = 0.0
ϕ2o = 0.0
p1o = 4.0
p2o = 2.0
Title = "Poincare: phi1 = "*string(ϕ1o)*", phi2 = "*string(ϕ2o)*", p1 = "*string(p1o)*
    ", p2 = "*string(p2o)

# Creating arrays to hold ϕ1 and ϕ2, p1 and p2
Φ = [[] for i = 1:2]
P = [[] for i = 1:2]

# Pushing initial values to the arrays
push!(Φ[1], ϕ1o)
push!(Φ[2], ϕ2o)
push!(P[1], p1o)
push!(P[2], p2o)

function F(ϕ1, ϕ2, p1, p2)
    global m1
    global m2
    global l1
    global l2
    global g
    c = cos(ϕ1 - ϕ2)
    s = sin(ϕ1 - ϕ2)

    a = 1/(l1 * l2)
    M = (m1 + m2)
    P = p1*p2
    q = (1/(m1 + m2*s^2))

    b = a * P * s * q
    d = 1/2 * a^2 * q^2 * (l2^2*m2*p1^2 + l1^2*M*p2^2 - l1*l2*m2*P*c*sin(2(ϕ1-ϕ2)))

    f1 = q * (l2*p1 - l1*p2*c)/(l1^2*l2)
    f2 = q * (l1*M*p2 - l2*m2*p1*c)/(l1*l2^2*m2)
    f3 = -M * g * l1 * sin(ϕ1) - b + d
    f4 = -m2 * g * l2 * sin(ϕ2) + b - d

    return f1, f2, f3, f4
end

# Function to do the integration
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
    return Ynew # Returns an array of [ϕ1, ϕ2, p1, p2]
end


let
    # I don't remember what pp stands for, but it contains ϕ and p values
    pp = [ϕ1o, ϕ2o, p1o, p2o]

    # x and y positions of the masses on the pendulum
    x = [l1 * sin(ϕ1o), l2 * sin(ϕ2o)]
    y = [-l1 * cos(ϕ1o), -l1 * cos(ϕ1o) - l2 * cos(ϕ2o)]

    # x and y positions for the strings that hold the values
    xstring = [[] for i = 1:2]
    ystring = [[] for i = 1:2]

    # x and y values to trace the path of the masses
    xTrace = [[] for i = 1:2]
    yTrace = [[] for i = 1:2]

    # Initialize the arrays for the strings and traces
    push!(xstring[1], 0, x[1])
    push!(xstring[2], x[1], x[2])
    push!(ystring[1], 0, y[1])
    push!(ystring[2], y[1], y[2])

    push!(xTrace[1], x[1])
    push!(xTrace[2], x[2])
    push!(yTrace[1], y[1])
    push!(yTrace[2], y[2])

    # Custom layout for the plot window to show all three graphs in one window
    lay = @layout [[a{0.5h}; b{0.5h}] c]
    poincare = [[] for i in 1:2]
    # Loop over each step
    for i = 1:N
        # Calculating the ϕ and p values and pushing them to the appropriate array
        pp = RK4(pp[1], pp[2], pp[3], pp[4])
        push!(Φ[1], pp[1])
        push!(Φ[2], pp[2])
        push!(P[1], pp[3])
        push!(P[2], pp[4])

        # Updating the values that get plotted for the real space plot
        x1 = l1 * sin(Φ[1][i])
        x2 = x1 + l2 * sin(Φ[2][i])
        y1 = -l1 * cos(Φ[1][i])
        y2 = y1 - l2 * cos(Φ[2][i])

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

        # Cut the trace off after 5000 steps so it doesn't fill the whole plot
        # window
        if i > 5000
            deleteat!(xTrace[1], 1)
            deleteat!(xTrace[2], 1)
            deleteat!(yTrace[1], 1)
            deleteat!(yTrace[2], 1)
        end
    end

    conds = ["phi1="*string(ϕ1o) "phi2="*string(ϕ2o) "p1="*string(p1o) "p2="*string(p2o)]
    # Plot the results using the custom layout
    pend = plot(xTrace, yTrace,size=(1280,720), label=conds)
    plot!(pend, xstring, ystring, seriescolor = :black,label=[conds[3] conds[4]])
    scatter!(pend, x, y, label="")

    phi = plot(Φ[1], Φ[2], xlabel = "phi1", ylabel = "phi2", legend = false,size=(1280,720))
    plt = plot(P[1], P[2], xlabel = "p1", ylabel = "p2", legend = false,size=(1280,720))

    double_pend = plot(phi, plt, pend, layout = lay)
    display(double_pend)

    #savefig(double_pend, "Different_m_and_l")

    println("Doneish!\n")
end
