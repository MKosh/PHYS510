#=--------------------------------------
Time Dependent Heat Equation

Using:
    - Explicit Euler Scheme

Current task:
    - Optimization
        - Putting everything into a function removes globals and speeds things up
        - The for loop in main allocates almost 70GB of data for 1million time
        steps and takes about 15 seconds to run. Adding the @. macro cuts allocations
        in half, but there's probably more I could do.
--------------------------------------=#

using Plots, Profile, BenchmarkTools

#-------------------- Start Function Declarations --------------------#
# Nothin' to see here, keep on scrolling
#-------------------- End Function Declarations ---------------------#

#-------------------- Start Main --------------------#
function main()
    k = 1.0
    dt = 0.5
    L = 100
    T0 = 0.0
    TN = 2.0
    N = 1000
    h = 1.05
    c = k*dt/h^2
    n = 1000000

    times = [0, Int64(floor(n/6)), Int64(floor(2n/6)), Int64(floor(3n/6)), Int64(floor(4n/6)), Int64(floor(5n/6)), n]

    T, Tnew = zeros(N), zeros(N)
    T[1], Tnew[1] = T0, T0
    T[N], Tnew[N] = TN, TN
    Temps = []

    for t in 0:n
        @. Tnew[2:N-1] = T[2:N-1] + c*(T[1:N-2] + T[3:N] - 2*T[2:N-1])
        T = copy(Tnew)
        if t == times[1] || t == times[2] || t == times[3] || t == times[4] || t == times[5] || t == times[6] || t == times[7]
            push!(Temps, T)
        end
    end

    x = collect(range(0, stop=L, length=N))

    plot(x, Temps, legend = false, xlabel = "x", ylabel = "T(x)", size = (1280, 720))
    #savefig("Heat_Equation")
    gui()
end

main()
#-------------------- End Main --------------------#

println("\nDone\n")
