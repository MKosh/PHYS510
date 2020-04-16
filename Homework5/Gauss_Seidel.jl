#=---------------------------------------
Poisson Equation Solver

Using:
    - Gauss-Seidel Method relaxation

---------------------------------------=#

using Plots


#-------------------- Start Function Declarations --------------------#
function ρ(x, y) # Charge Density Function
    if (40<x<=50 && 40<y<=50) || (50<x<=60 && 40<y<=50) # Ω1 & Ω2
        return 1 # using 1 and -1 instead of 50 like the book had
    elseif (40<x<=50 && 50<y<=60)|| (50<x<=60 && 50<y<=60) # Ω3 & Ω4
        return -1
    else
        return 0
    end
end

function grad(L, M) # Calculate the gradient at each point
    dfdx = zeros(M+1,M+1)
    dfdy = zeros(M+1,M+1)
    dfdx[2:M, 2:M] = L[2+1:M+1,2:M] - L[2-1:M-1,2:M]
    dfdy[2:M, 2:M] = L[2:M,2+1:M+1] - L[2:M,2-1:M-1]
    return (dfdx, dfdy)
end

function deconstruct!(L) # deconstruct the matrix L into an array
    temp = Float64[]
    for i in eachindex(L)
        push!(temp, L[i])
    end
    return temp
end

function meshReduce!(x, y, U, V) # Reduce the number of lines in the E plot
    Indices = []

    for i in 1:length(x)
        if x[i]%5.0 == 0.0 && y[i]%5.0 == 0.0 # plot lines every 5 x&y points
            push!(Indices, i)
        end
    end

    x = x[Indices]
    y = y[Indices]
    U = U[Indices]
    V = V[Indices]

    return x, y, U, V
end

meshgrid(x, y) = (repeat(x, outer=length(y)), repeat(y, inner=length(x)))

function Plot(x, y, U, V, phi)
    Eplot = quiver(x, y, quiver=(U,V), title = "Dipole Electic field")
    ϕplot = wireframe(phi, title = "Potential")
    plot(Eplot, ϕplot, layout = 2, size=(1280,720))
    #savefig("Dipole_w_E")
    gui()
end
#-------------------- End Function Declarations ---------------------#

#-------------------- Start Main --------------------#
function main()

    M = 100
    ϵ = 1e-4
    N = 10000
    phi = zeros(M+1,M+1)
    x, y = meshgrid(1:M+1, 1:M+1)
    U = Float64[] # Empty arrays to hold the x, and y components of the E field
    V = Float64[]
    Last = zeros(9801) # Array to hold the last time step values for phi

    @inbounds for x in 1:N
        V = [] # Empty array to hold the current phi values (instead of a matrix)
        for j in 2:M
            for i in 2:M
                phi[i,j] = (phi[i+1, j] + phi[i-1,j] + phi[i,j+1] + phi[i, j-1])/4 + ρ(i, j)/4
                push!(V, phi[i,j])
            end
        end
        diff = abs.(V - Last)
        delta = maximum(diff)
        if delta < ϵ # end the run after the accuracy reaches the desired level
            println("Done after ", x, " iterations")
            break
        end
        Last = copy(V)
    end

    u, v = grad(phi, M)
    U = deconstruct!(u)
    V = deconstruct!(v)

    x, y, U, V = meshReduce!(x, y, U, V)

    Plot(x, y, U, V, phi)
end

main()
#-------------------- End Main --------------------#

println("\nDone!\n")
