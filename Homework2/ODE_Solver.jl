using Plots

let
# The DE we are solving
f(y, t) = t^2 + y^2

# Initializing
t = 0.0
dt = 0.01
N = 40
y = 0.0
Title = "ODE Solver, y(0) = "*string(y)
# Arrays for holding the values
Time = [t]
Y_values = [[] for i in 1:4]    # An array of 4 arrays, one for each method

# Initialize the y(0) for each array
for i in 1:4
    push!(Y_values[i], y)
end

# The Adams-Bashford method relies on y_{n-1} the only thing I could think
# to do was to put a second value in the array, also 1,
# for that first call where n = 1. Otherwise there's a bounds error
push!(Y_values[4], y)

for i = 1:N
    #Euler
    y1new = Y_values[1][i] + f(Y_values[1][i], t)*dt
    push!(Y_values[1], y1new)

    #RK4
    Y1 = Y_values[2][i]
    Y2 = Y1 + dt/2 * f(Y1, t)
    Y3 = Y1 + dt/2 * f(Y2, t + t/2)
    Y4 = Y1 + dt * f(Y3, t + dt/2)
    y2new = Y1 + dt/6 * (f(Y1, t) + 2 * f(Y2, t + dt/2) + 2 * f(Y3, t+ dt/2) + f(Y4, t))
    push!(Y_values[2], y2new)

    #Taylor
    f_prime = 2 * Y_values[3][i]
    f_dot = 2 * t
    fn = f(Y_values[3][i], t)
    y3new = Y_values[3][i] + dt*fn + (dt^2)/2 * (f_dot + f_prime*fn)
    push!(Y_values[3], y3new)

    #Adams-Bashford
    y4new = Y_values[4][i+1] + dt/2 * (3*f(Y_values[2][i+1], t) - f(Y_values[2][i], t - dt))
    push!(Y_values[4], y4new)
    t += dt
end

Styles = [:solid :solid :dash :dash]
Labels = ["Euler" "RK4" "Taylor" "Adams-Bashford"]
plot(Y_values, xlabel="Time",ylabel="Y",label=Labels,linestyle=Styles,legend=:topleft,title=Title,size=(1280,720))
#savefig("ODEy0")
end
