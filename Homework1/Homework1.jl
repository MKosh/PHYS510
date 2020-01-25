using Plots

# Initializing
a = 5.0
N = 10000
h = (a+a)/(N)
x = -a:h:a # Generate a step range from -a to a in steps of length h
# Initializing

# Starting function declarations
# Functions for each of the finite difference methods
# Inputs are the function we're finding the derivative of, and the point we're differentiating around
function forward_diff(func, x)
    global h # use the global h instead of manually putting it in for every function call
    (func(x + h) - func(x))/(h)
end

function back_diff(func, x)
    global h
    (func(x) - func(x - h))/(h)
end

function central_diff(func, x)
    global h
    (func(x + h) - func(x - h)) / 2h
end

function central_diff2(func, x)
    global h
    (func(x + h) - 2func(x) + func(x - h))/h^2
end

function relative_error(y1, y2)
    abs((y1 - y2) / y1)
end

g(x) = exp(-x^2)
f(x) = sin(x)

d1g_analytic(x) = -2x*exp(-x^2) # Analytic expression for dg/dx
d1f_analytic(x) = cos(x) # Analytic expression for df/dx
# End function declarations

# Generate an array, calculate the values of each function over the step range declared above, and store them in the array
g1 = @. g(x)
d1g = @. d1g_analytic(x)
d1g_forward = @. forward_diff(g, x)
d1g_backward = @. back_diff(g, x)
d1g_central = @. central_diff(g, x)

f1 = @. f(x)
d1f = @. d1f_analytic(x)
d1f_forward = @. forward_diff(f, x)
d1f_backward = @. back_diff(f, x)
d1f_central = @. central_diff(f, x)
d2f_central = @. central_diff2(f, x)

d1g_forward_error = @. relative_error(d1g, d1g_forward)
d1g_backward_error = @. relative_error(d1g, d1g_backward)
d1g_central_error = @. relative_error(d1g, d1g_central)

d1f_forward_error = @. relative_error(d1f, d1f_forward)
d1f_backward_error = @. relative_error(d1f, d1g_backward)
d1f_central_error = @. relative_error(d1f, d1g_central)

# Collect all of the results into one set of arrays for each function
g_values = hcat(g1, d1g, d1g_forward, d1g_backward, d1g_central)
f_values = hcat(f1, d1f, d1f_forward, d1f_backward, d1f_central, d2f_central)

g_errors = hcat(d1g_forward_error, d1g_backward_error, d1g_central_error)
f_errors = hcat(d1f_forward_error, d1f_backward_error, d1f_central_error)

lines = [:solid :solid :solid :solid :dash :solid]
g_labels = ["g(x)" "g'(x) analytic" "g'(x) forward" "g'(x) backward" "g'(x) central"]
f_labels = ["f(x)" "f'(x) analytic" "f'(x) forward" "f'(x) backward" "f'(x) central" "f''(x) central"]

g_error_labels = ["g'(x) forward error" "g'(x) backward error" "g'(x) central error"]
f_error_labels = ["f'(x) forward error" "f'(x) backward error" "f'(x) central error"]

g_plt = plot(x, g_values, label=g_labels,title="g(x) = exp(-x^2)", xlabel="x", ylabel="y", line=lines)
f_plt = plot(x, f_values, label=f_labels, title = "f(x) = sin(x)", xlabel="x", ylabel="y", )

g_error = plot(x, g_errors, label=g_error_labels, title = "Relative Error g'(x)")
f_error = plot(x, f_errors, label=f_error_labels, title = "Relative Error f'(x)")

display(g_plt)
display(f_plt)
display(g_error)
display(f_error)

println("Done!")
