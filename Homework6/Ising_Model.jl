#=--------------------------------------------------#
The 2 Dimensional Ising Model

A visualization of the lattice and graphs of the avg Magnetization and Energy
Using:
	- Monte Carlo technique and the Metropolis algorithm
	- Default Julia RNG Mersenne Twister

Method:

	Start with generating a 2D NxN lattice populated by random spins up or down (+1 or -1)

	Account for boundary conditions (In this case toroidal) if at the top row of lattice
	"up" is at the bottom of the lattice, if at rightmost column of lattice "right" is the
	left most column. Vice versa for the bottom row and left column

	Only certain values of energy changes (dE) are possible based on the surrounding spin states
	Precompute the exponents e^(-dE/T) as e4 and e8
	Create a function w(t) to determine which exponent to use instead of recalculating them for each MCstep

	Go to a random lattice site and determine the energy change associated with a spin flip
	dE is obtained by summing the nearest 4 neighbors and multiplying by 2*L[i, j] (the lattice site in question)
	if dE < 0 accept the change and calculate resulting values and graph the new lattice
	if dE >= 0 do not accept the flip, check the metropolis condition
	if r (a randomly generated number from 0 to 1) < w(t) accept the flip anyway (this is the metropolis part,
	the flip comes from the energy absorbed from the heat bath)
	otherwise the spin stays its original value

	Run the program over a large number of Monte Carlo Steps and watch the lattice go

---------------------------------------------------=#
using Plots, Profile
gr()


#-------------------- Start Function Declarations --------------------#

#function for possible energy changes
function w(t)
    if t == 8 || t == -8
        return e8
    elseif t == 4 || t == -4
        return e4
    elseif t == 0
        return 0
    end
end

#function for calculating the initial energy
function CalcEnergy(latt, N, J)
	energy = 0.0
    for i in 1:N
	    for j in 1:N
		energy -= DeltaE(latt, N, i, j, J)
	    end
     end
     return energy
end
#function for counting spin up vs spin down


#Function to sum the four nearest neighbors, accounting for toroidal boundaries
function DeltaE(latt, N, i, j, J)
    if i == 1
        left = latt[N, j]
    else
        left = latt[i-1,j]
    end
    if i == N
        right = latt[1, j]
    else
        right = latt[i+1, j]
    end
    if j == N
        up = latt[i, 1]
    else
        up = latt[i, j+1]
    end
    if j == 1
        down = latt[i, N]
    else
        down = latt[i, j-1]
    end
    return 2*latt[i,j]*(up + down + left + right)*J
end

function neighbors(latt, N, i, j)
    if i == 1
        left = latt[N, j]
    else
        left = latt[i-1,j]
    end
    if i == N
        right = latt[1, j]
    else
        right = latt[i+1, j]
    end
    if j == N
        up = latt[i, 1]
    else
        up = latt[i, j+1]
    end
    if j == 1
        down = latt[i, N]
    else
        down = latt[i, j-1]
    end
    return (up + down + left + right)
end

#Flip the spin
function FlipSpin!(latt, i, j)
    latt[i,j] = -latt[i,j]
end

# Do One MC Step
function MCmove!(L, N, iT, accept, M, E, J, Accum)
	for k in 1:N*N
		i = rand(1:N)
		j = rand(1:N)

		dE = DeltaE(L, N, i, j, J)  #energy change

		if dE <= 0 || rand() < exp(-dE*iT)
			accept += 1
			FlipSpin!(L, i, j)
			M += 2*L[i,j]
			E += dE
		end

	end
		Accum[1] += M
		Accum[2] += E
		return L, accept, M, E, Accum
end


function LatticePlot(L, N, accept, temp, k)
	heatmap(L, color = :grays, clims=(-1,1), legend=false, aspect_ratio = 1, title = temp,xlims=(0,N),xlabel=(string(accept)*" Flips, "*string(k)*" MC steps"))
	gui()
end

function Initialize(N, temp)
	accept = 0
	L = rand(-1:2:1, N,N)		#populate the lattice with random states of +1 or -1

	Magnetization = sum(L)
	Energy = 0

	for j in 1:N, i in 1:N
		Energy -= L[i,j] * neighbors(L, N, i, j)
	end

#	plt = heatmap(L, xlims=(0,N), aspect_ratio = 1, color = :grays, title = temp,xlabel="0 steps")
#	plot(plt)
#	gui()

	Accum = zeros(2)

	return L, accept, Magnetization, Energy, Accum
end
#-------------------- End Function Declarations --------------------#

#-------------------- Start Main --------------------#

function main(gif)
	#Initial Values
	J = 1.0
	mcsteps = 2000               # number of Monte Carlo Sweeps
	N = 256                         #Lattice size NxN
	steps = 10	                  #Number of steps per frame
	T =  1.0                        #Temperature
	e4 = exp(-4/T)                 #precomputed exponents for the metropolis algorithm dE = 2 or -2
	e8 = e4^2                      #precomputed exponents for the metropolis algorithm dE = 4 or -4
	temp = "T = "*string(T)         #string for showing temp on graph
	n = 1.0/(mcsteps*N)
	iT = 1.0/T

	L, accept, M, E, Accum = Initialize(N, temp)

	En = []
	Mg = []

	lay = @layout [[a; b] c]

if gif == true

	p = plot(plot([], ylabel="Energy"), plot([], ylabel="Magnetization"), plot(L,st=:heatmap,title=temp, legend=false, color=:grays, aspect_ratio=1), size=(1300,720), layout=lay)

	@gif for k in 1:mcsteps
		L, accept, M, E, Accum = MCmove!(L, N, iT, accept, M, E, J, Accum)
		push!(En, n*E)
		push!(Mg, n*M)

		plot!(p[1], En, legend=false, color=1)
		plot!(p[2], Mg, legend=false, color=2)
		p[3][1][:z] = L
	end every 10

	println(accept, " Spins are flipped!\n")

else

	for k in 1:mcsteps
		L, accept, M, E, Accum = MCmove!(L, N, iT, accept, M, E, J, Accum)
		LatticePlot(L, N, accept, temp, k)
	end

	println(accept, " Spings, are flipped!\n")
end

end

main(false)

println("\nDone!\n")

#-------------------- End Main --------------------#
