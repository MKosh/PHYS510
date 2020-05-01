#=--------------------------------------------------#
The 2 Dimensional Ising Model

Graphs of the Energy, Magnetization, Specific Heat, and Magnetic Susceptibility
as a function of Temperature
Using:
	- Monte Carlo technique and the Metropolis algorithm
	- Default Julia RNG Mersenne Twister

---------------------------------------------------=#
using Plots


#-------------------- Start Function Declarations --------------------#

#function for calculating the energy
function CalcEnergy(latt, N, J, h)
	energy = 0.0
    for i in 1:N
	    for j in 1:N
		energy -= DeltaE(latt, N, i, j, J, h)
	    end
     end
     return energy
end

#Function to sum the four nearest neighbors, accounting for toroidal boundaries
function DeltaE(latt, N, i, j, J, h)
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
    return 2*J*latt[i,j]*(up + down + left + right) + 2*h*latt[i,j]
end


#Flip the spin
function FlipSpin!(latt, i, j)
    latt[i,j] = -latt[i,j]
end


function MCmove!(L, N, beta, J, h)
	for k in 1:N*N
		i = rand(1:N)
		j = rand(1:N)
		dE = DeltaE(L, N, i, j, J, h)  #energy change
		if dE <= 0 || rand() < exp(-dE*beta)
			FlipSpin!(L, i, j)
		end
	end
	return L
end

function Initialize(N)
	accept = 0
	L = rand(-1:2:1, N,N)		#populate the lattice with random states of +1 or -1
	return L, accept
end

function Plot(E, M, C, X, T, N)
	lay = @layout [a b; c d]
	Eplot = plot(T, E, legend=false, color=1, xlabel="Temperature", ylabel="<E>")
	scatter!(Eplot, T, E, color=1, title="N = "*string(N))

	Mplot = plot(T, abs.(M), legend=false, color=2, xlabel="Temperature", ylabel="<|M|>")
	scatter!(Mplot, T, abs.(M), color=2)#, ylims=(0,1))

	Cplot = plot(T, C, legend=false, color=3, xlabel="Temperature", ylabel="C")
	scatter!(Cplot, T, C, color=3)#, ylims=(0,2))

	Xplot = plot(T, X, legend=false, color=4, xlabel="Temperature", ylabel="X")
	scatter!(Xplot, T, X, color=4)

	plot(Eplot, Mplot, Cplot, Xplot, layout = lay, size=(1280,720))
	#savefig("Plots_n_50_h_01")
	gui()
end
#-------------------- End Function Declarations --------------------#

#-------------------- Start Main --------------------#

function main()
	#Initial Values
	J = 0.5
	h = 0.1
	eqsteps = 2000				# Number of steps to equilibrium
	mcsteps = 2000              # Number of Monte Carlo Steps
	N = 20                      # Lattice size NxN
	numTemps = 100				# Number of temp steps
	n1 = 1.0/(mcsteps*N^2)		# Normalize
	n2 = 1.0/(mcsteps^2*N^2)	# Normalize

	T = range(0.3,stop=3.0,length=numTemps)
	E = zeros(numTemps)
	M = zeros(numTemps)
	C = zeros(numTemps)
	X = zeros(numTemps)

	for t in 1:numTemps

		L, accept = Initialize(N)
		E1, M1, E2, M2 = 0.0, 0.0, 0.0, 0.0
		Temp = 1.0/T[t]
		Temp2 = Temp^2

		for i in 1:eqsteps
			L = MCmove!(L, N, Temp, J, h)
		end

		for k in 1:mcsteps
			L = MCmove!(L, N, Temp, J, h)

			Energy = CalcEnergy(L, N, J, h)
			Magnetization = sum(L)

			E1 += Energy
			M1 += Magnetization
			E2 += Energy^2
			M2 += Magnetization^2
		end

		E[t] = n1*E1
		M[t] = n1*M1
		C[t] = (n1*E2 - n2*E1^2)*Temp2
		X[t] = (n1*M2 - n2*M1^2)*Temp
	end

	Plot(E, M, C, X, T, N)
end

main()

println("\nDone!\n")

#-------------------- End Main --------------------#
