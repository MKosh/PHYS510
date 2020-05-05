#=--------------------------------------------------#
The 2 Dimensional Ising Model

Graphs of the Energy, Magnetization, Specific Heat, and Magnetic Susceptibility
as a function of Temperature OR the log plots of the thermal properties vs reduced
temperature and their slopes for finding the critical exponents
Using:
	- Monte Carlo technique and the Metropolis algorithm
	- Default Julia RNG Mersenne Twister

---------------------------------------------------=#
using Plots

#-------------------- Start Function Declarations --------------------#
function CalcEnergy(latt, N, J, h)	# Calculating the energy in the lattice
	energy = 0.0
    @inbounds for i in 1:N, j in 1:N, k in 1:N
		energy -= DeltaE(latt, N, i, j, k, J, h)
     end
     return energy
end

function DeltaE(latt, N, i, j, k, J, h)	# Sum the 6 nearest neighbors, accounting for toroidal boundaries
    if i == 1
        left = latt[N, j, k]
    else
        left = latt[i-1, j, k]
    end
    if i == N
        right = latt[1, j, k]
    else
        right = latt[i+1, j, k]
    end
    if j == N
        up = latt[i, 1, k]
    else
        up = latt[i, j+1, k]
    end
    if j == 1
        down = latt[i, N, k]
    else
        down = latt[i, j-1, k]
    end
	if k == 1
		backward = latt[i, j, N]
	else
		backward = latt[i, j, k-1]
	end
	if k == N
		forward = latt[i, j, 1]
	else
		forward = latt[i, j, k+1]
	end
    return 2*J*latt[i,j,k]*(up + down + left + right + forward + backward) + 2*h*latt[i,j,k]
end

function FlipSpin!(latt, i, j, k)	# Flip the spin, pretty self explanatory
    latt[i, j, k] = -latt[i, j, k]
end

function MCmove!(L, N, beta, J, h)	# Iterate over all lattice sites and determine if the spin should flip
	@inbounds for iter in 1:N^3
		i = rand(1:N)
		j = rand(1:N)
		k = rand(1:N)
		dE = DeltaE(L, N, i, j, k, J, h)  #energy change
		if dE <= 0 || rand() < exp(-dE*beta)
			FlipSpin!(L, i, j, k)
		end
	end
	return L
end

function Initialize(N)	# Populate the lattice with random states of +1 or -1
	accept = 0
	L = rand(-1:2:1, N,N,N)
	return L, accept
end

function LinReg_Bias(X, Y) # Linear regression to find the slope of the plots
	n_obs = size(X, 1)
	X = hcat(ones(Float64, n_obs), X)	# hcat: horizontal concatination
	β = inv(X' * X) * X' * Y	# inv: inverse, X': transpose
	return β
end

function Plot(E, M, C, X, T, N)	# Plot the thermal properties
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
	#savefig("Plots of Thermal Properties")
	gui()
end

function Plot_Regs(E, M, C, X, tau, N)	# Plot the log-log plots and their slopes
	lnτ = log.(tau)
	lnM = log.(abs.(M))
	lnC = log.(C)
	lnX = log.(X)

	β = LinReg_Bias(lnτ, lnM)
	α = LinReg_Bias(lnτ, lnC)
	γ = LinReg_Bias(lnτ, lnX)

	beta_plot = scatter(lnτ, lnM, color=1, title=("N = "*string(N)) ,label="", xlabel="ln(|T - Tc|)", ylabel="ln(M)")
	plot!(beta_plot, x -> β[1] + β[2]*x, color=2, label=(string(round(β[2], digits=4))))

	alpha_plot = scatter(lnτ, lnC, color=3, label="", xlabel="ln(|T - Tc|)", ylabel="ln(C)")
	plot!(alpha_plot, x -> α[1] + α[2]*x, color=4, label=(string(round(α[2], digits=4))))

	gamma_plot = scatter(lnτ, lnX, color=5, label="", xlabel="ln(|T - Tc|)", ylabel="ln(X)")
	plot!(gamma_plot, x -> γ[1] + γ[2]*x, color=6, label=(string(round(γ[2], digits=4))))

	values_plot = scatter(framestyle=:none, [0.5, 0.5, 0.5], [3, 2, 1], markeralpha=0, legend=false,xlims=(0,1), ylims=(0,4),
		series_annotations=[Plots.text("beta = "*string(round(abs(β[2]),digits=4))*",  Accepted = 0.3264"),
		Plots.text("alpha = "*string(round(abs(α[2]),digits=4))*",  Accepted = 0.1101"), Plots.text("gamma = "*string(round(abs(γ[2]),digits=4))*",  Accepted = 1.2372")])

	#plot(gamma_plot, size=(1280, 720))
	plot(beta_plot, alpha_plot, gamma_plot, values_plot, size=(1280,720))
	#savefig("Exponents")
	gui()
end
#-------------------- End Function Declarations --------------------#

#-------------------- Start Main --------------------#

function main(plt)
	#Initial Values
	J = 0.5
	h = 0.0
	Tc = 2.21
	eqsteps = 1000				# Number of steps to equilibrium
	mcsteps = 2000            	# Number of Monte Carlo Steps
	N = 10                      # Lattice size NxN
	numTemps = 150				# Number of temp steps
	n1 = 1.0/(mcsteps*N^3)		# Normalize
	n2 = 1.0/(mcsteps^2*N^3)	# Normalize

	T = range(2.0,stop=2.15,length=numTemps)
	τ = zeros(numTemps)
	E = zeros(numTemps)
	M = zeros(numTemps)
	C = zeros(numTemps)
	X = zeros(numTemps)

	@progress for t in 1:numTemps	# @progress gives a progress bar in Atom/Juno

		L, accept = Initialize(N)
		E1, M1, E2, M2 = 0.0, 0.0, 0.0, 0.0
		τ[t] = abs(T[t]-Tc)
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

	if plt == true
		Plot(E, M, C, X, T, N)
	else
		Plot_Regs(E, M, C, X, τ, N)
	end
end

main(false)

println("\nDone!\n")

#-------------------- End Main --------------------#
