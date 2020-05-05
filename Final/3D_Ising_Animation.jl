#=--------------------------------------------------#
The 2 Dimensional Ising Model

Plots 3D visualization of the 3D Ising model using Makie
Various "volume" algorithms can be used (:iso, :absorption, :mip) in the
visualization to varying degrees of success. Beyond the critical temperature
there isn't much to see because it all becomes a big mess.
Using:
	- Monte Carlo technique and the Metropolis algorithm
	- Default Julia RNG Mersenne Twister
	- Makie Plotting package

---------------------------------------------------=#
using Makie	# Makie has more functionality than Plots (especially with 3D or interactive plots), but is also more complicated

#-------------------- Start Function Declarations --------------------#
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
	L = rand(-1:2:1, N,N,N)
	return L
end

#-------------------- End Function Declarations --------------------#

#-------------------- Start Main --------------------#

function main()
	#Initial Values
	J = 0.5
	h = 0.0
	Tc = 2.21
	eqsteps = 100				# Number of steps to equilibrium
	mcsteps = 500             # Number of Monte Carlo Steps
	N = 64                      # Lattice size NxN
	numTemps = 20				# Number of temp steps

	T = range(2.2,stop=1.0,length=numTemps)
	L = Initialize(N)
	Temp = 1.0/T[1]

	for i in 1:eqsteps
		L = MCmove!(L, N, Temp, J, h)
	end

	vis = volume(L, transparency=true, algorithm=:absorption, resolution=(1280,720))

	for t in 1:numTemps

		Temp = 1.0/T[t]

		record(vis, "3DIsing"*string(t)*".mp4", framerate=24, compression=0, 1:mcsteps) do i
			L = MCmove!(L, N, Temp, J, h)
			vis.plots[2][:volume] = L
			rotate_cam!(vis, 0.01,0.0,0.0)
		end
	end
end

main()

println("\nDone!\n")

#-------------------- End Main --------------------#
