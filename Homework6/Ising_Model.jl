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

#Initial Values
MCsteps = 20000000               #number of Monte Carlo Steps
N = 256                         #Lattice size NxN
Counts = zeros(N, 2)            #array to count the spin up and spin down states
steps = 1000                   #Number of steps per frame
T =  1.0                        #Temperature
ΔT = 0.1                       #Change in temperature
e4 = exp(-4/T)                 #precomputed exponents for the metropolis algorithm dE = 2 or -2
e8 = e4^2                      #precomputed exponents for the metropolis algorithm dE = 4 or -4
temp = "T = "*string(T)         #string for showing temp on graph
mag = energy = 0.0             #initialize magnetization and energy

#populate the lattice with random states of +1 or -1
L = rand(-1:2:1, N,N)
#L = ones(N,N)
#show the initial random lattice
plt = heatmap(L, xlims=(0,N), aspect_ratio = 1, color = :grays, title = temp,xlabel="0 steps", size=(1280,1280))
plot(plt)
gui()

#-------------------- Start Function Declarations --------------------#

#function for possible energy changes
function w(t)
    if t == 4 || t == -4
        return e8
    elseif t == 2 || t == -2
        return e4
    elseif t == 0
        return 0
    end
end

#function for calculating the initial energy
function CalcEnergy(latt, N)
	ener = 0.0
    for i in 1:N
	    for j in 1:N
		S = latt[i, j]
		    neigh = boundary(latt, N, i, j)
		    ener += -neigh*S
	    end
     end
     return ener/2.0
end

#function for counting spin up vs spin down

#Count spins
function spinCount(latt)
    global countUps = 0
    global countDowns = 0
    for i in 1:N
        for j in 1:N
            if latt[i, j] == 1
                countUps += 1
            end
            if latt[i, j] == -1
                countDowns += 1
            end
        end
    end
    return countUps, countDowns
end

#Function to sum the four nearest neighbors, accounting for toroidal boundaries
function boundary(latt, N, i, j)
    ΔE = 0
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
function FlipSpin(latt, i, j)
    latt[i,j] = -latt[i,j]
    return latt[i,j]
end
#-------------------- End Function Declarations --------------------#

#-------------------- Start Main --------------------#
let

	stepcount = 0
	#Calculate initial Magnetization and energy
	mag = sum(L)
	energy = CalcEnergy(L, N)
	Mcum = 0.0
	Ecum = 0.0

	#Monte Carlo Algorithm and graphing
	for k in 1:MCsteps

    	i = rand(1:N)    #pick a random x value to test
    	j = rand(1:N)    #pick a random y value to test

    	neigh = boundary(L, N, i, j)  #energy change

    	ΔE = 2*L[i,j]*neigh

    	if ΔE <= 0 || rand() < w(neigh)
        	stepcount += 1
        	FlipSpin(L, i, j)
        	energy += 2*L[i, j]*neigh
			mag -= 2*L[i, j]
        	if k%steps == 0    #graph the changes to the lattice every “steps” changes
            	heatmap(L, color = :grays, aspect_ratio = 1, title = temp,xlims=(0,N), size=(720,720),xlabel=(string(k)*" steps"))
				gui()
        	end
    	end
    #=Mcum += mag
    Mavg = Mcum/k
    Ecum += energy
    Eavg = Ecum/k
    data = [Mavg, Eavg]
    labels = ["Mavg", "Eavg"]
    display(plot(k, data, label = labels, color = [:red, :green])=#
	end
println("\nDone!\n", stepcount, " Spins are flipped!\n")
end

#-------------------- End Main --------------------#
