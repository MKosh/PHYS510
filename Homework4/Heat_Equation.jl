#=-----------------------------------------------
Stationary Heat Equation

Using:
    - Gaussian Elimination
    - Inhomogeneous heat (diffusion) equation


-----------------------------------------------=#


using Plots

N = 10
L = 10
k = 1
Θ = -0.4
l = 1
T0 = 0.0
Tn = 2.0

#-------------------- Start Function Declarations --------------------#
function Elimination(A, v) # Gaussian Elimination Ax = v
	N = size(v,1)

	if A[1,1] == 0
		let j = 1
			while j<=N
				if A[j,1] == 0
					j +=1
				else
					rowTemp = A[1, :]
					A[1, :] = A[j, :]
					A[j, :] = rowTemp
					vtemp = v[1]
					v[1] = v[j]
					v[j] = vtemp
					break
				end
			end
		end
	end

	for m in 1:N
		div = A[m,m]
		A[m,:] /= div
		v[m] /= div

		for i in m+1:N
			mult = A[i,m]
			A[i,:] -= mult*A[m,:]
			v[i] -= mult*v[m]
		end
	end

	x = zeros(4,1)
	for m in N:-1:1
		x[m] = v[m]
		for i in m+1:N
			x[m] -= A[m,i]*x[i]
		end
	end

	return x
end

#-------------------- End Function Declarations --------------------#

#-------------------- Start Main Loop --------------------#

#-------------------- End Main Loop --------------------#

println("\nDone!\n")
