# Simulates heat transfer across a flat plate - Solution to the Laplace equation

# Import Packages
cd(@__DIR__)
using Pkg
Pkg.activate(".")

# Add package here
using Plots
using LinearAlgebra

# Initialise spacial domain
x = 0:0.01:1
y = 0:0.01:1
T = zeros((length(x),length(y)))
s = size(T)

# Boundary conditions
T[:, 1] .= 1 # T = 1 @ x = 0

a = zeros(length(y))
for i in range(1, length(y))
    a[i] = cos(6*3*y[i]*π/2) + 1
end
T[:, s[2]] = a

T[1, :] .= 1
T[s[1], :] = 1 .+ x

T[76, 11] = -2.5
T[51, 21] = 2.5
T[26, 26] = -0.5

# Calculate temperature distribution
# construct matrices for system
# AT = (D + L + U)T
# Tᵏ⁺¹ = -D⁻¹(L + U)Tᵏ
# Tᵏ⁺¹ = D⁻¹(b - (L + U)Tᵏ)
d = collect(-2 * ones(Int, length(x)))
dl = du = collect(ones(Int, length(x)-1))
#du = ones(length(x))
A = Tridiagonal(dl, d, du)
D = Diagonal(A)

# Let C = L + U, ∴ C = A - D
C = A - D

# - D⁻¹(L + U) = B
B = inv(D)*C

tolerance = 1e-6
err = 1
iterator = 1

while tolerance < err
    Tᵢ = T
    T₀ = T

    # global T = inv(D) * (T₀ - (C * Tᵢ))

    for i in range(2, convert(Int, length(x))-1)
        for j in range(2, convert(Int, length(y))-1)
            T[i,j] = 0.25(T₀[i+1,j] + T₀[i,j+1] + T₀[i-1,j] + T₀[i,j-1])
        end
    end

    # Set Boundary Conditions
    T[:, 1] .= 1 # T = 1 @ x = 0
    a = zeros(length(y))
    for i in range(1, length(y))
        a[i] = cos(6*3*y[i]*π/2) + 1
    end
    T[:, s[2]] = a
    T[1, :] .= 1
    T[s[1], :] = 1 .+ x

    T[76, 11] = -2.5
    T[51, 21] = 2.5
    T[26, 26] = -0.5

    global err = norm(T - T₀, Inf)
    global iterator += 1
end

println(iterator)
contour(T)




