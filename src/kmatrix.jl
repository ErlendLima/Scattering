using ArgCheck
using LinearAlgebra
using DocStringExtensions
using .Scattering

"""
Contains parameters for running K-Matrix calculations.

Takes as input the mesh size, equivalent to the number of
iteration steps. The time required to solve KMatrix(N) goes as 𝒪(N²).
"""
struct KMatrix <: Method
    N::Int  # Number of iteration steps
end

"""
Compute the phase shift of S-waves by computing the reactance matrix.

$(SIGNATURES)
where `k₀` is the momentum for which to compute, `m` is the mass, and `P` the potential.

# Details

The goal is to obtain the S-wave phase shift of the potential P for the momentum `k`.
This is done by solving the matrix equation
             `K = A⁻¹V`
where K is the reactance matrix. V is the matrix of the potential evaluated at the mesh points,
and A describes the Lippmann-Schwinger equation for K. The phase shift is then related to
the last diagonal element of K as
             `K(k₀, k₀) = 1/mk₀ ⋅ tan(δₗ₌₀(k₀))`
"""
function (method::KMatrix)(k₀, m, P::Potential)
    A, V = createA(k₀, m, method.N, P)
    # Create the reactance matrix K using the matrix equation AK = V
    K = inv(A)*V
    @assert size(A) == size(V) == size(K) == (method.N+1, method.N+1)

    # Diagonal elements of K are related to the phase shifts δₗ
    # as K(k₀, k₀) = -tan(δₗ)/mk₀
    Kk0 = K[end, end]
    δᵢ = atan(-Kk0*m*k₀) / pi * 180
    #energy = k₀^2 / m * 197
    δᵢ
end

"""
Construct the A matrix involved in K-Matrix calculations.

$(SIGNATURES)

where `k₀` is the momentum for which to compute, `m` is the mass, `N` the
mesh size and `P` the potential.

# Details

The goal is the set up the matrices V and A of the matrix equation

    AK = V

In the first phase, the mesh points are constructed by Gaussian quadrature,
and V is constructed as the potential evaluated at P(kᵢ, kⱼ) ∀ i,j in the mesh.
The set up of A is more complicated, and is delegated to another function.
"""
function createA(k₀, m, N::Int, P::Potential)
    k, ω = QuadGK.gauss(N)
    transform!(k, ω)
    @assert k₀ ∉ k "k₀ can not be in the mesh points k"
    push!(k, k₀)
    V = [P(ki, kj) for ki in k, kj in k]
    A = createA(V, k, ω, m)
    return A, V
end

"""
Construct the A matrix involved in K-Matrix calculations.

$(SIGNATURES)

where `V` matrix of the potential evaluated at the mesh points, `k` is the momenta,
`ω` are the quadrature weights, and `m` is the mass.
# Details

The goal is the set up the matrices V and A of the matrix equation

    K = VA.

Here A is constructed as described in
https://manybodyphysics.github.io/NuclearForces/doc/pub/scatteringtheory/html/scatteringtheory.html
"""
function createA(V, k, ω, m)
    @argcheck length(k) - 1 == length(ω)
    @argcheck size(V, 1) == size(V, 2) == length(k)
    @argcheck m > 0

    N = length(k)-1
    k₀ = k[end]
    A = diagm(0 => repeat([1.0,], N+1))  # δᵢⱼ
    uⱼ = 0.0
    @inbounds for j in 1:N+1
        uⱼ = 0.0
        if j ≠ N+1
            uⱼ = 2/π * ω[j]*k[j]^2/((k₀^2 - k[j]^2)/m)
        else
            for n in 1:N
                uⱼ += ω[n]*k₀^2/((k₀^2 - k[n]^2)/m)
            end
            uⱼ *= -2/π
        end

        for i in 1:N+1
            A[i, j] -= V[i, j] * uⱼ
        end
    end
    return A
end

"""
Map the Gaussian quadrature points from [0, 1] to (-∞, ∞)

$(SIGNATURES)

"""
function transform!(xs, weights)
    @. weights = π/4 * weights / cos(π/4*(1.0 + xs))^2
    @. xs      = tan(π/4*(1.0 + xs))
end

