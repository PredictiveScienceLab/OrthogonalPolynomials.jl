export OrthpolBasis, input_dim, dim, degree, basis!, basis,
       grad_basis!, grad_basis

# A basis of orthonormal polynomials defined on the real numbers
# The polynomials are defined using the recurrence relation:
# ϕ[i](x) = ((x - α[i])ϕ[i-1] - β[i-1]ϕ[i-2]) / γ[i]
immutable OrthpolBasis{T<:Real}
    α::Vector{T}
    β::Vector{T}
    γ::Vector{T}
end
OrthpolBasis{T}(α::Vector{T}, β::Vector{T}) =
    OrthpolBasis{T}(α, β, ones(T, length(α)))
OrthpolBasis{T}(degree::Integer, q::QuadratureRule{T}) =
    OrthpolBasis(lancz(q.x, q.w, degree + 1)...)
OrthpolBasis{T}(degree::Integer,
                wf::Function, left::T, right::T, ncap::Integer)
    = OrthpolBasis(degree, QuadratureRule(wf, left, right, ncap))
#OrthpolBasis(degree::Integer,
#             wf::Function, left::Real=-1, right::Real=1,; ncap::Integer=500)
#    = OrthpolBasis(Float64, degree, wf, left, right; ncap=ncap)
#OrthpolBasis(T::DataType, degree::Integer) =
#    OrthpolBasis{T}(degree, identity)
#OrthpolBasis(degree::Integer) = OrthpolBasis(Float64, degree)

# Dimensionality of the input
input_dim(p::OrthpolBasis) = 1

# The dimension of the vector space spanned by the polynomial
dim(p::OrthpolBasis) = length(p.α)

# The degree of the polyomial
degree(p::OrthpolBasis) = dim(p) - 1


# Evaluate the orthogonal polynomial basis at x and write the result at ϕ
function basis!{T}(p::OrthpolBasis{T}, x::T, ϕ::Vector{T})
    α = p.α
    β = p.β
    γ = p.γ
    ϕ[1] = one(T) / γ[1]
    if dim(p) >= 1
        ϕ[2] = (x - α[1]) * ϕ[1] / γ[2]
    end
    for i=3:dim(p)
        ϕ[i] = ((x - α[i]) * ϕ[i-1] - β[i-1] * ϕ[i-2]) / γ[i]
    end
end


# Evaluate the polynomial at x (N vector) and write the result
# to Φ (output_dim matrix x N)
function basis!{T}(p::OrthpolBasis{T}, x::Vector{T}, Φ::Matrix{T})
    num_points = length(X)
    @assert output_dim(p) == size(Φ, 1)
    @assert num_points == size(Φ, 2)
    for i = 1:num_points
        basis!(p, x[i], Φ[:, i])
    end
end


# Evaluate the polynomial at x and return the result as a vector
function basis{T}(p::OrthpolBasis{T}, x::T)
    ϕ = Array(T, output_dim(p))
    basis!(p, x, ϕ)
    return ϕ
end


# Evaluate the gradient of a polynomial with respect to x
function grad_basis!{T}(p::OrthpolBasis{T}, x::T, ϕ, ∇ϕ)
    basis!(p, x, ϕ)
    α = p.α
    β = p.β
    γ = p.γ
    ∇ϕ[1] = 0.
    if dim(p) >= 1
        ∇ϕ[2] = ϕ[1] / γ[1]
    end
    for i = 3:dim(p)
        ∇ϕ[i] = (ϕ[i-1] + (x - α[i-1]) * ∇ϕ[i-1] - β[i-1] * ∇ϕ[i-2]) / γ[i]
    end
end


# Evaluate the gradient of a polynomial and return the result as a vector
function grad_basis{T}(p::OrthpolBasis{T}, x::T)
    ϕ = Array(T, dim(p))
    ∇ϕ = Array(T, dim(p))
    grad_basis!(p, x, ϕ, ∇ϕ)
    return ϕ, ∇ϕ
end
