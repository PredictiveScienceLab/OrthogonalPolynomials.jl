module OrthogonalPolynomials

export OrthPolBasis, input_dim, output_dim, degree, polyeval!, polyeval,
       ∂!, ∂

include("quadrature.jl")

immutable type OrthPolBasis{T<:Real}
    α::Vector{T}
    β::Vector{T}
    γ::Vector{T}
end


# Dimensionality of the input
input_dim(p::OrthPolBasis) = 1

# The dimensionality of the output of a polynomial
output_dim(p::OrthPolBasis) = length(p.α)

# The degree of the polyomial
degree(p::OrthPolBasis) = length(p.α) - 1


# Evaluate the polynomial at x and write the result at ϕ
function polyeval!{T}(p::OrthPolBasis{T}, x::T, ϕ::Vector{T})
    ρ = degree(p)
    @assert output_dim(p) == length(ϕ)
    ϕ[1] = 1 / p.γ[1]
    if ρ >= 1
        ϕ[2] = (x - p.α[1]) * (ϕ[1] / p.γ[2])
    end
    for i = 3:ρ
        ϕ[i] = ((x - p.α[i-1]) * ϕ[i-1] - p.β[i-1] * ϕ[i-2]) / p.γ[i]
    end
end


# Evaluate the polynomial at x (N vector) and write the result
# to Φ (output_dim matrix x N)
function polyeval!{T}(p::OrthPolBasis{T}, x::Vector{T}, Φ::Matrix{T})
    num_points = length(X)
    @assert output_dim(p) == size(Φ, 1)
    @assert num_points == size(Φ, 2)
    for i = 1:num_points
        polyeval!(p, x[i], Φ[:, i])
    end
end


# Evaluate the polynomial at x and return the result as a vector
function polyeval{T}(p::OrthPolBasis{T}, x::T)
    ϕ = Array(T, output_dim(p))
    polyeval!(p, x, ϕ)
    return ϕ
end


# Evaluate the gradient of a polynomial with respect to x
function ∂!{T}(p::OrthPolBasis{T}, x::T, ϕ, ∂ϕ)
    polyeval!(p, x, ϕ)
    ρ = degree(p)
    @assert output_dim(p) == length(ϕ)
    ∂ϕ[1] = 0.
    if ρ >= 1
        ∂ϕ[2] = ϕ[1] / p.γ[1]
    end
    for i = 3:ρ
        ∂ϕ[i] = (ϕ[i-1] + (x - p.α[i-1]) * ∂ϕ[i-1] - p.β[i-1] * ∂ϕ[i-2]) / p.γ[i]
    end
end


# Evaluate the gradient of a polynomial and return the result as a vector
function ∂{T}(p::OrthPolBasis{T}, x::T)
    m = output_dim(p)
    ϕ = Array(T, m)
    ∂ϕ = Array(T, m)
    ∂!(p, x, ϕ, ∂ϕ)
    return ϕ, ∂ϕ
end

end # module
