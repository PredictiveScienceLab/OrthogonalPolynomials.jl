module OrthogonalPolynomials

export OrthogonalPolynomial

immutable type OrthogonalPolynomial{T<:Number}
    α::Vector{T}
    β::Vector{T}
    γ::Vector{T}


# Dimensionality of the input
input_dim(p::OrthogonalPolynomial) = 1


# The dimensionality of the output of a polynomial
output_dim(p::OrthogonalPolynomial) = 1

# The degree of the polyomial
degree(p::OrthogonalPolynomial) = length(p.α)

# Evaluate the polynomial at x
function eval{T}(p::OrthogonalPolynomial{T}, x::Number)
    R = promote_type(T, typeof(x))
    lenp = length(p)


end # module
