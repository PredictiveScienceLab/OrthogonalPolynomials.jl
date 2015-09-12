export lancz!, lancz

# The method finds the recurrence formula for orthogonal polynomials
# using the Lanczos method. This is from:
#
#   W. J. Harrod, The numerically stable reconstruction of Jacobi matrices
#   from spectral data, Numer. Math. 44, 1984, 317-335
#
# The code is a modified version of the {s,d}lancz.f Fortran code from the
# orthpol package
#
# Arguments:
#   x   -   quadrature points  (ncap)
#   w   -   quadrature weights (ncap)
#   α   -   recurrence coefficient (n)
#   β   -   recurrence coefficient (n)
#   p0  -   work vector (ncap)
#   p1  -   work vector (ncap)
function lancz!{T<:Real}(x::Vector{T}, w::Vector{T}, α::Vector{T}, β::Vector{T},
                         p0::Vector{T}, p1::Vector{T})
    ncap = length(x)
    n = length(α)
    @assert ncap == length(w)
    @assert n == length(α)
    @assert n == length(β)
    @assert ncap == length(p0)
    @assert ncap == length(p1)
    @assert n > 0
    @assert n <= ncap
    copy!(p0, x)
    fill(p1, 0)
    for i=1:ncap-1
        π = w[i+1]
        γ = one(T)
        σ = zero(T)
        t = zero(T)
        λ = x[i+1]
        for k=1:i+1
            ρ = p1[k] + π
            tmp = γ * ρ
            tσ = σ
            if ρ <= 0
                γ = one(T)
                σ = zero(T)
            else
                γ = p1[k] / ρ
                σ = π / ρ
            end
            tk = σ * (p0[k] - λ) - γ * t
            p0[k] -= tk - t
            t = tk
            if σ <= 0
                π = tσ * p1[k]
            else
                π = (t * t) / σ
            end
            tσ = σ
            p1[k] = tmp
        end
    end
    copy!(α, 1, p0, 1, n)
    copy!(β, 1, p1, 1, n)
end
lancz!{T}(x::Vector{T}, w::Vector{T}, α::Vector{T}, β::Vector{T}) = 
    lancz!(x, w, α, β, Array(T, length(x)), Array(T, length(x)))
function lancz{T}(x::Vector{T}, w::Vector{T}, n::Integer)
    α = Array(T, n)
    β = Array(T, n)
    lancz!(x, w, α, β)
    return α, β
end
