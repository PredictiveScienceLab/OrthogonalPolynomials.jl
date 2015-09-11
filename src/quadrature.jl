export QuadratureRule, fejer


import Base: length, one, zero, convert


immutable type QuadratureRule{T<:Real}
    x::Vector{T}
    w::Vector{T}

    function QuadratureRule(x::Vector{T}, w::Vector{T})
        @assert length(x) == length(w)
        for i = 1:length(w)
            @assert w[i] > 0
        end
        new(x, w)
    end
end
QuadratureRule{T<:Real}(x::Vector{T}) = QuadratureRule{T}(x,
                                            ones(T, length(x)) / length(x))
function QuadratureRule{T<:Real}(wf::Function, left::T=-1., right::T=1.,
                                 n::Integer=500)
    @assert left < right
    x, w = fejer(T, n)
    if left == -Inf && right == Inf
        ϕ, dϕ = symtr(x) 
    elseif right == Inf
        ϕ, dϕ = tr(x)
        ϕ += left
    elseif left == -Inf
        ϕ, dϕ = tr(-x)
        ϕ = right - ϕ
    else
        ϕ = convert(T, 0.5) * ((right - left) * x + right + left)
        dϕ = 0.5 * (right - left)
    end
    μ = w .* map(wf, x) .* dϕ
    return QuadratureRule{T}(ϕ, μ)
end
QuadratureRule{T<:Real}(left::T=-1., right::T=1., n::Integer=500) = 
                       QuadratureRule(identity, left, right, n)

length(q::QuadratureRule) = length(q.x)


function integrate{T<:Real}(f::Function, q::QuadratureRule{T})
    s = convert(T, 0)
    for i = 1:length(q)
        s += f(q.x[i]) * q.w[i]
    end
    return s
end


# Implements a transformation from [-1, 1] ⟶   [-∞, ∞]
function symtr{T<:Real}(t::T)
    t2 = t * t
    dϕ = 1 - t2
    ϕ = t / dϕ
    dϕ *= dϕ
    dϕ = (t2 + 1) / dϕ
    return ϕ, dϕ
end


@vectorize_1arg Real symtr


# Implements a transformation of [-1, 1] ⟶  [0, ∞]
function tr{T<:Real}(t::T)
    dϕ = 1 - t
    ϕ = (1 + t) / dϕ
    dϕ *= dϕ
    dϕ = 2 / dϕ
    return ϕ, dϕ
end


@vectorize_1arg Real tr


# This is the fejer quadrature rule
function fejer!{T<:Real}(n::Integer, x::Vector{T}, w::Vector{T})
    one = convert(T, 1)
    half = convert(T, 0.5)
    dpi = 4 * atan(one)
    nh = n / 2
    np1h = (n + 1) / 2
    dn = convert(T, n)
    for k = 1:nh
        x[n+1-k] = cos(half * convert(T, 2 * k - 1) * dpi / dn)
        x[k] = x[n+1-k]
    end
    if 2 * nh != n
        x[np1h] = 0
    end
    for k = 1:np1h
        dc1 = one
        dc0 = 2 * x[k] * x[k] - 1
        dt = 2 * dc0
        dsum = dc0 / 3
        for m = 2:nh
            dc2 = dc1
            dc1 = dc0
            dc0 = dt * dc1 - dc2
            dsum += dc0 / (4 * m * m - 1)
        end
        w[k] = 2 * (1 - 2 * dsum) / dn
        w[n+1-k] = w[k]
    end
end


function fejer(T, n::Integer)
    x = Array(T, n)
    w = Array(T, n)
    fejer!(n, x, w)
    return x, w
end


fejer(n::Integer) = fejer(Float64, n)
