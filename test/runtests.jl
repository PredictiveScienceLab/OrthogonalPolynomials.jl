using OrthogonalPolynomials
using Base.Test

# Quadrature rules
q = QuadratureRule()
#quit()

#p = OrthpolBasis(lancz(q.x, q.w, 10)...)
p = OrthpolBasis(10, q)
print("α: ", p.α, "\n")
print("β: ", p.β, "\n")
print("γ: ", p.γ, "\n")
quit()
#p = OrthpolBasis([1., 1.], [0.2, 0.3], [1., 1.])

ϕ = Array(Float64, 2)

basis!(p, 0.4, ϕ)

print("ϕ: ", ϕ, "\n")

ϕ, ∇ϕ = grad_basis(p, 0.5)

print("∇ϕ: ", ∇ϕ, "\n")

#∂!(p, 0.2, ϕ, ∂ϕ)

#print("ϕ: ", ϕ, "\n")
#print("∂ϕ/∂x: ", ∂ϕ, "\n")

#print(∂(p, 0.2), "\n")


# write your own tests here
#@test 1 == 1
