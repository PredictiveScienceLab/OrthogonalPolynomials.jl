using OrthogonalPolynomials
using Base.Test

# Quadrature rules
q = QuadratureRule()
quit()

p = OrthPolBasis([1., 1.], [0.2, 0.3], [1., 1.])

print("degree: ", degree(p), "\n")

ϕ = zeros(2)
∂ϕ = zeros(2)

∂!(p, 0.2, ϕ, ∂ϕ)

print("ϕ: ", ϕ, "\n")
print("∂ϕ/∂x: ", ∂ϕ, "\n")

print(∂(p, 0.2), "\n")


# write your own tests here
#@test 1 == 1
