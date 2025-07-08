using LinearAlgebra

# Parámetros del sistema (puedes ajustarlos)
p1, p2, p3, p4, p5 = 0.02, 0.03, 0.5, 0.1, 0.4
τ1, τ2 = 5.0, 15.0

# Matrices del modelo linealizado
A0 = [-p1 -1.0 0.0; 0.0 -p2 0.0; 0.0 0.0 -p4]
A1 = [0.0 0.0 0.0; 0.0 0.0 p3; 0.0 0.0 0.0]
A2 = [0.0 0.0 0.0; 0.0 0.0 0.0; p5 0.0 0.0]

# Construcción de matrices con expansión de Taylor de orden 2
B0 = A0 + A1 + A2
B1 = I + A1 * τ1 + A2 * τ2
B2 = (A1 * τ1^2 + A2 * τ2^2) / 2

# Determinante característico aproximado como función de λ
P(λ) = det(-B0 + λ * B1 - λ^2 * B2)

using Roots

λ0 = -0.01 + 0.3im  # Punto inicial (elige según tu heatmap o intuición)

# Derivada
h = 1e-6
dP(λ) = (-P(λ + 2h) + 8P(λ + h) - 8P(λ - h) + P(λ - 2h)) / (12h)

# Newton-Raphson
λ_raiz = find_zero((P, dP), λ0, Roots.Newton())
println("Raíz encontrada: λ = $λ_raiz")