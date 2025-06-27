using DifferentialEquations, Plots

# Parámetros
p₁, p₂, p₃ = 0.01, 0.02, 0.03
p₄, p₅ = 0.04, 0.05
τ₁, τ₂ = 10.0, 20.0

# Funciones no lineales
f(x) = 0.1 * x
σ(i) = i / (50 + i)
φ(g) = g^2 / (10000 + g^2)

# Historia: definida para t ∈ [-τ, 0]
function historia(t)
    [90.0, 0.0, 15.0]  # G(t), X(t), I(t)
end

# Modelo DDE con retardo doble
function modelo!(du, u, hist, p, t)
    G, X, I = u
    τ₁, τ₂ = p[6], p[7]
    
    Gτ = hist(t - τ₂)[1]  # G(t - τ₂)
    Iτ = hist(t - τ₁)[3]  # I(t - τ₁)

    du[1] = -p[1]*G - f(X)*G + 0       # D(t) = 0
    du[2] = -p[2]*X + p[3]*σ(Iτ)
    du[3] = -p[4]*I + p[5]*φ(Gτ)
end

# Paquete de parámetros (sin conflicto con historia)
p = (p₁, p₂, p₃, p₄, p₅, τ₁, τ₂)

# Declaración del problema con retardos constantes
lags = [τ₁, τ₂]
prob = DDEProblem(modelo!, historia, (0.0, 100.0), p; constant_lags=lags)

# Solución numérica
sol = solve(prob, MethodOfSteps(Tsit5()))

# Visualización
plot(sol, label=["G(t)" "X(t)" "I(t)"], lw=2, xlabel="Tiempo (min)", ylabel="Concentraciones", title="Modelo con Retardos y No Linealidades")
