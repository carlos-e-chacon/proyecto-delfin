# dde_cobelli_debug.jl - Simulación con chequeos de parámetros

using DifferentialEquations, Plots

# 1. Definición de parámetros fisiológicos
p1, p2, p3, p4, p5 = 0.02, 0.03, 0.5, 0.1, 0.4
τ1, τ2 = 10.0, 15.0
#Tortilla
#A, τd = 25 * 1000, 30.0

#Coca-Cola 600 ml
A, τd = 63 * 1000, 15.0

# 2. Funciones auxiliares
D_func(t) = t < 0 ? 0.0 : (A / τd) * t * exp(-t / τd)   # asegúrate de paréntesis para evitar precedencia

# 3. Historia inicial
G0, X0, I0 = 90.0, 0.0, 10.0
historia(t, _) = [G0, X0, I0]  # acepta segundo argumento requerido por HistoryFunction

# 4. Paquete de parámetros con chequeo directo
params = (
    p1 = p1,
    p2 = p2,
    p3 = p3,
    p4 = p4,
    p5 = p5,
    τ1 = τ1,
    τ2 = τ2,
    D = D_func
)

# 5. Definición del sistema Cobelli con chequeo de errores
function cobelli!(du, u, h, p, t)
    G, X, I = u
    ε = 1e-8  # Evitar evaluación justo en el borde de t = 0

    G_τ = G

    I_τ = I

    f(x) = x / (1 + x)
    σ(i) = i / (1 + i)
    φ(g) = g / (1 + g / 100)

    D_val = p.D(t)
    

    du[1] = -p.p1 * G - f(X) * G + D_val
    du[2] = -p.p2 * X + p.p3 * σ(I_τ)
    du[3] = -p.p4 * I + p.p5 * φ(G_τ)
end

# 6. Tiempo y retardos
lags = [params.τ1, params.τ2]
tspan = (0.0, 300.0)

# 7. Construcción del problema
prob = DDEProblem(cobelli!, historia, tspan, params; constant_lags = lags)

# 8. Resolución del problema
sol = solve(prob, MethodOfSteps(Tsit5()))

# 9. Visualización de tortilla
#plot(sol, xlabel="Tiempo (min)", ylabel="Concentración", lw=2, label=["Glucosa G(t)" "X(t)" "Insulina I(t)"], legend=:topright)

#Visualización de Coca-Cola 600 ml
plot(t -> D_func(t)/1000, 0, 60, xlabel="Tiempo (min)", ylabel="D(t)", lw=2, title="Entrada de glucosa por Coca-Cola")

#plot(sol, xlabel="Tiempo (min)", ylabel="Concentración", lw=2, label=["Glucosa G(t)" "X(t)" "Insulina I(t)"], legend=:topright)
