# dde_cobelli_debug.jl - Simulación con estimación de retardos

using DifferentialEquations, Plots, Optim, Random, Statistics

# 1. Definición de parámetros fisiológicos (conocidos)
p1, p2, p3, p4, p5 = 0.02, 0.03, 0.5, 0.1, 0.4
τ1_true, τ2_true = 10.0, 15.0
A, τd = 25 * 1000, 30.0

# 2. Funciones auxiliares
D_func(t) = t < 0 ? 0.0 : (A / τd) * t * exp(-t / τd)  # respuesta de ingesta

# 3. Historia inicial
G0, X0, I0 = 90.0, 0.0, 10.0
historia(t, _) = [G0, X0, I0]  # acepta segundo argumento requerido por HistoryFunction

# 4. Paquete de parámetros (con τ1, τ2 conocidos para simular datos)
params_true = (
    p1 = p1,
    p2 = p2,
    p3 = p3,
    p4 = p4,
    p5 = p5,
    τ1 = τ1_true,
    τ2 = τ2_true,
    D = D_func
)

# 5. Definición del sistema Cobelli
function cobelli!(du, u, h, p, t)
    G, X, I = u
    ε = 1e-8

    G_τ = try h(t - p.τ2 - ε, nothing)[1] catch; G end
    I_τ = try h(t - p.τ1 - ε, nothing)[3] catch; I end

    f(x) = x / (1 + x)
    σ(i) = i / (1 + i)
    φ(g) = g / (1 + g / 100)

    D_val = try p.D(t) catch; 0.0 end

    du[1] = -p.p1 * G - f(X) * G + D_val
    du[2] = -p.p2 * X + p.p3 * σ(I_τ)
    du[3] = -p.p4 * I + p.p5 * φ(G_τ)
end

# 6. Simular datos reales sintéticos
tspan = (0.0, 300.0)
lags = [τ1_true, τ2_true]
prob = DDEProblem(cobelli!, historia, tspan, params_true; constant_lags = lags)
sol = solve(prob, MethodOfSteps(Tsit5()), saveat=1.0)

# 7. Datos simulados con ruido realista
Random.seed!(123)
σG, σI = 5.0, 2.0  # desviación típica del error experimental
G_obs = sol[1, :] .+ randn(length(sol.t)) .* σG
I_obs = sol[3, :] .+ randn(length(sol.t)) .* σI

# 8. Función de pérdida: error cuadrático entre datos simulados y predichos
function loss(τ)
    τ1, τ2 = τ
    if τ1 <= 0 || τ2 <= 0
        return Inf
    end
    p_est = merge(params_true, (τ1 = τ1, τ2 = τ2))
    try
        prob_est = remake(prob; p=p_est, constant_lags=[τ1, τ2])
        sol_est = solve(prob_est, MethodOfSteps(Tsit5()), saveat=1.0, sensealg=InterpolatingAdjoint())
        G_pred = sol_est[1, :]
        I_pred = sol_est[3, :]
        return mean((G_obs .- G_pred).^2 .+ (I_obs .- I_pred).^2)
    catch
        return Inf
    end
end

# 9. Estimación de retardos con Optim
τ0 = [5.0, 5.0]  # inicial
result = optimize(loss, τ0, NelderMead(); iterations=100)
τ_est = Optim.minimizer(result)
println("\nEstimación de retardos:")
println("τ1 estimado = $(round(τ_est[1], digits=2)) (real: $τ1_true)")
println("τ2 estimado = $(round(τ_est[2], digits=2)) (real: $τ2_true)")

# 10. Visualización final
plot(sol.t, G_obs, label="Glucosa observada", lw=2)
plot!(sol.t, I_obs, label="Insulina observada", lw=2)
plot!(sol, vars=(1,), label="Glucosa verdadera", lw=2, ls=:dash)
plot!(sol, vars=(3,), label="Insulina verdadera", lw=2, ls=:dash)
plot!(xlabel="Tiempo (min)", ylabel="Concentración", title="Simulación y estimación de retardos", legend=:topright)
