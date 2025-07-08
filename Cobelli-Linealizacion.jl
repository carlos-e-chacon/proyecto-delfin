using LinearAlgebra, PlotlyJS

# Parámetros del modelo
p1, p2, p3, p4, p5 = 0.02, 0.03, 0.5, 0.1, 0.4
τ1, τ2 = 5.0, 15.0

A0 = [-p1 -1.0 0.0; 0.0 -p2 0.0; 0.0 0.0 -p4]
A1 = [0.0 0.0 0.0; 0.0 0.0 p3; 0.0 0.0 0.0]
A2 = [0.0 0.0 0.0; 0.0 0.0 0.0; p5 0.0 0.0]

function M(λ)
    B0 = A0 + A1 + A2
    B1 = I + A1 * τ1 + A2 * τ2
    B2 = (A1 * τ1^2 + A2 * τ2^2) / 2
    return -B0 + λ * B1 - λ^2 * B2
end

P(λ) = det(M(λ))

# Rango del plano complejo
λ_re = LinRange(-0.1, 0.1, 200)
λ_im = LinRange(-1.5, 1.5, 200)

Z = [log10(abs(P(x + im*y))) for y in λ_im, x in λ_re]

trace = heatmap(
    x = λ_re,
    y = λ_im,
    z = Z,
    colorscale = "Viridis",
    colorbar = attr(title="log₁₀|P(λ)|")
)

layout = Layout(
    title = "Plano complejo – Magnitud logarítmica del determinante",
    xaxis_title = "Re(λ)",
    yaxis_title = "Im(λ)",
    width = 800,
    height = 700
)

plot([trace], layout)