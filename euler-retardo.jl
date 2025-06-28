using DifferentialEquations
using PlotlyJS
using DataFrames
# 1. Definimos la ecuación con retardo
function f!(du, u, h, p, t)
    τ = 0.5
    du[1] = -2u[1] + h(p, t - τ)[1]
end

# 2. Condición inicial como función constante
historia(p, t) = [1.0]

# 3. Intervalo de integración y problema
tspan = (0.0, 10.0)
prob = DDEProblem(f!, historia, tspan, nothing; constant_lags=[1.0])

# 4. Solución numérica
sol = solve(prob, MethodOfSteps(Tsit5()), saveat=0.05)

# 5. Extraer solución y valores del retardo
ts = sol.t
ys = sol.(ts)         # y(t)
ys_tau = sol.(ts .- 1.0)  # y(t - τ)

# 7. Extraer y(t) y y(t - τ)
yt = [y[1] for y in ys]
ytau = [y[1] for y in ys_tau]

# 8. Crear gráfico en PlotlyJS
trace = PlotlyJS.scatter(
    x = ts,
    y = yt,
    mode = "lines+markers",
    marker = attr(
        color = ytau,
        colorscale = "Viridis",
        colorbar = attr(title = "y(t - τ)"),
        size = 8,
        line = attr(width = 0)
    ),
    line = attr(color = "rgba(0,0,0,0.3)"),
    name = "y(t)"
)

layout = Layout(
    title = "Solución de DDE con color según y(t - τ)",
    xaxis = attr(title = "t"),
    yaxis = attr(title = "y(t)"),
    plot_bgcolor = "white"
)

# 9. Mostrar gráfico
PlotlyJS.plot(trace, layout)

#-----------------------

using DifferentialEquations
using PlotlyJS
using PlotlyBase
historia(p,t) = [1.0]
tspan = (0.0, 10.0)
save_times = 0.0:0.05:10.0
taus = 0.1:0.1:2.0

frames = []
for τ in taus
    function f!(du,u,h,p,t)
        du[1] = -2u[1] + h(p,t - τ)[1]
    end
    prob = DDEProblem(f!, historia, tspan, nothing; constant_lags=[τ])
    sol = solve(prob, MethodOfSteps(Tsit5()), saveat=save_times)

    ts = sol.t
    ys = sol.(ts)
    ys_tau = sol.(ts .- τ)
    yt = [y[1] for y in ys]
    ytau = [y[1] for y in ys_tau]

    trace = PlotlyJS.scatter(
        x=ts,
        y=yt,
        mode="lines+markers",
        marker=attr(
            color=ytau,
            colorscale="Viridis",
            size=6,
            colorbar=nothing
        ),
        line=attr(color="rgba(0,0,0,0.3)")
    )
    push!(frames, (data=[trace], name=string(round(τ,digits=2))))
end

layout = Layout(
    title="Animación de y(t) con color según y(t - τ)",
    xaxis=attr(title="t"),
    yaxis=attr(title="y(t)"),
    updatemenus=[attr(
        type="buttons",
        showactive=false,
        y=1.1,
        x=1.05,
        buttons=[attr(
            label="Play",
            method="animate",
            args=[nothing, attr(frame=attr(duration=200, redraw=true), fromcurrent=true, mode="immediate")]
        )]
    )],
    sliders=[attr(
        steps=[attr(
            method="animate",
            label=string(round(τ,digits=2)),
            args=[[string(round(τ,digits=2))], attr(mode="immediate", frame=attr(duration=0, redraw=true), transition=attr(duration=0))]
        ) for τ in taus],
        transition=attr(duration=0),
        x=0.1,
        xanchor="left",
        y=-0.3,
        yanchor="top",
        currentvalue=attr(font=attr(size=16), prefix="τ = ", visible=true),
        len=0.9
    )]
)

PlotlyJS.plot([frames[1].data[1]], layout, frames)