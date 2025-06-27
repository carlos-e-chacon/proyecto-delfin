using Plots

# Parámetros
a = 1.0
τ = 2.0
t0, tf, h = 0.0, 20.0, 0.01
times = t0:h:tf

# Condición inicial: historia en t < 0
history(t) = 1.0  # constante

# Vector solución
u = zeros(length(times))

# Rellenar con historia en [−τ, 0)
for (i, t) in enumerate(times)
    if t <= 0
        u[i] = history(t)
    end
end

# Función para buscar e interpolar u(t - τ)
function get_delayed_value(t_now, u, times)
    t_delay = t_now - τ
    if t_delay <= 0
        return history(t_delay)
    else
        # Interpolación lineal entre puntos
        idx = findlast(t -> t <= t_delay, times)
        if idx === nothing || idx == length(times)
            return u[idx]
        end
        t1, t2 = times[idx], times[idx + 1]
        u1, u2 = u[idx], u[idx + 1]
        return u1 + (u2 - u1) * (t_delay - t1) / (t2 - t1)
    end
end

# Simulación explícita de Euler
for i in 2:length(times)
    t = times[i - 1]
    u_delay = get_delayed_value(t, u, times)
    du = -a * u_delay
    u[i] = u[i - 1] + h * du
end

plot(times, u, xlabel="t", ylabel="u(t)", title="ED con retardo sin DDEProblem", lw=2)