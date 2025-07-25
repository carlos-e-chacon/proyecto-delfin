using LinearAlgebra, Roots, PlotlyJS

function generar_grafica()
    τ1_vals = 10.0     # gran retardo de insulina
    τ2_vals = 15.0     # gran retardo en producción de insulina

    # Rango de retardos
    τ1_vals = 0.0:5.0:20.0
    τ2_vals = 0.0:5.0:20.0

    #Algoritmo detectó parámetros que generan inestabilidad
    p1 = 0.01027805603962313
    p2 = 0.010200427281483475
    p3 = 0.027081209428250527
    p4 = 0.010000921134245718
    p5 = 0.014363624547753803
    #τ1_vals = 1.521496870453525
    #τ2_vals = 2.970903119673876


    #Parámetro que tienden a generar inestabilidad (no exactos)
    p1 = 0.02     # captación de glucosa
    p2 = 0.03     # degradación de X
    p3 = 0.05     # MUY baja respuesta a insulina
    p4 = 0.2      # alta degradación de insulina
    p5 = 0.05     # muy baja secreción de insulina

    # Evaluar Re(λ) dominante para cada par (τ1, τ2)
    function raiz_dominante(τ1, τ2)
        A0 = [-p1 -1.0 0.0; 0.0 -p2 0.0; 0.0 0.0 -p4]
        A1 = [0.0 0.0 0.0; 0.0 0.0 p3; 0.0 0.0 0.0]
        A2 = [0.0 0.0 0.0; 0.0 0.0 0.0; p5 0.0 0.0]

        B0 = A0 + A1 + A2
        B1 = I + A1 * τ1 + A2 * τ2
        B2 = (A1 * τ1^2 + A2 * τ2^2) / 2

        P(λ) = det(-B0 + λ * B1 - λ^2 * B2)
        h = 1e-6
        dP(λ) = (P(λ + h) - P(λ - h)) / (2h)
        λ0 = -0.05 + 0.1im

        try
            λ = find_zero((P, dP), λ0, Roots.Newton())
            return real(λ)
        catch
            return NaN
        end
    end

    # Recolectar datos para 3D
    x, y, z = Float64[], Float64[], Float64[]
    for τ2 in τ2_vals, τ1 in τ1_vals
        push!(x, τ1)
        push!(y, τ2)
        push!(z, raiz_dominante(τ1, τ2))
    end

    push!(x, 1.521496870453525)
    push!(y, 2.970903119673876)
    push!(z, raiz_dominante(1.521496870453525, 2.970903119673876))

    # Visualización 3D
    trace = scatter3d(
        x = x, y = y, z = z,
        mode = "markers",
        marker = attr(
            size = 5,
            color = z,
            colorscale = "RdBu",
            colorbar = attr(title = "Re(λ)")
        )
    )

    layout = Layout(
        title = "Visualización 3D de Re(λ) dominante según τ₁ y τ₂",
        scene = attr(
            xaxis = attr(title = "τ₁"),
            yaxis = attr(title = "τ₂"),
            zaxis = attr(title = "Re(λ)")
        ),
        width = 800,
        height = 700
    )

    return plot(trace, layout)
end
grafica = generar_grafica()
grafica