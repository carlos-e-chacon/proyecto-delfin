using LinearAlgebra
using BlackBoxOptim
using Polynomials

# (1) Función para construir el polinomio característico
function matriz_caracteristica(λ, p1, p2, p3, p4, p5, τ1, τ2)
    # Aproximación de los retardos: e^{-λτ} ≈ 1 - λτ + (λτ)^2/2
    e1 = 1 - λ*τ1 + (λ*τ1)^2/2
    e2 = 1 - λ*τ2 + (λ*τ2)^2/2

    M = [-p1 - λ        -1          0;
          0         -p2 - λ     p3 * e1;
         p5 * e2        0       -p4 - λ]

    return det(M)
end

# (2) Encuentra la raíz dominante del polinomio característico
function raiz_dominante(p1, p2, p3, p4, p5, τ1, τ2)
    # Define el polinomio como función de λ
    function P(λ)
        return matriz_caracteristica(λ, p1, p2, p3, p4, p5, τ1, τ2)
    end

    # Muestreo inicial sobre parte real negativa y parte imaginaria moderada
    λs = ComplexF64[]
    for re in range(-1.0, stop=1.0, length=100)
        for im in range(-2.0, stop=2.0, length=100)
            λ₀ = re + im*im
            val = abs(P(λ₀))
            if val < 1e-1
                push!(λs, λ₀)
            end
        end
    end

    if isempty(λs)
        error("No se encontraron raíces aproximadas")
    end

    # Elegir la de mayor parte real
    λ_dominante = argmax(real.(λs))
    return λs[λ_dominante]
end

# (3) Función objetivo para optimización
function objetivo(p)
    p1, p2, p3, p4, p5, τ1, τ2 = p

    try
        λ = raiz_dominante(p1, p2, p3, p4, p5, τ1, τ2)
        return -real(λ)  # minimizar parte real
    catch
        return 1e6  # penaliza si falla
    end
end

# Rango de búsqueda y optimización
res = bboptimize(objetivo;
    SearchRange = [(0.01, 0.1), (0.01, 0.1), (0.01, 1.0), (0.01, 0.3), (0.01, 1.0), (0.0, 20.0), (0.0, 20.0)],
    MaxSteps = 5000,
    Method = :de_rand_1_bin
)

best_params = best_candidate(res)
println("Mejores parámetros encontrados: ", best_params)
