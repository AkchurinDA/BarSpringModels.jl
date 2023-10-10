abstract type LoadDeformationCurve end

# Load-deformation curve with constant stiffness:
struct CLDC <: LoadDeformationCurve
    K::Real # Stiffness
end

# Piecewise linear load-deformation curve: 
struct PLDC <: LoadDeformationCurve
    K::Real # Initial stiffness
    α::Union{Real,Vector{<:Real}} # Post-yield slopes
    Δ::Vector{Float64} # Deformations
    F::Vector{Float64} # Loads

    function PLDC(K::Real, α::Union{Real,Vector{<:Real}}, Δ::Union{Real,Vector{<:Real}})
        # Compute the number of inflection points:
        NIP = length(Δ)

        # Preallocate:
        F = zeros(NIP)

        # Compute corresponding loads:
        for i in 1:NIP
            if i == 1
                F[i] = K * Δ[i]
            else
                F[i] = F[i-1] + α[i] * K * (Δ[i] - Δ[i-1])
            end
        end

        # Clean up:
        Δ = float(vec(vcat(0, Δ)))
        F = float(vec(vcat(0, F)))

        # Return the result:
        new(K, α, Δ, F)
    end
end