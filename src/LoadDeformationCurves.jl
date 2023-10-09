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
        NumPoints = length(Δ)

        # Preallocate:
        F = zeros(NumPoints, 1)

        # Compute corresponding loads:
        for i in 1:NumPoints
            if i == 1
                F[i] = K * Δ[i]
            else
                F[i] = F[i-1] + α[i] * K * (Δ[i] - Δ[i-1])
            end
        end

        # Add origin:
        Δ = vcat(0, Δ)
        F = vcat(0, F)

        # Clean up:
        Δ = float(vec(Δ))
        F = float(vec(F))

        # Return the updated fields:
        new(K, α, Δ, F)
    end
end