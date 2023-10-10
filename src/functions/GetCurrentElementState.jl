# Get current deformations of each spring of a connection element:
function getCurrentElementΔ(EC::Vector{<:Element}, NE::Integer, U::Vector{Float64})
    # Preallocate:
    CurrentElementΔ = zeros(NE, 6)

    # Loop through each element:
    for i in 1:NE
        # Check the type of element:
        if isa(EC[i], CE)
            # Extract node indices:
            NI = EC[i].NI
            NJ = EC[i].NJ

            # Extract nodal translations and rotations:
            u_x_1 = U[6*NI-5]
            u_y_1 = U[6*NI-4]
            u_z_1 = U[6*NI-3]
            θ_x_1 = U[6*NI-2]
            θ_y_1 = U[6*NI-1]
            θ_z_1 = U[6*NI-0]
            u_x_2 = U[6*NJ-5]
            u_y_2 = U[6*NJ-4]
            u_z_2 = U[6*NJ-3]
            θ_x_2 = U[6*NJ-2]
            θ_y_2 = U[6*NJ-1]
            θ_z_2 = U[6*NJ-0]

            # Compute current element deformation:
            CurrentElementΔ[i, 1] = u_x_2 - u_x_1
            CurrentElementΔ[i, 2] = u_y_2 - u_y_1
            CurrentElementΔ[i, 3] = u_z_2 - u_z_1
            CurrentElementΔ[i, 4] = θ_x_2 - θ_x_1
            CurrentElementΔ[i, 5] = θ_y_2 - θ_y_1
            CurrentElementΔ[i, 6] = θ_z_2 - θ_z_1
        end
    end

    # Return the result:
    return CurrentElementΔ
end

# Get current stiffnesses of each spring of a connection element:
function getCurrentElementK(EC::Vector{<:Element}, NE::Integer, CurrentElementΔ::Matrix{Float64})
    # Preallocate:
    CurrentElementK = zeros(NE, 6)

    # Loop through each element:
    for i in 1:NE
        # Check the type of element:
        if isa(EC[i], CE)
            # Loop through each spring:
            for j in 1:6
                # Check the type of the load-deformation curve:
                if !isa(EC[i].LDC[j], CLDC) && !isa(EC[i].LDC[j], PLDC) # Error-catching
                    error("Incorrect load-deformation curve type!")
                elseif isa(EC[i].LDC[j], CLDC) # Load-deformation curve with constant stiffness
                    # Compute the current stiffness:
                    CurrentElementK[i, j] = EC[i].LDC[j].K
                elseif isa(EC[i].LDC[j], PLDC) # Piecewise linear load-deformation curve
                    # Compute the current location index:
                    Index = findlast(x -> x ≤ abs(CurrentElementΔ[i, j]), EC[i].LDC[j].Δ)

                    # Compute the current stiffness:
                    if Index == 1
                        CurrentElementK[i, j] = EC[i].LDC[j].K
                    else
                        CurrentElementK[i, j] = EC[i].LDC[j].α[Index-1] * EC[i].LDC[j].K
                    end
                end
            end
        end
    end

    # Return the result:
    return CurrentElementK
end

# Get current internal loads of each spring of a connection element:
function getCurrentElementF(EC::Vector{<:Element}, NE::Integer, CurrentElementΔ::Matrix{Float64}, CurrentElementK::Matrix{Float64})
    # Preallocate:
    CurrentElementF = zeros(NE, 6)

    # Loop through each element:
    for i in 1:NE
        # Check the type of element:
        if EC[i] isa CE
            # Loop through each spring:
            for j in 1:6
                # Check the type of the load-deformation curve:
                if !isa(EC[i].LDC[j], CLDC) && !isa(EC[i].LDC[j], PLDC) # Error-catching
                    error("Incorrect load-deformation curve type!")
                elseif isa(EC[i].LDC[j], CLDC) # Load-deformation curve with constant stiffness
                    # Compute the current internal load:
                    CurrentElementF[i, j] = CurrentElementΔ[i, j] * CurrentElementK[i, j]
                elseif isa(EC[i].LDC[j], PLDC) # Piecewise linear load-deformation curve
                    # Compute the current location index:
                    Index = findlast(x -> x ≤ abs(CurrentElementΔ[i, j]), EC[i].LDC[j].Δ)

                    # Compute the current internal load:
                    if Index == 1
                        CurrentElementF[i, j] = CurrentElementK[i, j] * CurrentElementΔ[i, j]
                    else
                        CurrentElementF[i, j] = EC[i].LDC[j].F[Index] + CurrentElementK[i, j] * (CurrentElementΔ[i, j] - EC[i].LDC[j].Δ[Index])
                    end
                end
            end
        end
    end

    # Return the result:
    return CurrentElementF
end