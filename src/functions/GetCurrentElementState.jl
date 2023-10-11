# Get current deformations of each spring of a connection element:
function getCurrentElementΔ(NC::Matrix{Float64}, EC::Vector{<:Element}, NE::Integer, U::Vector{Float64})
    # Preallocate:
    CurrentElementΔ = zeros(NE, 6)

    # Loop through each element:
    for i in 1:NE
        # Check the type of element:
        if !isa(EC[i], RL) && !isa(EC[i], CE1) && !isa(EC[i], CE2)
            error("Incorrect element type!")
        elseif isa(EC[i], CE1)
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
        elseif isa(EC[i], CE2)
            # Extract node indices:
            NI = EC[i].NI
            NJ = EC[i].NJ

            # Extract nodal coordinates:
            x_1 = NC[NI, 1]
            y_1 = NC[NI, 2]
            z_1 = NC[NI, 3]
            x_2 = NC[NJ, 1]
            y_2 = NC[NJ, 2]
            z_2 = NC[NJ, 3]

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

            # Compute undeformed and deformed lengths:
            L_1 = sqrt((x_2 - x_1)^2 + (y_2 - y_1)^2 + (z_2 - z_1)^2) # Undeformed
            L_2 = sqrt((x_2 + u_x_2 - x_1 - u_x_1)^2 + (y_2 + u_y_2 - y_1 - u_y_1)^2 + (z_2 + u_z_2 - z_1 - u_z_1)^2) # Deformed

            # Compute current element deformation:
            CurrentElementΔ[i, 1] = L_2 - L_1
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
        if !isa(EC[i], RL) && !isa(EC[i], CE1) && !isa(EC[i], CE2)
            error("Incorrect element type!")
        elseif isa(EC[i], CE1)
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
        elseif isa(EC[i], CE2)
            # Loop through each spring:
            for j in 1:4
                if j == 1 # Translational spring
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
                else # Rotational springs
                    # Adjust the index:
                    J = j + 2

                    # Check the type of the load-deformation curve:
                    if !isa(EC[i].LDC[j], CLDC) && !isa(EC[i].LDC[j], PLDC) # Error-catching
                        error("Incorrect load-deformation curve type!")
                    elseif isa(EC[i].LDC[j], CLDC) # Load-deformation curve with constant stiffness
                        # Compute the current stiffness:
                        CurrentElementK[i, J] = EC[i].LDC[j].K
                    elseif isa(EC[i].LDC[j], PLDC) # Piecewise linear load-deformation curve
                        # Compute the current location index:
                        Index = findlast(x -> x ≤ abs(CurrentElementΔ[i, J]), EC[i].LDC[j].Δ)

                        # Compute the current stiffness:
                        if Index == 1
                            CurrentElementK[i, J] = EC[i].LDC[j].K
                        else
                            CurrentElementK[i, J] = EC[i].LDC[j].α[Index-1] * EC[i].LDC[j].K
                        end
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
        if !isa(EC[i], RL) && !isa(EC[i], CE1) && !isa(EC[i], CE2)
            error("Incorrect element type!")
        elseif isa(EC[i], CE1)
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
        elseif isa(EC[i], CE2)
            # Loop through each spring:
            for j in 1:4
                if j == 1 # Translational spring
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
                else # Rotational springs
                    # Adjust the index:
                    J = j + 2

                    # Check the type of the load-deformation curve:
                    if !isa(EC[i].LDC[J], CLDC) && !isa(EC[i].LDC[J], PLDC) # Error-catching
                        error("Incorrect load-deformation curve type!")
                    elseif isa(EC[i].LDC[J], CLDC) # Load-deformation curve with constant stiffness
                        # Compute the current internal load:
                        CurrentElementF[i, J] = CurrentElementΔ[i, J] * CurrentElementK[i, J]
                    elseif isa(EC[i].LDC[J], PLDC) # Piecewise linear load-deformation curve
                        # Compute the current location index:
                        Index = findlast(x -> x ≤ abs(CurrentElementΔ[i, J]), EC[i].LDC[J].Δ)

                        # Compute the current internal load:
                        if Index == 1
                            CurrentElementF[i, J] = CurrentElementK[i, J] * CurrentElementΔ[i, J]
                        else
                            CurrentElementF[i, J] = EC[i].LDC[J].F[Index] + CurrentElementK[i, J] * (CurrentElementΔ[i, J] - EC[i].LDC[J].Δ[Index])
                        end
                    end
                end
            end
        end
    end

    # Return the result:
    return CurrentElementF
end