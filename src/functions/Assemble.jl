# Assemble global stiffness matrix:
function assembleKGlobal(NC::Matrix{Float64}, EC::Vector{<:Element}, NN::Integer, NE::Integer, U::Vector{Float64}, CurrentElementK::Matrix{Float64})
    # Preallocate:
    KGlobal = zeros(6 * NN, 6 * NN)

    # Loop through each element:
    for i in 1:NE
        # Check the element type:
        if !isa(EC[i], RL) && !isa(EC[i], CE1) && !isa(EC[i], CE2)
            error("Incorrect element type!")
        elseif isa(EC[i], CE1) # Connection element with decoupled translational DOFs
            # Extract the nodes indices:
            NI = EC[i].NI
            NJ = EC[i].NJ

            # Loop through each spring:
            for j in 1:6
                # Extract the current spring stiffness:
                K = CurrentElementK[i, j] * [1 -1; -1 1]

                # Place entries into the global stiffness matrix:
                KGlobal[6*NI-(6-j), 6*NI-(6-j)] = KGlobal[6*NI-(6-j), 6*NI-(6-j)] + K[1, 1]
                KGlobal[6*NI-(6-j), 6*NJ-(6-j)] = KGlobal[6*NI-(6-j), 6*NJ-(6-j)] + K[1, 2]
                KGlobal[6*NJ-(6-j), 6*NI-(6-j)] = KGlobal[6*NJ-(6-j), 6*NI-(6-j)] + K[2, 1]
                KGlobal[6*NJ-(6-j), 6*NJ-(6-j)] = KGlobal[6*NJ-(6-j), 6*NJ-(6-j)] + K[2, 2]
            end
        elseif isa(EC[i], CE2) # Connection element with coupled translational DOFs
            # Extract the nodes indices:
            NI = EC[i].NI
            NJ = EC[i].NJ

            # Loop through each spring:
            for j in 1:4
                if j == 1 # Translational spring
                    # Extract nodal coordinates:
                    x_1 = NC[NI, 1]
                    y_1 = NC[NI, 2]
                    z_1 = NC[NI, 3]
                    x_2 = NC[NJ, 1]
                    y_2 = NC[NJ, 2]
                    z_2 = NC[NJ, 3]

                    # Extract nodal translations:
                    u_x_1 = U[6*NI-5]
                    u_y_1 = U[6*NI-4]
                    u_z_1 = U[6*NI-3]
                    u_x_2 = U[6*NJ-5]
                    u_y_2 = U[6*NJ-4]
                    u_z_2 = U[6*NJ-3]

                    # Compute deformed length:
                    L = sqrt((x_2 + u_x_2 - x_1 - u_x_1)^2 + (y_2 + u_y_2 - y_1 - u_y_1)^2 + (z_2 + u_z_2 - z_1 - u_z_1)^2)

                    # Compute local-to-global transformation matrix:
                    C_x = (x_2 + u_x_2 - x_1 - u_x_1) / L
                    C_y = (y_2 + u_y_2 - y_1 - u_y_1) / L
                    C_z = (z_2 + u_z_2 - z_1 - u_z_1) / L
                    T = [C_x C_y C_z 0 0 0; 0 0 0 C_x C_y C_z]

                    # Compute element stiffness matrix in global coordinate system:
                    K = transpose(T) * (CurrentElementK[i, j] * [1 -1; -1 1]) * T

                    # Place entries into the global stiffness matrix:
                    KGlobal[(6*NI-5):(6*NI-3), (6*NI-5):(6*NI-3)] = KGlobal[(6*NI-5):(6*NI-3), (6*NI-5):(6*NI-3)] + K[1:3, 1:3]
                    KGlobal[(6*NI-5):(6*NI-3), (6*NJ-5):(6*NJ-3)] = KGlobal[(6*NI-5):(6*NI-3), (6*NJ-5):(6*NJ-3)] + K[1:3, 4:6]
                    KGlobal[(6*NJ-5):(6*NJ-3), (6*NI-5):(6*NI-3)] = KGlobal[(6*NJ-5):(6*NJ-3), (6*NI-5):(6*NI-3)] + K[4:6, 1:3]
                    KGlobal[(6*NJ-5):(6*NJ-3), (6*NJ-5):(6*NJ-3)] = KGlobal[(6*NJ-5):(6*NJ-3), (6*NJ-5):(6*NJ-3)] + K[4:6, 4:6]
                else # Rotational springs
                    # Adjust the index:
                    J = j + 2

                    # Extract the current spring stiffness:
                    K = CurrentElementK[i, j] * [1 -1; -1 1]

                    # Place entries into the global stiffness matrix:
                    KGlobal[6*NI-(6-J), 6*NI-(6-J)] = KGlobal[6*NI-(6-J), 6*NI-(6-J)] + K[1, 1]
                    KGlobal[6*NI-(6-J), 6*NJ-(6-J)] = KGlobal[6*NI-(6-J), 6*NJ-(6-J)] + K[1, 2]
                    KGlobal[6*NJ-(6-J), 6*NI-(6-J)] = KGlobal[6*NJ-(6-J), 6*NI-(6-J)] + K[2, 1]
                    KGlobal[6*NJ-(6-J), 6*NJ-(6-J)] = KGlobal[6*NJ-(6-J), 6*NJ-(6-J)] + K[2, 2]
                end
            end
        end
    end

    # Return the result:
    return KGlobal
end

# Assemble global internal load vector:
function assembleFGlobal(EC, NN, NE, CurrentElementF)
    # Preallocate:
    FGlobal = zeros(6 * NN)

    # Loop thorugh each element:
    for i in 1:NE
        # Check the element type:
        if !isa(EC[i], RL) && !isa(EC[i], CE1) && !isa(EC[i], CE2)
            error("Incorrect element type!")
        elseif isa(EC[i], CE1)
            # Extract the element nodes:
            NI = EC[i].NI
            NJ = EC[i].NJ

            # Loop through each spring:
            for j in 1:6
                # Convert local internal loads into global internal loads:
                F = CurrentElementF[i, j] * [-1, 1] # 2 x 1

                # Place entries into the global internal load vector:
                FGlobal[6*NI-(6-j)] = FGlobal[6*NI-(6-j)] + F[1]
                FGlobal[6*NJ-(6-j)] = FGlobal[6*NJ-(6-j)] + F[2]
            end
        elseif isa(EC[i], CE2)
            # Extract the element nodes:
            NI = EC[i].NI
            NJ = EC[i].NJ

            # Loop through each spring:
            for j in 1:4
                if j == 1 # Translational spring
                    # Extract nodal coordinates:
                    x_1 = NC[NI, 1]
                    y_1 = NC[NI, 2]
                    z_1 = NC[NI, 3]
                    x_2 = NC[NJ, 1]
                    y_2 = NC[NJ, 2]
                    z_2 = NC[NJ, 3]

                    # Extract nodal translations:
                    u_x_1 = U[6*NI-5]
                    u_y_1 = U[6*NI-4]
                    u_z_1 = U[6*NI-3]
                    u_x_2 = U[6*NJ-5]
                    u_y_2 = U[6*NJ-4]
                    u_z_2 = U[6*NJ-3]

                    # Compute deformed length:
                    L = sqrt((x_2 + u_x_2 - x_1 - u_x_1)^2 + (y_2 + u_y_2 - y_1 - u_y_1)^2 + (z_2 + u_z_2 - z_1 - u_z_1)^2)

                    # Compute local-to-global transformation matrix:
                    C_x = (x_2 + u_x_2 - x_1 - u_x_1) / L
                    C_y = (y_2 + u_y_2 - y_1 - u_y_1) / L
                    C_z = (z_2 + u_z_2 - z_1 - u_z_1) / L
                    T = [C_x C_y C_z 0 0 0; 0 0 0 C_x C_y C_z]

                    # Convert local internal loads into global internal loads:
                    F = transpose(T) * (CurrentElementF[i, j] * [-1, 1]) # 6 x 1

                    # Place entries into the global internal load vector:
                    FGlobal[6*NI-5] = FGlobal[6*NI-5] + F[1]
                    FGlobal[6*NI-4] = FGlobal[6*NI-4] + F[2]
                    FGlobal[6*NI-3] = FGlobal[6*NI-3] + F[3]
                    FGlobal[6*NJ-5] = FGlobal[6*NJ-5] + F[4]
                    FGlobal[6*NJ-4] = FGlobal[6*NJ-4] + F[5]
                    FGlobal[6*NJ-3] = FGlobal[6*NJ-3] + F[6]
                else # Rotational springs
                    # Adjust the index:
                    J = j + 2

                    # Convert local internal loads into global internal loads:
                    F = CurrentElementF[i, J] * [-1, 1] # 2 x 1

                    # Place entries into the global internal load vector:
                    FGlobal[6*NI-(6-J)] = FGlobal[6*NI-(6-J)] + F[1]
                    FGlobal[6*NJ-(6-J)] = FGlobal[6*NJ-(6-J)] + F[2]
                end
            end
        end
    end

    # Return the result:
    return FGlobal
end