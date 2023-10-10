# Apply boundary conditions to the global stiffness matrix:
function applyBC(KGlobal::Matrix{Float64}, BC::Vector{Float64}, Penalty::Real)
    # Preallocate:
    KGlobalUpdated = deepcopy(KGlobal)

    # Compute the number of DOFs:
    NDOF = length(BC)

    # Apply boundary conditions using penalty method:
    for i in 1:NDOF
        # Check if the current DOF is constrained:
        if BC[i] == 1
            KGlobalUpdated[i, i] = KGlobalUpdated[i, i] + Penalty
        end
    end

    # Return the result:
    return KGlobalUpdated
end

# Apply constraint equations to the global stiffness matrix:
function applyCE(KGlobal::Matrix{Float64}, NC::Matrix{Float64}, EC::Vector{<:Element}, NE::Integer, U::Vector{Float64}, Penalty::Real)
    # Preallocate:
    KGlobalUpdated = deepcopy(KGlobal)

    # Loop through each element:
    for i in 1:NE
        # Check the element type:
        if EC[i] isa RB
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

            # Extract nodal translations:
            u_x_1 = U[6*NI-5]
            u_y_1 = U[6*NI-4]
            u_z_1 = U[6*NI-3]
            u_x_2 = U[6*NJ-5]
            u_y_2 = U[6*NJ-4]
            u_z_2 = U[6*NJ-3]

            # Compute lengths:
            L_x = x_2 + u_x_2 - x_1 - u_x_1
            L_y = y_2 + u_y_2 - y_1 - u_y_1
            L_z = z_2 + u_z_2 - z_1 - u_z_1

            # Equation #1: u_x_1 + L_z * θ_y_1 - L_y * θ_z_1 - u_x_2 = 0
            V = [1, 0, 0, 0, L_z, -L_y, -1, 0, 0, 0, 0, 0]
            A = V * transpose(V)
            KGlobalUpdated[(6*NI-5):(6*NI-0), (6*NI-5):(6*NI-0)] = KGlobalUpdated[(6*NI-5):(6*NI-0), (6*NI-5):(6*NI-0)] + Penalty * A[1:6, 1:6]
            KGlobalUpdated[(6*NI-5):(6*NI-0), (6*NJ-5):(6*NJ-0)] = KGlobalUpdated[(6*NI-5):(6*NI-0), (6*NJ-5):(6*NJ-0)] + Penalty * A[1:6, 7:12]
            KGlobalUpdated[(6*NJ-5):(6*NJ-0), (6*NI-5):(6*NI-0)] = KGlobalUpdated[(6*NJ-5):(6*NJ-0), (6*NI-5):(6*NI-0)] + Penalty * A[7:12, 1:6]
            KGlobalUpdated[(6*NJ-5):(6*NJ-0), (6*NJ-5):(6*NJ-0)] = KGlobalUpdated[(6*NJ-5):(6*NJ-0), (6*NJ-5):(6*NJ-0)] + Penalty * A[7:12, 7:12]

            # Equation #2: u_y_1 - L_z * θ_x_1 + L_x * θ_z_1 - u_y_2 = 0
            V = [0, 1, 0, -L_z, 0, L_x, 0, -1, 0, 0, 0, 0]
            A = V * transpose(V)
            KGlobalUpdated[(6*NI-5):(6*NI-0), (6*NI-5):(6*NI-0)] = KGlobalUpdated[(6*NI-5):(6*NI-0), (6*NI-5):(6*NI-0)] + Penalty * A[1:6, 1:6]
            KGlobalUpdated[(6*NI-5):(6*NI-0), (6*NJ-5):(6*NJ-0)] = KGlobalUpdated[(6*NI-5):(6*NI-0), (6*NJ-5):(6*NJ-0)] + Penalty * A[1:6, 7:12]
            KGlobalUpdated[(6*NJ-5):(6*NJ-0), (6*NI-5):(6*NI-0)] = KGlobalUpdated[(6*NJ-5):(6*NJ-0), (6*NI-5):(6*NI-0)] + Penalty * A[7:12, 1:6]
            KGlobalUpdated[(6*NJ-5):(6*NJ-0), (6*NJ-5):(6*NJ-0)] = KGlobalUpdated[(6*NJ-5):(6*NJ-0), (6*NJ-5):(6*NJ-0)] + Penalty * A[7:12, 7:12]

            # Equation #3: u_z_1 + L_y * θ_x_1 - L_x * θ_y_1 - u_z_2 = 0
            V = [0, 0, 1, L_y, -L_x, 0, 0, 0, -1, 0, 0, 0]
            A = V * transpose(V)
            KGlobalUpdated[(6*NI-5):(6*NI-0), (6*NI-5):(6*NI-0)] = KGlobalUpdated[(6*NI-5):(6*NI-0), (6*NI-5):(6*NI-0)] + Penalty * A[1:6, 1:6]
            KGlobalUpdated[(6*NI-5):(6*NI-0), (6*NJ-5):(6*NJ-0)] = KGlobalUpdated[(6*NI-5):(6*NI-0), (6*NJ-5):(6*NJ-0)] + Penalty * A[1:6, 7:12]
            KGlobalUpdated[(6*NJ-5):(6*NJ-0), (6*NI-5):(6*NI-0)] = KGlobalUpdated[(6*NJ-5):(6*NJ-0), (6*NI-5):(6*NI-0)] + Penalty * A[7:12, 1:6]
            KGlobalUpdated[(6*NJ-5):(6*NJ-0), (6*NJ-5):(6*NJ-0)] = KGlobalUpdated[(6*NJ-5):(6*NJ-0), (6*NJ-5):(6*NJ-0)] + Penalty * A[7:12, 7:12]

            # Equation #4: θ_x_1 - θ_x_2 = 0
            V = [0, 0, 0, 1, 0, 0, 0, 0, 0, -1, 0, 0]
            A = V * transpose(V)
            KGlobalUpdated[(6*NI-5):(6*NI-0), (6*NI-5):(6*NI-0)] = KGlobalUpdated[(6*NI-5):(6*NI-0), (6*NI-5):(6*NI-0)] + Penalty * A[1:6, 1:6]
            KGlobalUpdated[(6*NI-5):(6*NI-0), (6*NJ-5):(6*NJ-0)] = KGlobalUpdated[(6*NI-5):(6*NI-0), (6*NJ-5):(6*NJ-0)] + Penalty * A[1:6, 7:12]
            KGlobalUpdated[(6*NJ-5):(6*NJ-0), (6*NI-5):(6*NI-0)] = KGlobalUpdated[(6*NJ-5):(6*NJ-0), (6*NI-5):(6*NI-0)] + Penalty * A[7:12, 1:6]
            KGlobalUpdated[(6*NJ-5):(6*NJ-0), (6*NJ-5):(6*NJ-0)] = KGlobalUpdated[(6*NJ-5):(6*NJ-0), (6*NJ-5):(6*NJ-0)] + Penalty * A[7:12, 7:12]

            # Equation #5: θ_y_1 - θ_y_2 = 0
            V = [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, -1, 0]
            A = V * transpose(V)
            KGlobalUpdated[(6*NI-5):(6*NI-0), (6*NI-5):(6*NI-0)] = KGlobalUpdated[(6*NI-5):(6*NI-0), (6*NI-5):(6*NI-0)] + Penalty * A[1:6, 1:6]
            KGlobalUpdated[(6*NI-5):(6*NI-0), (6*NJ-5):(6*NJ-0)] = KGlobalUpdated[(6*NI-5):(6*NI-0), (6*NJ-5):(6*NJ-0)] + Penalty * A[1:6, 7:12]
            KGlobalUpdated[(6*NJ-5):(6*NJ-0), (6*NI-5):(6*NI-0)] = KGlobalUpdated[(6*NJ-5):(6*NJ-0), (6*NI-5):(6*NI-0)] + Penalty * A[7:12, 1:6]
            KGlobalUpdated[(6*NJ-5):(6*NJ-0), (6*NJ-5):(6*NJ-0)] = KGlobalUpdated[(6*NJ-5):(6*NJ-0), (6*NJ-5):(6*NJ-0)] + Penalty * A[7:12, 7:12]

            # Equation #6: θ_z_1 - θ_z_2 = 0
            V = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, -1]
            A = V * transpose(V)
            KGlobalUpdated[(6*NI-5):(6*NI-0), (6*NI-5):(6*NI-0)] = KGlobalUpdated[(6*NI-5):(6*NI-0), (6*NI-5):(6*NI-0)] + Penalty * A[1:6, 1:6]
            KGlobalUpdated[(6*NI-5):(6*NI-0), (6*NJ-5):(6*NJ-0)] = KGlobalUpdated[(6*NI-5):(6*NI-0), (6*NJ-5):(6*NJ-0)] + Penalty * A[1:6, 7:12]
            KGlobalUpdated[(6*NJ-5):(6*NJ-0), (6*NI-5):(6*NI-0)] = KGlobalUpdated[(6*NJ-5):(6*NJ-0), (6*NI-5):(6*NI-0)] + Penalty * A[7:12, 1:6]
            KGlobalUpdated[(6*NJ-5):(6*NJ-0), (6*NJ-5):(6*NJ-0)] = KGlobalUpdated[(6*NJ-5):(6*NJ-0), (6*NJ-5):(6*NJ-0)] + Penalty * A[7:12, 7:12]
        end
    end

    # Return the result:
    return KGlobalUpdated
end