function getEquivalentLoad(Load::Vector{Float64}, NC::Matrix{Float64}, EC::Vector{<:Element}, NE::Integer, U::Vector{Float64})
    # Preallocate:
    UpdatedLoad = deepcopy(Load)

    # Loop through each element:
    for i in 1:NE
        # Check the element type:
        if EC[i] isa RL
            # Extract node indices:
            NI = EC[i].NI # Master node
            NJ = EC[i].NJ # Slave node

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

            # Extract nodal loads:
            F_x_1 = UpdatedLoad[6*NI-5]
            F_y_1 = UpdatedLoad[6*NI-4]
            F_z_1 = UpdatedLoad[6*NI-3]
            M_x_1 = UpdatedLoad[6*NI-2]
            M_y_1 = UpdatedLoad[6*NI-1]
            M_z_1 = UpdatedLoad[6*NI-0]
            F_x_2 = UpdatedLoad[6*NJ-5]
            F_y_2 = UpdatedLoad[6*NJ-4]
            F_z_2 = UpdatedLoad[6*NJ-3]
            M_x_2 = UpdatedLoad[6*NJ-2]
            M_y_2 = UpdatedLoad[6*NJ-1]
            M_z_2 = UpdatedLoad[6*NJ-0]

            # Compute the equivalent loads:
            F_x_eq = -(F_x_1 + F_x_2)
            F_y_eq = -(F_y_1 + F_y_2)
            F_z_eq = -(F_z_1 + F_z_2)
            M_x_eq = -(M_x_1 + M_x_2 - L_z * F_y_2 + L_y * F_z_2)
            M_y_eq = -(M_y_1 + M_y_2 - L_x * F_z_2 + L_z * F_x_2)
            M_z_eq = -(M_z_1 + M_z_2 - L_y * F_x_2 + L_x * F_y_2)

            # Place entries into the load vector:
            UpdatedLoad[6*NI-5] = -F_x_eq
            UpdatedLoad[6*NI-4] = -F_y_eq
            UpdatedLoad[6*NI-3] = -F_z_eq
            UpdatedLoad[6*NI-2] = -M_x_eq
            UpdatedLoad[6*NI-1] = -M_y_eq
            UpdatedLoad[6*NI-0] = -M_z_eq
            UpdatedLoad[6*NJ-5] = 0
            UpdatedLoad[6*NJ-4] = 0
            UpdatedLoad[6*NJ-3] = 0
            UpdatedLoad[6*NJ-2] = 0
            UpdatedLoad[6*NJ-1] = 0
            UpdatedLoad[6*NJ-0] = 0
        end
    end

    # Return the result:
    return UpdatedLoad
end