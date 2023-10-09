# Assemble global stiffness matrix:
function assembleKGlobal(EC::Vector{<:Element}, NN::Integer, NE::Integer, CurrentElementK::Matrix{<:AbstractFloat})
    # Preallocate:
    KGlobal = zeros(6 * NN, 6 * NN)

    # Loop through each element:
    for i in 1:NE
        # Check the element type:
        if EC[i] isa ConnectionElement
            # Extract the nodes indices:
            NI = EC[i].NI
            NJ = EC[i].NJ

            # Loop through each spring:
            for j in 1:6
                # Extract the current spring stiffness:
                K = CurrentElementK[i, j]

                # Place entries into the global stiffness matrix:
                KGlobal[6*NI-(6-j), 6*NI-(6-j)] = KGlobal[6*NI-(6-j), 6*NI-(6-j)] + K
                KGlobal[6*NI-(6-j), 6*NJ-(6-j)] = KGlobal[6*NI-(6-j), 6*NJ-(6-j)] - K
                KGlobal[6*NJ-(6-j), 6*NI-(6-j)] = KGlobal[6*NJ-(6-j), 6*NI-(6-j)] - K
                KGlobal[6*NJ-(6-j), 6*NJ-(6-j)] = KGlobal[6*NJ-(6-j), 6*NJ-(6-j)] + K
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
        if EC[i] isa ConnectionElement
            # Extract the element nodes:
            NI = EC[i].NI
            NJ = EC[i].NJ

            # Place entries into the global internal load vector:
            FGlobal[6*NI-5] = FGlobal[6*NI-5] + CurrentElementF[i, 1]
            FGlobal[6*NI-4] = FGlobal[6*NI-4] + CurrentElementF[i, 2]
            FGlobal[6*NI-3] = FGlobal[6*NI-3] + CurrentElementF[i, 3]
            FGlobal[6*NI-2] = FGlobal[6*NI-2] + CurrentElementF[i, 4]
            FGlobal[6*NI-1] = FGlobal[6*NI-1] + CurrentElementF[i, 5]
            FGlobal[6*NI-0] = FGlobal[6*NI-0] + CurrentElementF[i, 6]
            FGlobal[6*NJ-5] = FGlobal[6*NJ-5] - CurrentElementF[i, 1]
            FGlobal[6*NJ-4] = FGlobal[6*NJ-4] - CurrentElementF[i, 2]
            FGlobal[6*NJ-3] = FGlobal[6*NJ-3] - CurrentElementF[i, 3]
            FGlobal[6*NJ-2] = FGlobal[6*NJ-2] - CurrentElementF[i, 4]
            FGlobal[6*NJ-1] = FGlobal[6*NJ-1] - CurrentElementF[i, 5]
            FGlobal[6*NJ-0] = FGlobal[6*NJ-0] - CurrentElementF[i, 6]
        end
    end

    # Return the result:
    return FGlobal
end