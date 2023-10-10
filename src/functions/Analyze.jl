function analyze(NC::Matrix{Float64}, EC::Vector{<:Element}, BC::Matrix{Float64}, P̃::Matrix{Float64}, iₘ::Integer, jₘ::Integer, ControlDOF::Integer, ControlFactor::Real, ϵ::Real)
    # Compute the number of nodes and elements:
    NN = size(NC, 1) # Number of nodes
    NE = size(EC, 1) # Number of elements

    # Flatten out matrices into vectors:
    BC = Vector(vec(transpose(BC)))
    P̃ = Vector(vec(transpose(P̃)))

    # Initialize:
    λ = 0
    P = zeros(6 * NN) # Total applied load vector
    U = zeros(6 * NN) # Total nodal displacement vector

    # Compute the equivalent applied reference load vector:
    P̃Equivalent = getEquivalentLoad(P̃, NC, EC, NE, U)

    # Initialize termination flags:
    F1 = false # Singularity of the global stiffness matrix
    F2 = false # Maximum number of iterations (j) is reached
    F3 = false # Maximum number of increments (i) is reached

    # Preallocate:
    λStore = zeros(iₘ,)
    PStore = zeros(6 * NN, iₘ)
    UStore = zeros(6 * NN, iₘ)

    # Loop through increments (i):
    for i in 1:iₘ
        # Initialize:
        R = zeros(6 * NN) # Residual load vector

        # Loop through iterations (j):
        for j in 1:jₘ
            # Compute current element deformations:
            CurrentElementΔ = getCurrentElementΔ(EC, NE, U)

            # Compute current element stiffnesses:
            CurrentElementK = getCurrentElementK(EC, NE, CurrentElementΔ)

            # Assemble global stiffness matrix:
            KGlobal = assembleKGlobal(EC, NN, NE, CurrentElementK)

            # Apply boundary conditions and constraint equations to the global stiffness matrix:
            KGlobal = applyBC(KGlobal, BC, 10^8)
            KGlobal = applyCE(KGlobal, NC, EC, NE, U, 10^8)

            # Check if the global stiffness matrix is singular:
            if abs(det(KGlobal)) ≤ eps()
                println("Global stiffness matrix is singular! Process is terminated!")
                F1 = true
                break
            end

            # Compute nodal displacement increments due to applied reference and residual loads:
            δUP = KGlobal \ P̃Equivalent
            δUR = KGlobal \ R

            # Compute the load factor increment:
            if j == 1
                δλ = ControlFactor / δUP[ControlDOF]
            else
                δλ = -δUR[ControlDOF] / δUP[ControlDOF]
            end

            # Update:
            λ = λ + δλ
            P = P + δλ * P̃Equivalent
            U = U + δλ * δUP + δUR

            # Compute current element deformations:
            CurrentElementΔ = getCurrentElementΔ(EC, NE, U)

            # Compute current element stiffnesses:
            CurrentElementK = getCurrentElementK(EC, NE, CurrentElementΔ)

            # Compute current element internal loads:
            CurrentElementF = getCurrentElementF(EC, NE, CurrentElementΔ, CurrentElementK)

            # Assemble global internal load vector:
            FGlobal = assembleFGlobal(EC, NN, NE, CurrentElementF)

            # Compute the equivalent global internal load vector:
            FGlobalEquivalent = getEquivalentLoad(FGlobal, NC, EC, NE, U)

            # Compute the residual load vector:
            R = P - FGlobalEquivalent
            # display(hcat(P̃, P̃Equivalent, P, FGlobal, FGlobalEquivalent, R, BC))
            # display(norm(R[BC.==0.0]))

            # Check for convergance:
            if norm(R[BC.==0.0]) ≤ ϵ
                λStore[i] = λ
                PStore[:, i] = P
                UStore[:, i] = U
                break
            end

            # Check if the maximum iteration is reached:
            if j == jₘ
                println("Maximum number of iterations (j) is reached! Process is terminated!")
                F2 = true
                break
            end
        end

        # Check if the maximum increment is reached:
        if i == iₘ
            println("Maximum number of increments (i) is reached! Process is terminated!")
            F3 = true
        end

        # Check flags:
        if F1 || F2 || F3
            λStore = λStore[1:i-1]
            PStore = PStore[:, 1:i-1]
            UStore = UStore[:, 1:i-1]
            break
        end
    end

    # Return the result:
    return λStore, PStore, UStore
end