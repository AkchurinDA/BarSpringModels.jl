using BarSpringModels
using CairoMakie, MathTeXEngine

# Define the number of nodes and elements in the model:
NN = 4
NE = 3

# Define nodal coordinates:
NC = [
    0 0 0
    1.44 72 0
    1.44 72 0
    0 144 0
]

# Define element connectivity:
EC = [
    RigidLink(1, 2)
    ConnectionElement(2, 3, [CLDC(10^8), CLDC(10^8), CLDC(10^8), CLDC(10^8), CLDC(10^8), PLDC(314 / 0.055, 0, 0.055)])
    RigidLink(3, 4)
]

# Define boundary conditions:
BC = zeros(NN, 6)
BC[1, :] = [1, 1, 1, 1, 1, 0]
BC[2, :] = [0, 0, 1, 1, 1, 0]
BC[3, :] = [0, 0, 1, 1, 1, 0]
BC[4, :] = [1, 0, 1, 1, 1, 0]

# Define applied reference load vector:
P̃ = zeros(NN, 6)
P̃[4, :] = [0, -1, 0, 0, 0, 0]

# Solve the model:
Solution = analyze(NC, EC, BC, P̃, 10, 100, 20, 10^(-3), 10^(-3))

# Plot the results:
P = vcat(0, abs.(Solution[1])) ./ 155.974
U = (1.44 .+ vcat(0, abs.(Solution[3][7, :]))) ./ 144

begin
    F = Figure(resolution=72 .* (8, 6), fonts=(; regular=texfont()))
    A = Axis(F[1, 1],
        xlabel=L"$\delta_0/l$", ylabel=L"$P/P_{\text{cr}}$",
        xminorticksvisible=true, yminorticksvisible=true,
        xminorgridvisible=true, yminorgridvisible=true)
    scatterlines!(U, P,
        markersize=6)
    xlims!(0, nothing)
    ylims!(0, nothing)
    F
end
