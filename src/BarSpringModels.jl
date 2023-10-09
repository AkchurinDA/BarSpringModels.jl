module BarSpringModels
# Dependencies:
using LinearAlgebra

# Custom types:
include("LoadDeformationCurves.jl")
export CLDC, PLDC
include("Elements.jl")
export RigidLink
export ConnectionElement

# Custom functions:
include("GetCurrentElementState.jl")
export getCurrentElementK, getCurrentElementΔ, getCurrentElementF
include("Assemble.jl")
export assembleKGlobal, assembleFGlobal
include("Apply.jl")
export applyBC, applyCE
include("Equivalence.jl")
export getEquivalentLoad
include("Analyze.jl")
export analyze
end