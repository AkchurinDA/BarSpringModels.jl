module BarSpringModels
# Dependencies:
using LinearAlgebra

# Types:
include("types/LoadDeformationTypes.jl")
include("types/Elements.jl")
export CLDC
export PLDC
export RB
export CE

# Functions:
include("functions/Analyze.jl")
include("functions/Assemble.jl")
include("functions/Apply.jl")
include("functions/GetCurrentElementState.jl")
include("functions/GetEquivalentLoads.jl")
export analyze
export assembleKGlobal
export assembleFGlobal
export applyBC
export applyCE
export getCurrentElementΔ
export getCurrentElementK
export getCurrentElementF
export getEquivalentLoad
end