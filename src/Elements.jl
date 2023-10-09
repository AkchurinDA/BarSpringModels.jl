abstract type Element end

struct RigidLink <: Element
    NI::Integer # Node (i)
    NJ::Integer # Node (j)
end

struct ConnectionElement <: Element
    NI::Integer # Node (i)
    NJ::Integer # Node (j)
    LDC::Vector{<:LoadDeformationCurve} # Load-deformation curves
end