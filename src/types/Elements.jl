abstract type Element end

# Rigid bar:
struct RB <: Element
    NI::Integer
    NJ::Integer
end

# Connection element:
struct CE <: Element
    NI::Integer
    NJ::Integer
    LDC::Vector{<:LoadDeformationCurve}
end