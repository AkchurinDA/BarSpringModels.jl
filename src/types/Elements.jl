abstract type Element end

# Rigid link:
struct RL <: Element
    NI::Integer
    NJ::Integer
end

# Connection element with decoupled translational DOFs:
struct CE1 <: Element
    NI::Integer
    NJ::Integer
    LDC::Vector{<:LoadDeformationCurve} # [T_x, T_y, T_z, R_x, R_y, R_z]
end

# Connection element with coupled translational DOFs:
struct CE2 <: Element
    NI::Integer
    NJ::Integer
    LDC::Vector{<:LoadDeformationCurve} # [T_x, R_x, R_y, R_z]
end