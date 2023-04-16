using Plots, ModelingToolkit, DifferentialEquations, Unitful

mutable struct BodyState{T<:Number}
    x::Union{T, Quantity{T}}
    y::Union{T, Quantity{T}}
    z::Union{T, Quantity{T}}
    ẋ::Union{T, Quantity{T}}
    ẏ::Union{T, Quantity{T}}
    ż::Union{T, Quantity{T}}
end

mutable struct OrbitalBody{T<:Number}
    name::String
    mass::Union{T, Quantity{T}}
    state::BodyState{T}
end

struct OrbitalSystem{T<:Number}
    type::String
    orbital_bodies::Dict{String, OrbitalBody}
    eqns::Vector{Equation}
    system::ODESystem
	problem::ODEProblem
    u₀::Dict{Num, Union{T, Quantity{T}}}
    params::Dict{Num, Union{T, Quantity{T}}}
    tspan::Vector{T}
end

function Base.show(io::IO, system::OrbitalSystem)
    out_string = """Orbital System Model: $(system.type)

    Initial Conditions:
    $(system.u₀)
    
    Parameters:
    $(system.params)
    """

    println(out_string)
end

function SolveOrbitalSystem(system::OrbitalSystem, solver::Union{DEAlgorithm,Nothing}; solver_args...)
	return solve(system.problem, solver; solver_args...)
end

function InterpolateOrbitalSolution(solution::ODESolution, vars::Vector{Num}, times::Union{StepRangeLen, Vector})
    sol_interp = solution(times)
    return [sol_interp[var] for var in vars]
end