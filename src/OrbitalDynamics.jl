# module OrbitalDynamics

using Plots, ModelingToolkit, DifferentialEquations, Unitful

# export OrbitalSystem, Base.show, SolveOrbitalSystem

mutable struct OrbitalSystem{T<:Number, Quant<:Number}
    type::String
    eqns::Vector{Equation}
    system::ODESystem
	problem::ODEProblem
    u₀::Dict{Num, Quant}
    params::Dict{Num, Quant}
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

# end