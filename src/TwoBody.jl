# module TwoBody

using Plots, ModelingToolkit, DifferentialEquations, Unitful
include("Utils.jl");
include("OrbitalDynamics.jl");

# export two_body_equations, diffeq_two_body_system, TwoBodySystem, two_body_example_u₀, two_body_example_p, two_body_xample_tspan

const G_val = 6.6743e-11u"N*m^2/kg^2"

@variables(begin
    t, [unit=u"s"],
	x₁(t), [unit=u"m"], 
	y₁(t), [unit=u"m"], 
	z₁(t), [unit=u"m"],
	ẋ₁(t), [unit=u"m/s"], 
	ẏ₁(t), [unit=u"m/s"], 
	ż₁(t), [unit=u"m/s"],
	x₂(t), [unit=u"m"], 
	y₂(t), [unit=u"m"], 
	z₂(t), [unit=u"m"],
	ẋ₂(t), [unit=u"m/s"], 
	ẏ₂(t), [unit=u"m/s"], 
	ż₂(t), [unit=u"m/s"],
    r(t), [unit=u"m"]
end)

D = Differential(t)

@parameters G [unit=u"N*m^2/kg^2"] m₁ [unit=u"kg"] m₂ [unit=u"kg"]

two_body_equations = [
	r ~ sqrt((x₂ - x₁)^2 + (y₂ - y₁)^2 + (z₂ - z₁)^2),
		
	D(x₁) ~ ẋ₁,
	D(y₁) ~ ẏ₁,
	D(z₁) ~ ż₁,

	D(ẋ₁) ~ G*m₂*(x₂ - x₁)/r^3,
	D(ẏ₁) ~ G*m₂*(y₂ - y₁)/r^3,
	D(ż₁) ~ G*m₂*(z₂ - z₁)/r^3,

	D(x₂) ~ ẋ₂,
	D(y₂) ~ ẏ₂,
	D(z₂) ~ ż₂,

	D(ẋ₂) ~ G*m₁*(x₁ - x₂)/r^3,
	D(ẏ₂) ~ G*m₁*(y₁ - y₂)/r^3,
	D(ż₂) ~ G*m₁*(z₁ - z₂)/r^3,
]

diffeq_two_body_system = structural_simplify(ODESystem(
	two_body_equations,
	t,
	name=:two_body_system
));

# example values
two_body_example_u₀ = Dict(
	x₁ => 0.0u"m",
	y₁ => 0.0u"m",
	z₁ => 0.0u"m",
	ẋ₁ => 10.0u"km/s",
	ẏ₁ => 20.0u"km/s",
	ż₁ => 30.0u"km/s",
	x₂ => 3000.0u"km",
	y₂ => 0.0u"m",
	z₂ => 0.0u"m",
	ẋ₂ => 0.0u"m/s",
	ẏ₂ => 40.0u"km/s",
	ż₂ => 0.0u"m/s"
)

two_body_example_p = Dict(
	G => G_val,
	m₁ => 10e26u"kg",
	m₂ => 10e26u"kg"
)

two_body_example_tspan = [0.0, 480.0]

function TwoBodySystem(u₀::Dict{Num, Quantity{T}}, params::Dict{Num, Quantity{T}}, tspan::Vector{T}) where {T<:Number}
    @assert valtype(u₀) == valtype(params) "Numeric values in u₀ and params must be the same type"

    prob = ODEProblem(
        diffeq_two_body_system, 
        remove_units(u₀), 
        tspan, 
        remove_units(params), 
        jac=true
    )

    return OrbitalSystem{T, Quantity{T}}("Two-Body System", two_body_equations, diffeq_two_body_system, prob, u₀, params, tspan)
end

function TwoBodySystem(u₀::Dict{Num, T}, params::Dict{Num, T}, tspan::Vector{T}) where {T<:Number}
    @assert valtype(u₀) == valtype(params) "Numeric values in u₀ and params must be the same type"

    println("Warning: It is recommended to use Unitful.jl units when working with SAT")

    prob = ODEProblem(
        diffeq_two_body_system, 
        u₀, 
        tspan, 
        params, 
        jac=true
    )

    return OrbitalSystem{T, T}("Two-Body System", two_body_equations, diffeq_two_body_system, prob, u₀, params, tspan)
end

# end