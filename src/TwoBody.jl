# module TwoBody

using Plots, ModelingToolkit, DifferentialEquations, Unitful
include("Utils.jl");
include("OrbitalDynamics.jl");

function LoadTwoBodyVariables()
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

	return t, x₁, y₁, z₁, ẋ₁, ẏ₁, ż₁, x₂, y₂, z₂, ẋ₂, ẏ₂, ż₂, r, D, G, m₁, m₂
end

function LoadTwoBodyVariablesUnitless()
	@variables(begin
		t, 
		x₁(t), 
		y₁(t), 
		z₁(t), 
		ẋ₁(t), 
		ẏ₁(t), 
		ż₁(t), 
		x₂(t), 
		y₂(t), 
		z₂(t), 
		ẋ₂(t), 
		ẏ₂(t), 
		ż₂(t), 
		r(t)
	end)

	D = Differential(t)

	@parameters G m₁ m₂

	return t, x₁, y₁, z₁, ẋ₁, ẏ₁, ż₁, x₂, y₂, z₂, ẋ₂, ẏ₂, ż₂, r, D, G, m₁, m₂
end

function LoadTwoBodyEquations()
	[
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
end

function TwoBodySystem(u₀::Dict{Num, Quantity{T}}, params::Dict{Num, Quantity{T}}, tspan::Vector{T}) where {T<:Number}
    @assert valtype(u₀) == valtype(params) "Numeric values in u₀ and params must be the same type"
	
	two_body_equations = LoadTwoBodyEquations()
	
	diffeq_two_body_system = structural_simplify(ODESystem(
		two_body_equations,
		t,
		name=:two_body_system
	))

    prob = ODEProblem(
        diffeq_two_body_system, 
        remove_units(u₀), 
        tspan, 
        remove_units(params), 
        jac=true
    )

	orbital_bodies = Dict(
		"Body 1" => OrbitalBody("Body 1", params[m₁], BodyState(u₀[x₁], u₀[y₁], u₀[z₁], u₀[ẋ₁], u₀[ẏ₁], u₀[ż₁],)),
		"Body 2" => OrbitalBody("Body 2", params[m₂], BodyState(u₀[x₂], u₀[y₂], u₀[z₂], u₀[ẋ₂], u₀[ẏ₂], u₀[ż₂],)),
	)

    return OrbitalSystem{T}("Two-Body System", orbital_bodies, two_body_equations, diffeq_two_body_system, prob, u₀, params, tspan)
end

function TwoBodySystem(u₀::Dict{Num, T}, params::Dict{Num, T}, tspan::Vector{T}) where {T<:Number}
    @assert valtype(u₀) == valtype(params) "Numeric values in u₀ and params must be the same type"

	two_body_equations = LoadTwoBodyEquations()
	
	diffeq_two_body_system = structural_simplify(ODESystem(
		two_body_equations,
		t,
		name=:two_body_system
	))

    println("Warning: It is recommended to use Unitful.jl units when working with SAT")

    prob = ODEProblem(
        diffeq_two_body_system, 
        u₀, 
        tspan, 
        params, 
        jac=true
    )

	orbital_bodies = Dict(
		"Body 1" => OrbitalBody("Body 1", params[m₁], BodyState(u₀[x₁], u₀[y₁], u₀[z₁], u₀[ẋ₁], u₀[ẏ₁], u₀[ż₁],)),
		"Body 2" => OrbitalBody("Body 2", params[m₂], BodyState(u₀[x₂], u₀[y₂], u₀[z₂], u₀[ẋ₂], u₀[ẏ₂], u₀[ż₂],)),
	)

    return OrbitalSystem{T}("Two-Body System", orbital_bodies, two_body_equations, diffeq_two_body_system, prob, u₀, params, tspan)
end

function OrbitalSystem(body_1::OrbitalBody{T}, body_2::OrbitalBody{T}, tspan::Vector{T}; G_val=6.6743e-11) where {T<:Number}
	orbital_bodies = Dict(
		body_1.name => body_1,
		body_2.name => body_2
	)

	two_body_equations = LoadTwoBodyEquations()
	
	diffeq_two_body_system = structural_simplify(ODESystem(
		two_body_equations,
		t,
		name=:two_body_system
	))
	
	u₀ = Dict(
		x₁ => body_1.state.x,
		y₁ => body_1.state.y,
		z₁ => body_1.state.z,
		ẋ₁ => body_1.state.ẋ,
		ẏ₁ => body_1.state.ẏ,
		ż₁ => body_1.state.ż,
		x₂ => body_2.state.x,
		y₂ => body_2.state.y,
		z₂ => body_2.state.z,
		ẋ₂ => body_2.state.ẋ,
		ẏ₂ => body_2.state.ẏ,
		ż₂ => body_2.state.ż
	)

	params = Dict(	
		G => G_val,
		m₁ => body_1.mass,
		m₂ => body_2.mass
	)

	# check if using units
	if (typeof(body_1.mass) == typeof(body_2.mass)) && typeof(body_1.mass) <: Quantity{T}
		# if units of G_val not defined by user, adding units to default value
		if !(typeof(G_val) <: Quantity{T})
			params[G] = G_val * 1u"N*m^2/kg^2"
		end

		prob = ODEProblem(
			diffeq_two_body_system, 
			remove_units(u₀), 
			tspan, 
			remove_units(params), 
			jac=true
		)
	else
		prob = ODEProblem(
			diffeq_two_body_system, 
			u₀, 
			tspan, 
			params, 
			jac=true
		)
	end

	return OrbitalSystem{T}("Two-Body System", orbital_bodies, two_body_equations, diffeq_two_body_system, prob, u₀, params, tspan)
end