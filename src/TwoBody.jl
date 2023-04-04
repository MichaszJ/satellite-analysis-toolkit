using Plots, ModelingToolkit, DifferentialEquations, Unitful

# helper functions for dealing with units in initial conditions and parameter dicts
# see: https://docs.sciml.ai/ModelingToolkit/dev/basics/Validation/#Parameter-and-Initial-Condition-Values
function remove_units(p::Dict)
    Dict(k => Unitful.ustrip(ModelingToolkit.get_unit(k), v) for (k, v) in p)
end

add_units(p::Dict) = Dict(k => v * ModelingToolkit.get_unit(k) for (k, v) in p)

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
example_u₀ = Dict(
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

example_p = Dict(
	G => G_val,
	m₁ => 10e26u"kg",
	m₂ => 10e26u"kg"
)

example_tspan = (0.0, 480.0)

example_two_body_problem = ODEProblem(diffeq_two_body_system, remove_units(example_u₀), example_tspan, remove_units(example_p), jac=true)