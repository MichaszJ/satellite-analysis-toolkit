using Plots, ModelingToolkit, DifferentialEquations, Unitful

function remove_units(p::Dict)
    Dict(k => Unitful.ustrip(ModelingToolkit.get_unit(k), v) for (k, v) in p)
end

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
    x₃(t), [unit=u"m"],
    y₃(t), [unit=u"m"],
    z₃(t), [unit=u"m"],
    ẋ₃(t), [unit=u"m/s"],
    ẏ₃(t), [unit=u"m/s"],
    ż₃(t), [unit=u"m/s"],
    r₁₂(t), [unit=u"m"],
    r₁₃(t), [unit=u"m"],
    r₂₃(t), [unit=u"m"]
end)

D = Differential(t)

@parameters G [unit=u"N*m^2/kg^2"] m₁ [unit=u"kg"] m₂ [unit=u"kg"] m₃ [unit=u"kg"]

three_body_equations = [
    r₁₂ ~ sqrt((x₂ - x₁)^2 + (y₂ - y₁)^2 + (z₂ - z₁)^2),
    r₁₃ ~ sqrt((x₃ - x₁)^2 + (y₃ - y₁)^2 + (z₃ - z₁)^2),
    r₂₃ ~ sqrt((x₃ - x₂)^2 + (y₃ - y₂)^2 + (z₃ - z₂)^2),
    
    D(x₁) ~ ẋ₁,
    D(y₁) ~ ẏ₁,
    D(z₁) ~ ż₁,
    D(ẋ₁) ~ G*m₂*(x₂ - x₁)/r₁₂^3 + G*m₃*(x₃ - x₁)/r₁₃^3,
    D(ẏ₁) ~ G*m₂*(y₂ - y₁)/r₁₂^3 + G*m₃*(y₃ - y₁)/r₁₃^3,
    D(ż₁) ~ G*m₂*(z₂ - z₁)/r₁₂^3 + G*m₃*(z₃ - z₁)/r₁₃^3,
    
    D(x₂) ~ ẋ₂,
    D(y₂) ~ ẏ₂,
    D(z₂) ~ ż₂,
    D(ẋ₂) ~ G*m₁*(x₁ - x₂)/r₁₂^3 + G*m₃*(x₃ - x₂)/r₂₃^3,
    D(ẏ₂) ~ G*m₁*(y₁ - y₂)/r₁₂^3 + G*m₃*(y₃ - y₂)/r₂₃^3,
    D(ż₂) ~ G*m₁*(z₁ - z₂)/r₁₂^3 + G*m₃*(z₃ - z₂)/r₂₃^3,
    
    D(x₃) ~ ẋ₃,
    D(y₃) ~ ẏ₃,
    D(z₃) ~ ż₃,
    D(ẋ₃) ~ G*m₁*(x₁ - x₃)/r₁₃^3 + G*m₂*(x₂ - x₃)/r₂₃^3,
    D(ẏ₃) ~ G*m₁*(y₁ - y₃)/r₁₃^3 + G*m₂*(y₂ - y₃)/r₂₃^3,
    D(ż₃) ~ G*m₁*(z₁ - z₃)/r₁₃^3 + G*m₂*(z₂ - z₃)/r₂₃^3,
]

diffeq_three_body_system = ODESystem(
    three_body_equations,
    t,
    name=:three_body_system
) |> structural_simplify

example_u₀ = Dict(
    x₁ => -0.97138u"m",
    y₁ => 0.0u"m",
    z₁ => 0.0u"m",
    ẋ₁ => 0.0u"m/s",
    ẏ₁ => -1.37584u"m/s",
    ż₁ => 0.0u"m/s",
    x₂ => 1.0u"m",
    y₂ => 0.0u"m",
    z₂ => 0.0u"m",
    ẋ₂ => 0.0u"m/s",
    ẏ₂ => -0.34528u"m/s",
    ż₂ => 0.0u"m/s",
    x₃ => 0.0u"m",
    y₃ => 0.0u"m",
    z₃ => 0.0u"m",
    ẋ₃ => 0.0u"m/s",
    ẏ₃ => 1.519362144u"m/s",
    ż₃ => 0.0u"m/s",
)

example_p = Dict(
    G => 1.0u"N*m^2/kg^2",
    m₁ => 0.5312u"kg",
    m₂ => 2.2837u"kg",
    m₃ => 1u"kg"
)

example_three_body_problem = ODEProblem(
    diffeq_three_body_system,
    remove_units(example_u₀),
    (0.0, 10.0),
    remove_units(example_p),
    jac=true
)

example_three_body_sol = solve(example_three_body_problem, Vern7())