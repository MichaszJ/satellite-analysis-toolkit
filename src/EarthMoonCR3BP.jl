using Plots, ModelingToolkit, DifferentialEquations

const mₑ = 5.97e24
const mₘ = 7.3459e22
const G = 6.6743e-11

const r₁₂ = 384400e3

const Ω = sqrt(G * (mₑ + mₘ)/r₁₂^3)
const μ₁ = G*mₑ
const μ₂ = G*mₘ

const π₁ = mₑ/(mₑ + mₘ)
const π₂ = mₘ/(mₑ + mₘ)

@variables(begin
    t,
    x(t),
    y(t),
    z(t),
    ẋ(t),
    ẏ(t),
    ż(t),
    r₁(t),
    r₂(t)
end)

D = Differential(t)

earth_moon_cr_three_body_equations = [
    r₁ ~ sqrt((x + π₂*r₁₂)^2 + y^2 + z^2),
    r₂ ~ sqrt((x - π₁*r₁₂)^2 + y^2 + z^2),
    
    D(x) ~ ẋ,
    D(y) ~ ẏ,
    D(z) ~ ż,
    
    D(ẋ) ~ 2*Ω*ẏ + x*Ω^2 - (x + π₂*r₁₂)*μ₁/r₁^3 - (x - π₁*r₁₂)*μ₂/r₂^3,
    D(ẏ) ~ y*Ω^2 - 2*Ω*ẋ - y*μ₁/r₁^3 - y*μ₂/r₂^3,
    D(ż) ~ -z*μ₁/r₁^3 - z*μ₂/r₂^3,
]

earth_moon_cr_three_body_system = ODESystem(
    earth_moon_cr_three_body_equations,
    t,
    name=:earth_moon_cr_three_body_system
) |> structural_simplify

example_u₀ = Dict(
	x => -4671e3,
    y => -6378e3 - 200e3,
    z => 0.0,
    ẋ => 10.9148e3 * cos(deg2rad(19)),
    ẏ => -10.9148e3 * sin(deg2rad(19)),
    ż => 0.0
)

example_earth_moon_cr_three_body_problem = ODEProblem(
    earth_moon_cr_three_body_system,
    example_u₀,
    (0.0, 3.4 * 24 * 60 * 60),
    [],
    jac=true
)

example_cr_three_body_sol = solve(example_earth_moon_cr_three_body_problem, Tsit5())