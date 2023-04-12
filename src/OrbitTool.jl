using GLMakie, GeoMakie, ModelingToolkit, DifferentialEquations, Unitful
GLMakie.activate!()

Makie.set_theme!(theme_dark())

include("TwoBody.jl")
# using .TwoBody

two_body_system = TwoBodySystem(two_body_example_u₀, two_body_example_p, two_body_example_tspan)
two_body_sol = SolveOrbitalSystem(two_body_system, Tsit5());

times = collect(0.0:0.1:480.0);
sol_interp = two_body_sol(times);

fig = Figure(resolution=(1500,1500)); display(fig);
ax1 = Axis3(fig[1,1])

# lines1 = Matrix(hcat(sol_interp[x₁], sol_interp[y₁], sol_interp[z₁])');
lines!(ax1, sol_interp[x₁], sol_interp[y₁], sol_interp[z₁], label="Mass 1")
lines!(ax1, sol_interp[x₂], sol_interp[y₂], sol_interp[z₂], label="Mass 2")

axislegend(ax1)

ax2 = fig[1,2] = GridLayout()

colsize!(fig.layout, 2, Relative(1/4))

Textbox(ax2[1, 1], placeholder = "Enter a string...")
