module SAT

include("GNC.jl");
export determine_eccentric_anomaly, geocentric_state_vector_transform, position_to_asc_dec, ground_track

include("Satellite.jl");
export Satellite

# include("TwoBody.jl");
# export remove_units, add_units, G_val, t, x₁, y₁, z₁, ẋ₁, ẏ₁, ż₁, x₂, y₂, z₂, ẋ₂, ẏ₂, ż₂, r
# export D, G, m₁, m₂, two_body_equations, diffeq_two_body_system
# export example_u₀, example_p, example_tspan, example_two_body_problem

end