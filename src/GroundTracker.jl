include("SAT.jl");

using GLMakie, GeoMakie, FileIO, LinearAlgebra, .SAT
using Downloads: download
GLMakie.activate!()

Makie.set_theme!(theme_dark())
update_theme!(fonts = (; regular = "iA Writer Quattro V", bold = "iA Writer Quattro V"))

using Colors, ColorSchemes

using ModelingToolkit, DifferentialEquations

mutable struct GroundTrackerSatellite
    name::String
    elements::Observable{Vector{Float64}}
    ground_track_coords::Union{Observable{Matrix{Float64}}, Nothing}
    pos_coords::Union{Observable{Matrix{Float64}}, Nothing}

    function GroundTrackerSatellite(name::String, elements::Observable{Vector{Float64}})
        return new(name, elements, nothing, nothing)
    end
end

# a, e, i, asc, arg, theta = orbital_elements
sats = Dict(
    "Sat 1" => GroundTrackerSatellite("Sat 1", Observable([8350, 0.19760, deg2rad(60), deg2rad(270), deg2rad(45), deg2rad(230)])),
    "Sat 2" => GroundTrackerSatellite("Sat 2", Observable([15000, 0.05, deg2rad(60), deg2rad(270), deg2rad(45), deg2rad(50)])),
)

sat_visibilities = Dict()

struct GroundStation
	name::String
	position::Vector{Float64}
	min_elevation::Float64
end

ground_stations = Dict(
    "GS 1" => GroundStation("GS 1", [2761.8, 4783.5, 3189.0], deg2rad(5.0))
)

function check_gs_visibility(station::GroundStation, sat_position::Vector{Float64})
    slant = sat_position - station.position
    slant_range = norm(slant)
    local_vert = station.position ./ norm(station.position)
    elevation = π/2 - acos(dot(local_vert, slant) / slant_range)

    return elevation >= station.min_elevation
end

colormaps = Dict(
    # "Sat 1" => range(colorant"red", stop=colorant"green", length=2),
    # "Sat 2" => range(colorant"blue", stop=colorant"yellow", length=2)
    # "Sat 1" => colormap("Blues"),
    # "Sat 2" => colormap("Greens"),
    # "Sat 1" => colormap("Greens"),
    # "Sat 2" => cgrad(ColorScheme([RGB{Float64}(i, 1.5i, 2i) for i in [0.25, 0.35, 0.5]]), 3, categorical=true)

    "Sat 1" => cgrad(ColorScheme([colorant"#7E1717", colorant"#068DA9", colorant"#E55807"]), 3, categorical=true),
    "Sat 2" => cgrad(ColorScheme([colorant"#146C94", colorant"#AFD3E2", colorant"#19A7CE"]), 3, categorical=true)
)

ωe = 0.261799387799149/(60*60)

J2 = 1.083e-3
Re = 6378.137
μ = 398600
f = 1/298.257223563

@variables t x(t) y(t) z(t) ẋ(t) ẏ(t) ż(t) r(t) E(t) F(t) G(t)
D = Differential(t)

eqs = [
    r ~ sqrt(x^2 + y^2 + z^2),
    
    D(x) ~ ẋ,
    D(y) ~ ẏ,
    D(z) ~ ż,

    D(ẋ) ~ -x*μ*(2*r^2 + 3*J2*Re^2*(5*z^2/r^2 - 1)) / (2*r^5),
    D(ẏ) ~ -y*μ*(2*r^2 + 3*J2*Re^2*(5*z^2/r^2 - 1)) / (2*r^5),
    D(ż) ~ -z*μ*(2*r^4 + 9*J2*r^2 * Re^2 - 15*J2*Re^2 * z^2) / (2*r^7),

    E ~ x*cos(ωe*t) - y*sin(ωe*t),
    F ~ y*cos(ωe*t) + x*sin(ωe*t),
    G ~ z
]
    
J2_system = structural_simplify(ODESystem(
    eqs,
    t,
    name=:J2_system
))

function elements_to_state(a, e, i, Ω, ω, θ; μ=398600)
    h = sqrt(-a*μ*(e^2 - 1))

    r_w = h^2 / μ / (1 + e * cos(θ)) .* [cos(θ) sin(θ) 0];
    v_w = μ / h .* [-sin(θ) e + cos(θ) 0];

    R1 = [cos(-ω) -sin(-ω) 0; sin(-ω) cos(-ω) 0; 0 0 1];
    R2 = [1 0 0; 0 cos(-i) -sin(-i); 0 sin(-i) cos(-i)];
    R3 = [cos(-Ω) -sin(-Ω) 0; sin(-Ω) cos(-Ω) 0; 0 0 1];
    r_rot = r_w * R1 * R2 * R3;
    v_rot = v_w * R1 * R2 * R3;
    
    return vec(r_rot), vec(v_rot)
end

const a = 6378.137
const e = sqrt(2*fe - fe^2)

function LLA_from_EFG_approx(E, F, G)
    ϕ = atan(G, sqrt(E^2 + F^2)*(1 - f)^2)
    λ = atan(F, E)
    h = sqrt(E^2 + F^2 + G^2) - a*sqrt((1 - e^2)/(1 - e^2 * cos(atan(G, sqrt(E^2 + F^2)))^2))

    return [ϕ, λ, h]
end

function normal_distance(ϕ)
    a / sqrt(1 - (e * sin(ϕ))^2)
end

function EFG_from_LLA(ϕ, λ, h)
    E = (normal_distance(ϕ) + h) * cos(λ) * cos(ϕ)
    F = (normal_distance(ϕ) + h) * sin(λ) * cos(ϕ)
    G = (normal_distance(ϕ) * (1 - e^2) + h)*sin(ϕ)

    return [E, F, G]
end

function LLA_from_EFG(E, F, G; tol=0.01)
    pos = [E, F, G]
    fakePos = copy(pos)
    result = LLA_from_EFG_approx(fakePos...)
    error = EFG_from_LLA(result...) .- pos

    while norm(error) > tol
        fakePos .-= error
        result = LLA_from_EFG_approx(fakePos...)
        error = EFG_from_LLA(result...) .- pos
    end

    return reverse(result[1:2])
end

function GroundTrack2(elements, tspan, num_steps)
    states = elements_to_state(elements...)

    u0 = Dict(
        x => states[1][1], 
        y => states[1][2], 
        z => states[1][3], 
        ẋ => states[2][1], 
        ẏ => states[2][2], 
        ż => states[2][3]
    )

    prob = ODEProblem(
        J2_system, 
        u0, 
        tspan, 
        [],
        jac=true
    )

    sol = solve(prob)

    times = LinRange(0, tspan[end], num_steps)
    interp = sol(times)

    ground_track_coords = rad2deg.(reduce(hcat,  LLA_from_EFG.(interp[E], interp[F], interp[G])))
    orbit_coords = Matrix(hcat(interp[E], interp[F], interp[G])')

    return ground_track_coords, orbit_coords
end

begin
    # source: https://beautiful.makie.org/dev/examples/generated/2d/geo/blue_marble/
    earth_img = load(download("https://upload.wikimedia.org/wikipedia/commons/5/56/Blue_Marble_Next_Generation_%2B_topography_%2B_bathymetry.jpg"));
    n = 1024 ÷ 4 # 2048
    θe = LinRange(0, π, n);
    φe = LinRange(0, 2π, 2 * n);
    xe = [cos(φe) * sin(θe) for θe in θe, φe in φe] .* 6378.137;
    ye = [sin(φe) * sin(θe) for θe in θe, φe in φe] .* 6378.137;
    ze = [cos(θe) for θe in θe, φe in φe] .* 6378.137;

    fig = Figure(resolution=(2000,1500), fonts = (; regular = "iA Writer Quattro V")); display(fig);
    ax1 = Axis3(fig[1,1], aspect=:data, viewmode=:fitzoom, title="Satellite ECEF Frame Orbit")
    surf = surface!(
        ax1, xe, ye, ze;
        color = earth_img,
        shading = false,
        lightposition = Vec3f(-2, -3, -3),
        ambient = Vec3f(0.8, 0.8, 0.8),
        backlight = 1.5f0,
    )
    Makie.rotate!(surf, Quaternionf(0.0, 0.0, 1, 0));

    colsize!(fig.layout, 1, Relative(1/3))

    ax3 = fig[2,1] = GridLayout()
    rowsize!(fig.layout, 2, Relative(1/5))

    num_steps_obsv = Observable(250)
    num_steps_sl = Slider(ax3[1,2], range = 50:750, startvalue = to_value(num_steps_obsv), horizontal = true)
    connect!(num_steps_obsv, num_steps_sl.value)
    Label(ax3[1,1], @lift("Number of Computed Points: $($num_steps_obsv)"))

    final_time_obsv = Observable(24) # in hours
    final_time_sl = Slider(ax3[2,2], range = 1:50, startvalue = to_value(final_time_obsv), horizontal = true)
    connect!(final_time_obsv, final_time_sl.value)
    Label(ax3[2,1], @lift("Orbit Time: $($final_time_obsv) hr"))

    projections = [
        "+proj=natearth", "+proj=natearth2", "+proj=moll", "+proj=weren", "+proj=times", "+proj=wink1", "+proj=wink2","+proj=sinu", "+proj=crast",
    ]

    proj_obsv = Observable("")
    Label(ax3[3,1], "Map Projection")
    proj_menu = Menu(ax3[3,2], options=projections)
    connect!(proj_obsv, proj_menu.selection)

    on(proj_menu.selection) do s
        proj_obsv[] = s
    end

    proj = lift(proj_obsv) do s
        s
    end

    ax2 = GeoAxis(
        fig[1, 2:3]; 
        title="Satellite Groundtrack",
        coastlines = true,
        coastline_attributes = (; color = "white"),
        dest=proj,
    );

    orbit_lines, geo_points = Dict(), Dict()

    for (key, sat) in sats
        ground_track_coords, orbit_coords =  GroundTrack2(to_value(sat.elements), [0.0, to_value(final_time_obsv)*60^2], to_value(num_steps_obsv))
        sat.ground_track_coords = Observable(ground_track_coords);
    
        sat.pos_coords = Observable(orbit_coords);

        visibility = Observable([check_gs_visibility(ground_stations["GS 1"], p) for p in [[orbit_coords[1,i], orbit_coords[2,i], orbit_coords[3,i]] for i in 1:size(orbit_coords, 2)]])

        sat_visibilities[key] = visibility

        orbit_lines[key] = lines!(ax1, sat.pos_coords, linestyle=:dash, color=visibility, colormap=colormaps[key])
        geo_points[key] = scatter!(ax2, sat.ground_track_coords, marker=:cross, color=visibility, colormap=colormaps[key])
    end

    for (key, station) in ground_stations
        longitude, latitude = position_to_asc_dec(station.position)
        scatter!(ax2, rad2deg.([longitude latitude]), marker='⧇', markersize=40, label=key)
    end

    axislegend(ax2)

    colorbar_ax = fig[2,3] = GridLayout()    

    for (i, (key, sat)) in enumerate(sats)
        if i == 1
            Colorbar(colorbar_ax[i,1], orbit_lines[key], label=key, vertical = false, ticks=([0.17, 0.5, 0.82], ["Not Visible", "Default", "Visble"]))
        else
            Colorbar(colorbar_ax[i,1], orbit_lines[key], label=key, vertical = false, ticklabelsvisible=false)
        end
    end

    on(num_steps_sl.value) do ns
        for (key, sat) in sats
            ground_track_coords, orbit_coords =  GroundTrack2(to_value(sat.elements), [0.0, to_value(final_time_obsv)*60^2], to_value(num_steps_obsv))
            sat.ground_track_coords[] = ground_track_coords;
        
            sat.pos_coords[] = orbit_coords

            sat_visibilities[key][] = [check_gs_visibility(ground_stations["GS 1"], p) for p in [[orbit_coords[1,i], orbit_coords[2,i], orbit_coords[3,i]] for i in 1:size(orbit_coords,2)]]
        end
        autolimits!(ax1)
    end

    on(final_time_sl.value) do tf
        for (key, sat) in sats
            ground_track_coords, orbit_coords =  GroundTrack2(to_value(sat.elements), [0.0, to_value(final_time_obsv)*60^2], to_value(num_steps_obsv))
            sat.ground_track_coords[] = ground_track_coords;
        
            sat.pos_coords[] = orbit_coords

            sat_visibilities[key][] = [check_gs_visibility(ground_stations["GS 1"], p) for p in [[orbit_coords[1,i], orbit_coords[2,i], orbit_coords[3,i]] for i in 1:size(orbit_coords,2)]]
        end
        autolimits!(ax1)
    end

    rowsize!(fig.layout, 1, Relative(0.6))

    sats_list = [key for (key, sat) in sats]
    selected_sat_obsv = Observable(sats_list[1])

    Label(ax3[4,1], "Select Satellite")
    sat_menu = Menu(ax3[4, 2], options=sats_list,)

    on(sat_menu.selection) do s
        selected_sat_obsv[] = s

        active_sat_elements = to_value(sats[s].elements)

        a_textbox.stored_string[] = active_sat_elements[1] |> string
        e_textbox.stored_string[] = active_sat_elements[2] |> string
        i_textbox.stored_string[] = rad2deg(active_sat_elements[3]) |> string
        Ω_textbox.stored_string[] = rad2deg(active_sat_elements[4]) |> string
        ω_textbox.stored_string[] = rad2deg(active_sat_elements[5]) |> string
        θ_textbox.stored_string[] = rad2deg(active_sat_elements[6]) |> string

        a_textbox.displayed_string[] = active_sat_elements[1] |> string
        e_textbox.displayed_string[] = active_sat_elements[2] |> string
        i_textbox.displayed_string[] = rad2deg(active_sat_elements[3]) |> string
        Ω_textbox.displayed_string[] = rad2deg(active_sat_elements[4]) |> string
        ω_textbox.displayed_string[] = rad2deg(active_sat_elements[5]) |> string
        θ_textbox.displayed_string[] = rad2deg(active_sat_elements[6]) |> string
    end

    ax4 = fig[2,2] = GridLayout()    
    
    row1 = ax4[1,1] = GridLayout()
    Label(row1[1, 1], "Semi-Major Axis", tellwidth=false, halign=:left)
    a_validator(str) = tryparse(Float64, str) ≠ nothing && (7000 ≤ parse(Float64, str) ≤ 45000) ? true : false
    a_textbox = Textbox(row1[1, 2], validator=a_validator, halign=:left, tellwidth=false, width=250, stored_string=to_value(sats[sats_list[1]].elements)[1] |> string)

    row2 = ax4[2,1] = GridLayout()
    Label(row2[1, 1], "Eccentricity", tellwidth=false, halign=:left)
    e_validator(str) = tryparse(Float64, str) ≠ nothing && (0 ≤ parse(Float64, str) < 1) ? true : false
    e_textbox = Textbox(row2[1, 2], validator=e_validator, halign=:left, tellwidth=false, width=250, stored_string=to_value(sats[sats_list[1]].elements)[2] |> string)

    row3 = ax4[3,1] = GridLayout()
    Label(row3[1, 1], "Inclination", tellwidth=false, halign=:left)
    i_validator(str) = tryparse(Float64, str) ≠ nothing && (0 ≤ parse(Float64, str) ≤ 180) ? true : false
    i_textbox = Textbox(row3[1, 2], validator=i_validator, halign=:left, tellwidth=false, width=250, stored_string=to_value(sats[sats_list[1]].elements)[3] |> rad2deg |> string)

    row4 = ax4[4,1] = GridLayout()
    Label(row4[1, 1], "RAAN", tellwidth=false, halign=:left)
    degree_validator(str) = tryparse(Float64, str) ≠ nothing && (0 ≤ parse(Float64, str) ≤ 360) ? true : false
    Ω_textbox = Textbox(row4[1, 2], validator=degree_validator, halign=:left, tellwidth=false, width=250, stored_string=to_value(sats[sats_list[1]].elements)[4] |> rad2deg |> string)

    row5 = ax4[5,1] = GridLayout()
    Label(row5[1, 1], "Argument of Periapsis", tellwidth=false, halign=:left)
    ω_textbox = Textbox(row5[1, 2], validator=degree_validator, halign=:left, tellwidth=false, width=250, stored_string=to_value(sats[sats_list[1]].elements)[5] |> rad2deg |> string)

    row6 = ax4[6,1] = GridLayout()
    Label(row6[1, 1], "True Anomaly", tellwidth=false, halign=:left)
    θ_textbox = Textbox(row6[1, 2], validator=degree_validator, halign=:left, tellwidth=false, width=250, stored_string=to_value(sats[sats_list[1]].elements)[6] |> rad2deg |> string)

    textboxes = lift(a_textbox.stored_string, e_textbox.stored_string, i_textbox.stored_string, Ω_textbox.stored_string, ω_textbox.stored_string, θ_textbox.stored_string) do slvalues...
        a, e, i, asc, arg, theta  = parse.(Float64, [slvalues...])

        active_sat = sats[to_value(selected_sat_obsv)]

        active_sat.elements[] = [a, e, deg2rad(i), deg2rad(asc), deg2rad(arg), deg2rad(theta)]

        ground_track_coords, orbit_coords =  GroundTrack2(to_value(active_sat.elements), [0.0, to_value(final_time_obsv)*60^2], to_value(num_steps_obsv))
        active_sat.ground_track_coords[] = ground_track_coords;
    
        active_sat.pos_coords[] = orbit_coords

        sat_visibilities[active_sat.name][] = [check_gs_visibility(ground_stations["GS 1"], p) for p in [[orbit_coords[1,i], orbit_coords[2,i], orbit_coords[3,i]] for i in 1:size(orbit_coords,2)]]
        autolimits!(ax1)

        [slvalues...]
    end

    Label(ax3[5,1], "Show GS Visibility")
    visibility_toggle = Toggle(ax3[5,2], active = true, tellwidth=false)

    on(visibility_toggle.active) do state
        if state == true
            for (key, lines) in orbit_lines
                lines.color = to_value(sat_visibilities[key])
            end

            for (key, points) in geo_points
                points.color = to_value(sat_visibilities[key])
            end
        else
            for (key, lines) in orbit_lines
                lines.color = [0.5 for _ in 1:length(to_value(sat_visibilities[key]))]
            end

            for (key, points) in geo_points
                points.color = [0.5 for _ in 1:length(to_value(sat_visibilities[key]))]
            end
        end
    end
end