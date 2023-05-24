include("SAT.jl");

using GLMakie, GeoMakie, FileIO, LinearAlgebra, .SAT
using Downloads: download
GLMakie.activate!()

Makie.set_theme!(theme_dark())

using Colors, ColorSchemes

# a, e, i, asc, arg, theta = orbital_elements
# mutable struct Elements{T<:Number}
#     a::Union{T, Quantity{T}}
#     e::T
#     i::Union{T, Quantity{T}}
#     Ω::Union{T, Quantity{T}}
#     ω::Union{T, Quantity{T}}
#     θ::Union{T, Quantity{T}}
# end

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
    "Sat 1" => GroundTrackerSatellite("Sat 1", Observable([10000, 0.19760, deg2rad(60), deg2rad(270), deg2rad(45), deg2rad(230)])),
    "Sat 2" => GroundTrackerSatellite("Sat 2", Observable([8350, 0.05, deg2rad(60), deg2rad(270), deg2rad(45), deg2rad(50)])),
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

begin
    # source: https://beautiful.makie.org/dev/examples/generated/2d/geo/blue_marble/
    earth_img = load(download("https://upload.wikimedia.org/wikipedia/commons/5/56/Blue_Marble_Next_Generation_%2B_topography_%2B_bathymetry.jpg"));
    n = 1024 ÷ 4 # 2048
    θ = LinRange(0, π, n);
    φ = LinRange(0, 2π, 2 * n);
    x = [cos(φ) * sin(θ) for θ in θ, φ in φ] .* 6378.137;
    y = [sin(φ) * sin(θ) for θ in θ, φ in φ] .* 6378.137;
    z = [cos(θ) for θ in θ, φ in φ] .* 6378.137;

    fig = Figure(resolution=(2000,1500)); display(fig);
    ax1 = Axis3(fig[1,1], aspect=:data, viewmode=:fitzoom, title="Satellite ECI Frame Orbit")
    surf = surface!(
        ax1, x, y, z;
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



    # clamp_lon(lon) = lon > 180.0 || lon < -180.0 ? lon % 180 - sign(lon)*floor(abs(lon/180))*180 : lon
    # clamp_lat(lat) = lat > 90.0 || lat < -90.0 ? lat % 90 - sign(lat)*floor(abs(lat/90))*90 : lat

    # make more elegant
    function clamp_lon(lon)
        if lon > 180.0
            while lon > 180.0
                lon -= 360
            end
            return lon
        elseif lon < -180.0
            while lon < 180.0
                lon += 360
            end
            return lon
        else
            return lon
        end
    end

    function clamp_lat(lat)
        if lat > 90.0
            while lat > 90.0
                lat -= 180
            end
            return lat
        elseif lat < -90.0
            while lat < 90.0
                lat += 180
            end
            return lat
        else
            return lat
        end
    end

    const eq_rad = 6378.137 
    const ecc_sq = 6.69437999014e-3

    # gt_coords_dict = Dict()
    # pos_coords_dict = Dict()

    # for (key, sat) in sats

    orbit_lines, geo_points = Dict(), Dict()

    for (key, sat) in sats
        lon, lat, alt, times = ground_track(to_value(sat.elements), 0.0, to_value(final_time_obsv)*60^2, num_steps=to_value(num_steps_obsv), return_times=true);

        sat.ground_track_coords = Observable(Matrix(hcat(clamp_lon.(rad2deg.(Float64.(lon))), clamp_lat.(rad2deg.(Float64.(lat))))'));
    
        N = eq_rad ./ sqrt.(1 .- ecc_sq .* sin.(lat).^2);
    
        pos_x = (N .+ alt) .* cos.(lat) .* cos.(lon);
        pos_y = (N .+ alt) .* cos.(lat) .* sin.(lon);
        pos_z = ((1 - ecc_sq) .* N .+ alt) .* sin.(lat);
    
        sat.pos_coords = Observable(Matrix(hcat(pos_x, pos_y, pos_z)'));

        pos = [[pos_x[i], pos_y[i], pos_z[i]] for i in eachindex(pos_x)]
        visibility = Observable([check_gs_visibility(ground_stations["GS 1"], p) for p in pos])

        sat_visibilities[key] = visibility
    
        # push!(geo_points, scatter!(ax2, sat.ground_track_coords, marker=:cross, label=key, color=visibility))
        # push!(orbit_lines, lines!(ax1, sat.pos_coords, linestyle=:dash, label=key, color=visibility))

        orbit_lines[key] = lines!(ax1, sat.pos_coords, linestyle=:dash, color=visibility, colormap=colormaps[key])
        geo_points[key] = scatter!(ax2, sat.ground_track_coords, marker=:cross, color=visibility, colormap=colormaps[key])
    end

    # ground station plotting
    for (key, station) in ground_stations
        longitude, latitude = position_to_asc_dec(station.position)
        scatter!(ax2, rad2deg.([longitude latitude]), marker='⧇', markersize=40, label=key)
        # meshscatter!(ax1, station.position..., markersize=1000, label=key)
    end

    # axislegend(ax1)
    axislegend(ax2)

    colorbar_ax = fig[2,3] = GridLayout()    
    # colsize!(fig.layout, 3, Relative(0.175))

    for (i, (key, sat)) in enumerate(sats)
        if i == 1
            Colorbar(colorbar_ax[i,1], orbit_lines[key], label=key, vertical = false, ticks=([0.17, 0.5, 0.82], ["Not Visible", "Default", "Visble"]))
        else
            Colorbar(colorbar_ax[i,1], orbit_lines[key], label=key, vertical = false, ticklabelsvisible=false)
        end
    end

    on(num_steps_sl.value) do ns
        # for (key, sat) in sats
        for (key, sat) in sats
            lon, lat, alt, times = ground_track(to_value(sat.elements), 0.0, to_value(final_time_obsv)*60^2, num_steps=ns, return_times=true)
            sat.ground_track_coords[] = Matrix(hcat(clamp_lon.(rad2deg.(Float64.(lon))), clamp_lat.(rad2deg.(Float64.(lat))))')

            N = eq_rad ./ sqrt.(1 .- ecc_sq .* sin.(lat).^2)

            pos_x = (N .+ alt) .* cos.(lat) .* cos.(lon)
            pos_y = (N .+ alt) .* cos.(lat) .* sin.(lon)
            pos_z = ((1 - ecc_sq) .* N .+ alt) .* sin.(lat)

            pos = [[pos_x[i], pos_y[i], pos_z[i]] for i in eachindex(pos_x)]
            sat_visibilities[key][] = [check_gs_visibility(ground_stations["GS 1"], p) for p in pos]

            sat.pos_coords[] = Matrix(hcat(pos_x, pos_y, pos_z)')
        end
        autolimits!(ax1)
    end

    on(final_time_sl.value) do tf
        # for (key, sat) in sats
        for (key, sat) in sats
            lon, lat, alt, times = ground_track(to_value(sat.elements), 0.0, tf*60^2, num_steps=to_value(num_steps_obsv), return_times=true)
            sat.ground_track_coords[] = Matrix(hcat(clamp_lon.(rad2deg.(Float64.(lon))), clamp_lat.(rad2deg.(Float64.(lat))))')

            N = eq_rad ./ sqrt.(1 .- ecc_sq .* sin.(lat).^2)

            pos_x = (N .+ alt) .* cos.(lat) .* cos.(lon)
            pos_y = (N .+ alt) .* cos.(lat) .* sin.(lon)
            pos_z = ((1 - ecc_sq) .* N .+ alt) .* sin.(lat)

            pos = [[pos_x[i], pos_y[i], pos_z[i]] for i in eachindex(pos_x)]
            sat_visibilities[key][] = [check_gs_visibility(ground_stations["GS 1"], p) for p in pos]

            sat.pos_coords[] = Matrix(hcat(pos_x, pos_y, pos_z)')
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

        lon, lat, alt, times = ground_track(to_value(active_sat.elements), 0.0, to_value(final_time_obsv)*60^2, num_steps=to_value(num_steps_obsv), return_times=true)
        active_sat.ground_track_coords[] = Matrix(hcat(clamp_lon.(rad2deg.(Float64.(lon))), clamp_lat.(rad2deg.(Float64.(lat))))')

        N = eq_rad ./ sqrt.(1 .- ecc_sq .* sin.(lat).^2)

        pos_x = (N .+ alt) .* cos.(lat) .* cos.(lon)
        pos_y = (N .+ alt) .* cos.(lat) .* sin.(lon)
        pos_z = ((1 - ecc_sq) .* N .+ alt) .* sin.(lat)

        pos = [[pos_x[i], pos_y[i], pos_z[i]] for i in eachindex(pos_x)]
        sat_visibilities[active_sat.name][] = [check_gs_visibility(ground_stations["GS 1"], p) for p in pos]

        active_sat.pos_coords[] = Matrix(hcat(pos_x, pos_y, pos_z)')
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