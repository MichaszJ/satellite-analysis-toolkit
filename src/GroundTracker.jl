include("SAT.jl");

using GLMakie, GeoMakie, FileIO, LinearAlgebra, .SAT
using Downloads: download
GLMakie.activate!()

Makie.set_theme!(theme_dark())

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

begin
    # source: https://beautiful.makie.org/dev/examples/generated/2d/geo/blue_marble/
    earth_img = load(download("https://upload.wikimedia.org/wikipedia/commons/5/56/Blue_Marble_Next_Generation_%2B_topography_%2B_bathymetry.jpg"));
    n = 1024 ÷ 4 # 2048
    θ = LinRange(0, π, n);
    φ = LinRange(0, 2π, 2 * n);
    x = [cos(φ) * sin(θ) for θ in θ, φ in φ] .* 6378.137;
    y = [sin(φ) * sin(θ) for θ in θ, φ in φ] .* 6378.137;
    z = [cos(θ) for θ in θ, φ in φ] .* 6378.137;

    fig = Figure(resolution=(2000,1000)); display(fig);
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
    Label(ax3[1,3], "Map Projection")
    proj_menu = Menu(ax3[2,3], options=projections)
    connect!(proj_obsv, proj_menu.selection)

    on(proj_menu.selection) do s
        proj_obsv[] = s
    end

    proj = lift(proj_obsv) do s
        s
    end

    ax2 = GeoAxis(
        fig[1, 2]; 
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
    
        scatter!(ax2, sat.ground_track_coords, marker=:cross, label=key, color=visibility)
        lines!(ax1, sat.pos_coords, linestyle=:dash, label=key, color=visibility)
    end

    # ground station plotting
    for (key, station) in ground_stations
        longitude, latitude = position_to_asc_dec(station.position)
        scatter!(ax2, rad2deg.([longitude latitude]), marker='⧇', markersize=40, label=key)
        # meshscatter!(ax1, station.position..., markersize=1000, label=key)
    end

    axislegend(ax1)
    axislegend(ax2)

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

    ax4 = fig[2,2] = GridLayout()

    sats_list = [key for (key, sat) in sats]

    elements_sg = SliderGrid(
        ax4[2,1],
        (label = "Semi-Major Axis", range = 7000:45000, format = "{:.2f} km", startvalue = to_value(sats[sats_list[1]].elements)[1]),
        (label = "Eccentricity", range = 0:0.001:0.999, format = "{:.3f}", startvalue = to_value(sats[sats_list[1]].elements)[2]),
        (label = "Inclination", range = 0:180, format = "{:.2}°", startvalue = rad2deg(to_value(sats[sats_list[1]].elements)[3])),
        (label = "RAAN", range = 0:360, format = "{:.2}°", startvalue = rad2deg(to_value(sats[sats_list[1]].elements)[4])),
        (label = "Argument of Periapsis", range = 0:360, format = "{:.2}°", startvalue = rad2deg(to_value(sats[sats_list[1]].elements)[5])),
        (label = "True Anomaly", range = 0:360, format = "{:.2}°", startvalue = rad2deg(to_value(sats[sats_list[1]].elements)[6])),
    )

    sliderobservables = [s.value for s in elements_sg.sliders]
    
    selected_sat_obsv = Observable(sats_list[1])
    sat_menu = Menu(ax4[1, 1], options=sats_list)

    on(sat_menu.selection) do s
        selected_sat_obsv[] = s

        elements_sg.sliders[1].value[] = to_value(sats[s].elements)[1]
        elements_sg.sliders[2].value[] = to_value(sats[s].elements)[2]
        elements_sg.sliders[3].value[] = rad2deg(to_value(sats[s].elements)[3])
        elements_sg.sliders[4].value[] = rad2deg(to_value(sats[s].elements)[4])
        elements_sg.sliders[5].value[] = rad2deg(to_value(sats[s].elements)[5])
        elements_sg.sliders[6].value[] = rad2deg(to_value(sats[s].elements)[6])

        sliderobservables = [slider.value for slider in elements_sg.sliders]
    end

    bars = lift(sliderobservables...) do slvalues...
        a, e, i, asc, arg, theta  = [slvalues...]

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
end