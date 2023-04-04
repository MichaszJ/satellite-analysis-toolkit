include("../src/SAT.jl");

using GLMakie, GeoMakie, FileIO, .SAT
using Downloads: download
GLMakie.activate!()

Makie.set_theme!(theme_dark())

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

    # a, e, i, asc, arg, theta = orbital_elements
    elements = Observable([
        8350,
        0.19760,
        deg2rad(60),
        deg2rad(270),
        deg2rad(45),
        deg2rad(230)
    ]);

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

    lon, lat, alt, times = ground_track(to_value(elements), 0.0, to_value(final_time_obsv)*60^2, num_steps=to_value(num_steps_obsv), return_times=true);

    ground_track_coords = Observable(Matrix(hcat(clamp_lon.(rad2deg.(Float64.(lon))), clamp_lat.(rad2deg.(Float64.(lat))))'));

    const eq_rad = 6378.137 
    const ecc_sq = 6.69437999014e-3
    N = eq_rad ./ sqrt.(1 .- ecc_sq .* sin.(lat).^2);

    pos_x = (N .+ alt) .* cos.(lat) .* cos.(lon);
    pos_y = (N .+ alt) .* cos.(lat) .* sin.(lon);
    pos_z = ((1 - ecc_sq) .* N .+ alt) .* sin.(lat);

    pos_coords = Observable(Matrix(hcat(pos_x, pos_y, pos_z)'));

    scatter!(ax2, ground_track_coords, marker=:cross, color=:seagreen2)
    lines!(ax1, pos_coords, linestyle = :dash, color=:seagreen2)

    on(num_steps_sl.value) do ns
        lon, lat, alt, times = ground_track(to_value(elements), 0.0, to_value(final_time_obsv)*60^2, num_steps=ns, return_times=true)
        ground_track_coords[] = Matrix(hcat(clamp_lon.(rad2deg.(Float64.(lon))), clamp_lat.(rad2deg.(Float64.(lat))))')

        N = eq_rad ./ sqrt.(1 .- ecc_sq .* sin.(lat).^2)

        pos_x = (N .+ alt) .* cos.(lat) .* cos.(lon)
        pos_y = (N .+ alt) .* cos.(lat) .* sin.(lon)
        pos_z = ((1 - ecc_sq) .* N .+ alt) .* sin.(lat)

        pos_coords[] = Matrix(hcat(pos_x, pos_y, pos_z)')
        autolimits!(ax1)
    end

    on(final_time_sl.value) do tf
        lon, lat, alt, times = ground_track(to_value(elements), 0.0, tf*60^2, num_steps=to_value(num_steps_obsv), return_times=true)
        ground_track_coords[] = Matrix(hcat(clamp_lon.(rad2deg.(Float64.(lon))), clamp_lat.(rad2deg.(Float64.(lat))))')

        N = eq_rad ./ sqrt.(1 .- ecc_sq .* sin.(lat).^2)

        pos_x = (N .+ alt) .* cos.(lat) .* cos.(lon)
        pos_y = (N .+ alt) .* cos.(lat) .* sin.(lon)
        pos_z = ((1 - ecc_sq) .* N .+ alt) .* sin.(lat)

        pos_coords[] = Matrix(hcat(pos_x, pos_y, pos_z)')
        autolimits!(ax1)
    end

    ax4 = fig[2,2] = GridLayout()

    elements_sg = SliderGrid(
        ax4[1,1],
        (label = "Semi-Major Axis", range = 7000:45000, format = "{:.2f} km", startvalue = 8350),
        (label = "Eccentricity", range = 0:0.001:0.999, format = "{:.3f}", startvalue = 0.19760),
        (label = "Inclination", range = 0:180, format = "{:.2}°", startvalue = 60),
        (label = "RAAN", range = 0:360, format = "{:.2}°", startvalue = 270),
        (label = "Argument of Periapsis", range = 0:360, format = "{:.2}°", startvalue = 45),
        (label = "True Anomaly", range = 0:360, format = "{:.2}°", startvalue = 230),
    )

    sliderobservables = [s.value for s in elements_sg.sliders]
    bars = lift(sliderobservables...) do slvalues...
        a, e, i, asc, arg, theta  = [slvalues...]

        elements[] = [a, e, deg2rad(i), deg2rad(asc), deg2rad(arg), deg2rad(theta)]

        lon, lat, alt, times = ground_track(to_value(elements), 0.0, to_value(final_time_obsv)*60^2, num_steps=to_value(num_steps_obsv), return_times=true)
        ground_track_coords[] = Matrix(hcat(clamp_lon.(rad2deg.(Float64.(lon))), clamp_lat.(rad2deg.(Float64.(lat))))')

        N = eq_rad ./ sqrt.(1 .- ecc_sq .* sin.(lat).^2)

        pos_x = (N .+ alt) .* cos.(lat) .* cos.(lon)
        pos_y = (N .+ alt) .* cos.(lat) .* sin.(lon)
        pos_z = ((1 - ecc_sq) .* N .+ alt) .* sin.(lat)

        pos_coords[] = Matrix(hcat(pos_x, pos_y, pos_z)')
        autolimits!(ax1)

        [slvalues...]
    end
end