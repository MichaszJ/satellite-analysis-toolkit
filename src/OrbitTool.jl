using GLMakie, GeoMakie, FileIO, .SAT
using Downloads: download
GLMakie.activate!()

Makie.set_theme!(theme_dark())

earth_img = load(download("https://upload.wikimedia.org/wikipedia/commons/5/56/Blue_Marble_Next_Generation_%2B_topography_%2B_bathymetry.jpg"));
n = 1024 ÷ 4 # 2048
θ = LinRange(0, π, n);
φ = LinRange(0, 2π, 2 * n);
x = [cos(φ) * sin(θ) for θ in θ, φ in φ];
y = [sin(φ) * sin(θ) for θ in θ, φ in φ];
z = [cos(θ) for θ in θ, φ in φ];

fig = Figure()
ax1 = Axis3(fig[1,1], aspect=:data, viewmode=:fitzoom)
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

ax2 = GeoAxis(
    fig[1, 2]; 
    title="Satellite Groundtrack",
    coastlines = true,
    coastline_attributes = (; color = "white"),
)

# orbital_elements, r_apo, r_per, t_init, t_final;
# a, e, i, asc, arg, theta = orbital_elements

elements = Observable([
    26600,
    0.74,
    deg2rad(63.4),
    deg2rad(90),
    deg2rad(270),
    deg2rad(230)
]);

sl = Slider(fig[2, 2], range = 50:1000, startvalue = 100, horizontal = true)

const eq_rad = 6378.137 
const ecc_sq = 6.69437999014e-3

clamp_lon(lon) = lon % 180 - sign(lon)*180

# initial plot
lon, lat, alt, times = ground_track(elements, 0.0,  1*718*60.0, num_steps=100, return_times=true);
ground_track_coords = Observable(Matrix(hcat(clamp_lon.(rad2deg.(Float64.(lon))), rad2deg.(Float64.(lat)))'));

N = eq_rad ./ sqrt.(1 .- ecc_sq .* sin.(lat).^2);

pos_x = (N .+ alt) .* cos.(lat) .* cos.(lon) ./ eq_rad;
pos_y = (N .+ alt) .* cos.(lat) .* sin.(lon) ./ eq_rad;
pos_z = ((1 - ecc_sq) .* N .+ alt) .* sin.(lat) ./ eq_rad;

pos_coords = Observable(Matrix(hcat(pos_x, pos_y, pos_z)'));

scatter!(ax2, ground_track_coords, marker=:cross, color=:orange)
lines!(ax1, pos_coords, linestyle = :dash, color=:orange)


on(sl.value) do ns
    lon, lat, alt, times = ground_track(elements, 0.0, 2.5*718*60.0, num_steps=ns, return_times=true)
    ground_track_coords[] = Matrix(hcat(clamp_lon.(rad2deg.(Float64.(lon))), rad2deg.(Float64.(lat)))');

    N = eq_rad ./ sqrt.(1 .- ecc_sq .* sin.(lat).^2)

    pos_x = (N .+ alt) .* cos.(lat) .* cos.(lon) ./ eq_rad
    pos_y = (N .+ alt) .* cos.(lat) .* sin.(lon) ./ eq_rad
    pos_z = ((1 - ecc_sq) .* N .+ alt) .* sin.(lat) ./ eq_rad

    pos_coords[] = Matrix(hcat(pos_x, pos_y, pos_z)')

    # Observable(Matrix(hcat(rad2deg.(Float64.(lon)), rad2deg.(Float64.(lat)))')), Observable(Float64.(pos_x)), Observable(Float64.(pos_y)), Observable(Float64.(pos_z))
end



# sg = SliderGrid(fig[1, 1],
#     (label = "Amplitude", range = 0:0.1:10, startvalue = 5),
#     (label = "Frequency", range = 0:0.5:50, format = "{:.1f}Hz", startvalue = 10),
#     (label = "Phase", range = 0:0.01:2pi,
#         format = x -> string(round(x/pi, digits = 2), "π"))
# )

# on(sg.sliders[1].value) do val
#     # do something with `val`
# end

