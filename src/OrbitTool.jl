include("../src/SAT.jl");

using GLMakie, GeoMakie, FileIO, .SAT
using Downloads: download
GLMakie.activate!()

Makie.set_theme!(theme_dark())

# source: https://beautiful.makie.org/dev/examples/generated/2d/geo/blue_marble/
earth_img = load(download("https://upload.wikimedia.org/wikipedia/commons/5/56/Blue_Marble_Next_Generation_%2B_topography_%2B_bathymetry.jpg"));
n = 1024 ÷ 4 # 2048
θ = LinRange(0, π, n);
φ = LinRange(0, 2π, 2 * n);
x = [cos(φ) * sin(θ) for θ in θ, φ in φ];
y = [sin(φ) * sin(θ) for θ in θ, φ in φ];
z = [cos(θ) for θ in θ, φ in φ];

begin
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

# tb_num_steps = Textbox(ax3[1,2], stored_string="100", validator = Int, textcolor=:white)
num_steps_obsv = Observable(250)
num_steps_sl = Slider(ax3[1,2], range = 50:750, startvalue = to_value(num_steps_obsv), horizontal = true)
connect!(num_steps_obsv, num_steps_sl.value)
Label(ax3[1,1], @lift("Number of Computed Points: $($num_steps_obsv)"))

final_time_obsv = Observable(24) # in hours
final_time_sl = Slider(ax3[2,2], range = 1:50, startvalue = to_value(final_time_obsv), horizontal = true)
connect!(final_time_obsv, final_time_sl.value)
Label(ax3[2,1], @lift("Orbit Time: $($final_time_obsv) hr"))

# projections = ["+proj=adams_hemi", "+proj=adams_ws1", "+proj=adams_ws2",
#     "+proj=aea +lat_1=29.5 +lat_2=42.5", "+proj=aeqd", "+proj=airy", "+proj=aitoff",
#     "+proj=apian", "+proj=august", "+proj=bacon", "+proj=bertin1953", "+proj=bipc +ns",
#     "+proj=boggs", "+proj=bonne +lat_1=10", "+proj=cass", "+proj=cea",
#     "+proj=chamb +lat_1=10 +lon_1=30 +lon_2=40", "+proj=collg", "+proj=comill",
#     "+proj=crast", "+proj=denoy", "+proj=eck1", "+proj=eck2", "+proj=eck3",
#     "+proj=eck4", "+proj=eck5", "+proj=eck6", "+proj=eqc", "+proj=eqdc +lat_1=55 +lat_2=60",
#     "+proj=eqearth", "+proj=euler +lat_1=67 +lat_2=75", "+proj=fahey", "+proj=fouc", "+proj=fouc_s",
#     "+proj=gall", "+proj=geos +h=35785831.0 +lon_0=-60 +sweep=y", "+proj=gins8", "+proj=gn_sinu +m=2 +n=3",
#     "+proj=goode", "+proj=guyou", "+proj=hammer", "+proj=hatano",
#     "+proj=igh", "+proj=igh_o +lon_0=-160", "+proj=imw_p +lat_1=30 +lat_2=-40", "+proj=isea",
#     "+proj=kav5", "+proj=kav7", "+proj=laea", "+proj=lagrng", "+proj=larr", "+proj=lask",
#     "+proj=lcca +lat_0=35", "+proj=leac", "+proj=loxim",
#     "+proj=lsat +ellps=GRS80 +lat_1=-60 +lat_2=60 +lsat=2 +path=2", "+proj=mbt_s", "+proj=mbt_fps",
#     "+proj=mbtfpp", "+proj=mbtfpq", "+proj=mbtfps", "+proj=merc", "+proj=mill", "+proj=misrsom +path=1",
#     "+proj=moll", "+proj=murd1 +lat_1=30 +lat_2=50",
#     "+proj=murd3 +lat_1=30 +lat_2=50", "+proj=natearth", "+proj=natearth2",
#     "+proj=nell", "+proj=nell_h", "+proj=nicol",
#     "+proj=ob_tran +o_proj=mill +o_lon_p=40 +o_lat_p=50 +lon_0=60", "+proj=ocea", "+proj=oea +m=1 +n=2",
#     "+proj=omerc +lat_1=45 +lat_2=55", "+proj=ortel", "+proj=ortho", "+proj=patterson", "+proj=poly",
#     "+proj=putp1", "+proj=putp2", "+proj=putp3", "+proj=putp3p", "+proj=putp4p", "+proj=putp5",
#     "+proj=putp5p", "+proj=putp6", "+proj=putp6p", "+proj=qua_aut", "+proj=robin", "+proj=rouss",
#     "+proj=rpoly", "+proj=sinu", "+proj=times", "+proj=tissot +lat_1=60 +lat_2=65", "+proj=tmerc",
#     "+proj=tobmerc", "+proj=tpeqd +lat_1=60 +lat_2=65", "+proj=urm5 +n=0.9 +alpha=2 +q=4",
#     "+proj=urmfps +n=0.5", "+proj=vandg", "+proj=vandg2", "+proj=vandg3", "+proj=vandg4",
#     "+proj=vitk1 +lat_1=45 +lat_2=55", "+proj=wag1", "+proj=wag2", "+proj=wag3", "+proj=wag4",
#     "+proj=wag5", "+proj=wag6", "+proj=wag7", "+proj=webmerc +datum=WGS84", "+proj=weren",
#     "+proj=wink1", "+proj=wink2", "+proj=wintri"
# ]

# projections = ["+proj=weren", "+proj=crast",
#     "+proj=wink1", "+proj=wink2", "+proj=wintri", "+proj=poly",
#     "+proj=putp1", "+proj=putp2", "+proj=putp3", "+proj=putp3p", "+proj=putp4p", "+proj=putp5",
#     "+proj=putp5p", "+proj=putp6", "+proj=putp6p", "+proj=qua_aut", "+proj=robin", "+proj=rouss",
#     "+proj=rpoly", "+proj=sinu"
# ]

projections = [
    "+proj=natearth", "+proj=natearth2", "+proj=moll", "+proj=weren", "+proj=times", "+proj=wink1", "+proj=wink2","+proj=sinu", "+proj=crast",
]

proj_obsv = Observable("")
Label(ax3[1,3], "Map Projection")
proj_menu = Menu(ax3[2,3], options=projections)
connect!(proj_obsv, proj_menu.selection)

on(proj_menu.selection) do s
    proj_obsv[] = s
    # autolimits!(ax2)
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

# img = rotr90(GeoMakie.earth());
# image!(ax2, -180..180, -90..90, img; interpolate = true) # this must be included
# surface!(ax2,  -180..180, -90..90, ones(size(img)...); color = img, shading = false)

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

# initial plot
lon, lat, alt, times = ground_track(to_value(elements), 0.0, to_value(final_time_obsv)*60^2, num_steps=to_value(num_steps_obsv), return_times=true);

# println("Lat: $(rad2deg(maximum(lat))) $(rad2deg(minimum(lat)))")
# println("Lat Clamp: $(clamp_lat(rad2deg(maximum(lat)))) $(clamp_lat(rad2deg(minimum(lat))))")
# println("Lon: $(rad2deg(maximum(lon))) $(rad2deg(minimum(lon)))")
# println("Lon Clamp: $(clamp_lon(rad2deg(maximum(lon)))) $(clamp_lon(rad2deg(minimum(lon))))")

ground_track_coords = Observable(Matrix(hcat(clamp_lon.(rad2deg.(Float64.(lon))), clamp_lat.(rad2deg.(Float64.(lat))))'));

const eq_rad = 6378.137 
const ecc_sq = 6.69437999014e-3
N = eq_rad ./ sqrt.(1 .- ecc_sq .* sin.(lat).^2);

pos_x = (N .+ alt) .* cos.(lat) .* cos.(lon) ./ eq_rad;
pos_y = (N .+ alt) .* cos.(lat) .* sin.(lon) ./ eq_rad;
pos_z = ((1 - ecc_sq) .* N .+ alt) .* sin.(lat) ./ eq_rad;

pos_coords = Observable(Matrix(hcat(pos_x, pos_y, pos_z)'));

scatter!(ax2, ground_track_coords, marker=:cross, color=:orange)
lines!(ax1, pos_coords, linestyle = :dash, color=:orange)

on(num_steps_sl.value) do ns
    lon, lat, alt, times = ground_track(to_value(elements), 0.0, to_value(final_time_obsv)*60^2, num_steps=ns, return_times=true)
    ground_track_coords[] = Matrix(hcat(clamp_lon.(rad2deg.(Float64.(lon))), clamp_lat.(rad2deg.(Float64.(lat))))')

    N = eq_rad ./ sqrt.(1 .- ecc_sq .* sin.(lat).^2)

    pos_x = (N .+ alt) .* cos.(lat) .* cos.(lon) ./ eq_rad
    pos_y = (N .+ alt) .* cos.(lat) .* sin.(lon) ./ eq_rad
    pos_z = ((1 - ecc_sq) .* N .+ alt) .* sin.(lat) ./ eq_rad

    pos_coords[] = Matrix(hcat(pos_x, pos_y, pos_z)')
end

on(final_time_sl.value) do tf
    lon, lat, alt, times = ground_track(to_value(elements), 0.0, tf*60^2, num_steps=to_value(num_steps_obsv), return_times=true)
    ground_track_coords[] = Matrix(hcat(clamp_lon.(rad2deg.(Float64.(lon))), clamp_lat.(rad2deg.(Float64.(lat))))')

    N = eq_rad ./ sqrt.(1 .- ecc_sq .* sin.(lat).^2)

    pos_x = (N .+ alt) .* cos.(lat) .* cos.(lon) ./ eq_rad
    pos_y = (N .+ alt) .* cos.(lat) .* sin.(lon) ./ eq_rad
    pos_z = ((1 - ecc_sq) .* N .+ alt) .* sin.(lat) ./ eq_rad

    pos_coords[] = Matrix(hcat(pos_x, pos_y, pos_z)')
end

end

elements_sg = SliderGrid(
    ax3[1,4],
    (label = "Semi-Major Axis", range = 7000:45000, format = "{:.2f} km", startvalue = 8350),
    (label = "Eccentricity", range = 0:0.001:0.999, format = "{:.3f}", startvalue = 0.19760),
    (label = "Inclination", range = 0:180, format = "{:.2}°", startvalue = 60),
    (label = "RAAN", range = 0:360, format = "{:.2}°", startvalue = 270),
    (label = "Argument of Periapsis", range = 0:360, format = "{:.2}°", startvalue = 45),
    (label = "True Anomaly", range = 0:360, format = "{:.2}°", startvalue = 230),
    width = 750,
)

colsize!(ax3, 1, Relative(1/3))
colsize!(ax3, 2, Relative(1/3))

sliderobservables = [s.value for s in elements_sg.sliders]
bars = lift(sliderobservables...) do slvalues...
    a, e, i, asc, arg, theta  = [slvalues...]
    # a, e, i = [slvalues...]
    old_els = to_value(elements)

    old_els[1] = a
    old_els[2] = e
    old_els[3] = deg2rad(i)
    old_els[4] = deg2rad(asc)
    old_els[5] = deg2rad(arg)
    old_els[6] = deg2rad(theta)

    elements[] = old_els

    lon, lat, alt, times = ground_track(to_value(elements), 0.0, to_value(final_time_obsv)*60^2, num_steps=to_value(num_steps_obsv), return_times=true)
    ground_track_coords[] = Matrix(hcat(clamp_lon.(rad2deg.(Float64.(lon))), clamp_lat.(rad2deg.(Float64.(lat))))')

    N = eq_rad ./ sqrt.(1 .- ecc_sq .* sin.(lat).^2)

    pos_x = (N .+ alt) .* cos.(lat) .* cos.(lon) ./ eq_rad
    pos_y = (N .+ alt) .* cos.(lat) .* sin.(lon) ./ eq_rad
    pos_z = ((1 - ecc_sq) .* N .+ alt) .* sin.(lat) ./ eq_rad

    pos_coords[] = Matrix(hcat(pos_x, pos_y, pos_z)')

    [slvalues...]
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

# sg = SliderGrid(
#     fig[1, 2],
#     (label = "Voltage", range = 0:0.1:10, format = "{:.1f}V", startvalue = 5.3),
#     (label = "Current", range = 0:0.1:20, format = "{:.1f}A", startvalue = 10.2),
#     (label = "Resistance", range = 0:0.1:30, format = "{:.1f}Ω", startvalue = 15.9),
#     width = 350,
#     tellheight = false)

# sliderobservables = [s.value for s in sg.sliders]
# bars = lift(sliderobservables...) do slvalues...
#     [slvalues...]
# end