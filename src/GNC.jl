using LinearAlgebra, Dates

function determine_eccentric_anomaly(eccentricity, mean_anomaly; tolerance=1e-8)
    if mean_anomaly < pi
        eccentric_anomaly = mean_anomaly + eccentricity/2
    else
        eccentric_anomaly = mean_anomaly - eccentricity/2
    end

    function ratio(eccentricty, mean_anomaly, eccentric_anomaly)
        return (eccentric_anomaly - eccentricity*sin(eccentric_anomaly) - mean_anomaly) / (1 - eccentricity*cos(eccentric_anomaly))
    end

    current_ratio = ratio(eccentricity, mean_anomaly, eccentric_anomaly)

    while abs(current_ratio) > tolerance
        eccentric_anomaly = eccentric_anomaly - current_ratio
        current_ratio = ratio(eccentricity, mean_anomaly, eccentric_anomaly)
    end

    return eccentric_anomaly
end

function geocentric_state_vector_transform(angular_momentum, orbital_elements; mu=398600)
    h = angular_momentum
    e, i, asc, arg, theta = orbital_elements

    p⃗_perifocal = ((h^2 / mu) / (1 + e * cos(theta))) * [cos(theta); sin(theta); 0]
    v⃗_perifocal = (mu / h) * [-sin(theta); e + cos(theta); 0]

    transform = [
        -sin(asc) * cos(i) * sin(arg) + cos(asc) * cos(arg) -sin(asc) * cos(i) * cos(arg) - cos(asc) * sin(arg) sin(asc) * sin(i)
        cos(asc) * cos(i) * sin(arg) + sin(asc) * cos(arg) cos(asc) * cos(i) * cos(arg) - sin(asc) * sin(arg) -cos(asc) * sin(i)
        sin(i) * sin(arg) sin(i) * cos(arg) cos(i)
    ]

    pos_vec_geocentric = transform * p⃗_perifocal
    vel_vec_geocentric = transform * v⃗_perifocal

    return pos_vec_geocentric, vel_vec_geocentric
end

function position_to_asc_dec(position_vector)
    magnitude = sqrt(sum(position_vector.^2))

    direct_l = position_vector[1] / magnitude
    direct_m = position_vector[2] / magnitude
    direct_n = position_vector[3] / magnitude

    declination = asin(direct_n)

    direct_m > 0 ? ascension = acos(direct_l / cos(declination)) : ascension = 2*pi - acos(direct_l / cos(declination))

    return ascension, declination
end

function ground_track(orbital_elements, t_init, t_final; num_steps=100, return_times=false, mu=398600, radius=6378, j2=1.08263e-3, body_angular_vel=7.292124e-5)
    a, e, i, asc, arg, theta = orbital_elements
    r_per = a*(1 - e)

    common_term = (3 / 2) * ((sqrt(mu) * j2 * radius^2) / ((1 - e^2) * a^(7/2)))
    delta_right_ascension = -common_term * cos(i)
    delta_argument_perigee = -common_term * ((5 / 2) * sin(i)^2 - 2)

    angular_momentum = sqrt(mu * r_per * (1 + e))
    period = (a^(3/2)) * 2*pi / sqrt(mu)

    eccentric_anomaly_0 = 2 * atan(tan(theta / 2) * sqrt((1 - e) / (1 + e)))
    mean_anomaly_0 = eccentric_anomaly_0 - e * sin(eccentric_anomaly_0)
    t_0 = period * mean_anomaly_0 / (2 * pi)

    longitude_vec = []
    latitude_vec = []
    altitude_vec = []
    times = []

    earth_rot = 2 * pi / (24 * 60 * 60)

    step_size = (t_final - t_init) / num_steps
    for c=1:num_steps
        t_i = t_0 + c * step_size

        mean_anomaly_i = t_i * (2 * pi) / period
        eccentric_anomaly_i = determine_eccentric_anomaly(e, mean_anomaly_i)
        true_anomaly_i = 2 * atan(tan(eccentric_anomaly_i / 2) * sqrt((1 + e) / (1 - e)))

        right_ascension_i = asc + step_size * delta_right_ascension
        argument_perigee_i = arg + step_size * delta_argument_perigee

        pos_geo, vel_geo = geocentric_state_vector_transform(angular_momentum, [e, i, right_ascension_i, argument_perigee_i, true_anomaly_i])

        transform_theta = step_size * body_angular_vel
        rotation_matrix = [
            cos(transform_theta) sin(transform_theta) 0
            -sin(transform_theta) cos(transform_theta) 0
            0 0 1
        ]

        pos_geo_rotating = rotation_matrix * pos_geo[:,:]
        push!(altitude_vec, norm(pos_geo_rotating))

        longitude, latitude = position_to_asc_dec(pos_geo_rotating)

        push!(longitude_vec, longitude - (c * step_size) * earth_rot)
        push!(latitude_vec, latitude)
        push!(times, t_i)

    end

    if return_times 
        return longitude_vec, latitude_vec, altitude_vec, times 
    else 
        return longitude_vec, latitude_vec, altitude_vec
    end
end

function get_orbital_elements(r⃗, v⃗; μ=398600, elements_type="standard")
    r = norm(r⃗)
    v = norm(v⃗)

    h⃗ = r⃗ × v⃗
    h = norm(h⃗)

    i = acos(h⃗[3] / h)

    N⃗ = [0, 0, 1] × h⃗
    N = norm(N⃗)

    N⃗[2] >= 0 ? Ω = acos(N⃗[1] / N) : Ω = 2 * pi - acos(N⃗[1] / N)

    v_rad = (r⃗ ⋅ v⃗) / r
    e⃗ = (1 / μ) * ((v^2 - μ / r) * r⃗ - r * v_rad * v⃗)
    e = norm(e⃗)

    e⃗[3] >= 0 ? ω = acos(dot(N⃗ / N, e⃗ / e)) : ω = 2 * π - acos(dot(N⃗ / N, e⃗ / e))

    v_rad >= 0 ? θ = acos(dot(e⃗ / e, r⃗ / r)) : θ = 2 * π - acos(dot(e⃗ / e, r⃗ / r))

    a = (r * (1 + e * cos(θ))) / (1 - e^2)

    if elements_type == "standard"
        return [a, e, i, Ω, ω, θ]
    elseif elements_type == "curtis"
        return [h, i, Ω, e, ω, θ]
    end
end

function gibbs_method(r⃗₁, r⃗₂, r⃗₃; μ=398600, error_tolerance=1e-5, return_state_vector=false, elements_type="standard")
	r₁, r₂, r₃ = norm(r⃗₁), norm(r⃗₂), norm(r⃗₃)
	C⃗₁₂, C⃗₂₃, C⃗₃₁ = r⃗₁ × r⃗₂, r⃗₂ × r⃗₃, r⃗₃ × r⃗₁

    coplanar_error = abs((r⃗₁ ./ r₁) ⋅ (C⃗₂₃ ./ norm(C⃗₂₃)))
	if coplanar_error > error_tolerance
		error("Error: Vectors are not coplanar, |ûᵣ₁ ⋅ Ĉ₂₃| = $(round(coplanar_error, digits=8)) > $(error_tolerance)\nAdjust error_tolerance or choose different position vectors")
	
	else
		N⃗ = r₁*C⃗₂₃ + r₂*C⃗₃₁ + r₃*C⃗₁₂
		D⃗ = C⃗₁₂ + C⃗₂₃ + C⃗₃₁
		S⃗ = r⃗₁.*(r₂ - r₃) + r⃗₂.*(r₃ - r₁) + r⃗₃.*(r₁ - r₂)

		v⃗₂ = sqrt(μ / (norm(N⃗) * norm(D⃗))) * ((D⃗ × r⃗₂)./r₂ + S⃗)

        if return_state_vector
		    return [r⃗₂, v⃗₂], get_orbital_elements(r⃗₂, v⃗₂, μ=μ, elements_type=elements_type)
        else
            return get_orbital_elements(r⃗₂, v⃗₂, μ=μ, elements_type=elements_type)
        end
	end
end