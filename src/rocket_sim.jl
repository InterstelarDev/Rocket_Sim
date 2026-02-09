using Plots, Statistics, Random

#Fuels(Fuel Name, Specific, Dry mass, Fuel mass, Graph color, Diameter of rocket)
const fuels = [
    (name = "RP-1/LOX", Isp = 340, dry_mass = 120, fuel_mass = 1200, color = :blue, Diameter = 0.8),
    (name = "LCH4/LOX (Methane)", Isp = 375, dry_mass = 150, fuel_mass = 1400, color = :green, Diameter = 1.0),
    (name = "E-methane", Isp = 360, dry_mass = 150, fuel_mass = 1400, color = :orange, Diameter = 1.0),
    (name = "LH2/LOX (Hydrogen)", Isp = 452, dry_mass = 300, fuel_mass = 2500, color = :red, Diameter = 2.2),
    (name = "Green Hydrogen", Isp = 450, dry_mass = 300, fuel_mass = 2500, color = :black, Diameter = 2.2),
    (name = "UDMH/NTO", Isp = 315, dry_mass = 110, fuel_mass = 1000, color = :brown, Diameter = 0.7),
    (name = "Green Ammonia", Isp = 310, dry_mass = 135, fuel_mass = 950, color = :purple, Diameter = 0.85),
    (name = "Solid (HTPB)", Isp = 285, dry_mass = 100, fuel_mass = 2000, color = :gray, Diameter = 0.6)
]

# Calculates atmospheric density using altitude
# if altitude > 100000 there is no air, so no density to calculate
# if altitude < 11000, we are in the troposphere, where density decreases exponentially with altitude by around 8420 m
# if altitude < 20000, we are in the lower stratosphere, where density decreases exponentially with altitude by around 6350 m
# if altitude > 20000, we are in the upper stratosphere, where density decreases exponentially with altitude by around 6800 m 
function calc_density(altitude)
    if altitude > 100000
        return 0.0
    end
    if altitude < 11000
        return 1.225 * exp(-altitude / 8420)
    elseif altitude < 20000
        return 0.3639 * exp(-(altitude - 11000) / 6350)
    else
        return 0.0880 * exp(-(altitude - 20000) / 6800)
    end
end

# Calculates the drag coefficient using velocity and altitude
# Simple ternary operator, if altitude > 11000, speed of sound is 295 m/s, otherwise it is 340 m/s
# then the mach number is calculated by dividing the absolute value of velocity by the speed of sound
# Mach is the ratio of the speed of an object to the speed of sound of the medium it is traveling through.
# if mach < 0.8, the drag coefficient is 0.2
# if mach < 1.2, the drag coefficient is 0.2 + (mach - 0.8), which increases linearly from 0.2 at mach 0.8 to 0.6 at mach 1.2
# this is done to simulate the increase in drag as the rocket approaches the speed of sound.
# this is known as the transonic regime, where drag increases significantly due to shock waves forming around the rocket.
# if mach >= 1.2, the drag coefficient is 0.6 / mach, starting at 0.6 at mach 1.2 and decreasing as mach increases 
# this is done to simulate the decrease in drag as the rocket goes supersonic, as shock waves become more stable and the flow becomes more streamlined.
function calc_drag_coefficient(velocity, altitude)
    speed_of_sound = (altitude > 11000) ? 295.0 : 340.0
    mach = abs(velocity) / speed_of_sound
    if mach < 0.8
        return 0.2       
    elseif mach < 1.2
        return 0.2 + (mach - 0.8)  
    else
        return 0.6 / mach 
    end
end

# Simulates the rocket's ascent and returns time, altitude, maximum altitude, and maximum dynamic pressure
function rocket_simulation(Isp, dry_mass, fuel_mass, Diameter, burn_time=50)
    # defines gravity constant and time interval
    g0, ti = 9.81, 0.1
    # defines initial time steps for the burn phase as time_steps, which is an Int of burn_time/ti
    time_steps = Int(burn_time / ti)
    # defines time steps for the coast phase as coast_steps, which is an Int of 1000/ti, or 10000 steps
    coast_steps = Int(1000 / ti)
    # defines total_steps, which is the sum of time_steps and coast_steps
    total_steps = time_steps + coast_steps

    # creates a array from 0 to total_steps * ti with the length of total_steps, defined as time
    time = range(0, total_steps * ti, length=total_steps)
    # creates arrays with zeroes for velocity and altitude
    velocity, altitude = zeros(total_steps), zeros(total_steps)
    # defines max_q, which will be used to track the maximum dynamic pressure experienced by the rocket during ascent
    max_q = 0.0
    # defines a q_limit of 100000 Pa
    q_limit = 100000.0
    # defines a safety margin of 5000 Pa, which will be used to throttle the rockt if pressure is greater than this
    q_safety_margin = 5000.0

    # calculates the mdot(mass flow rate) using fuel_mass/burn_time
    mdot = fuel_mass / burn_time
    # calculates thrust using mdot * Isp * g0
    full_thrust = mdot * Isp * g0
    # calculates the cross-sectional area of the rocket
    Area = π * (Diameter / 2)^2
    # sets the payload diameter to 0.5 m, which will be used to calculate the drag area after burnout
    # this simulates the change in shape of the rocket as it sheds its fairing and exposes the payload
    payload_diameter = 0.5 

    # for loop, for every step:
    # calculates the atmospheric density at the current altitude using calc_density()
    # sets the velocity to the previous velocity, which will be updated in the loop
    # calculates the drag coefficient using calc_drag_coefficient()
    # calculates the dynamic pressure(q) using q = 0.5 * rho * v^2
    # updates max_q if the current q is greater than max_q
    # calculates the drag force using drag_force = q * Cd * Area
    # if max_q is larger than q_limit, the rocket breaks apart and the simulation ends
    # if the dynamic pressure is greater than the safety margin and we are still in the burn phase, it throttles the thrust to 60% of full thrust
    # else, it uses full thrust
    # if the altitude is greater than 60000 m, it assumes the rocket has shed its fairing and updates the drag area to be the payload diameter
    # if the current step is within the burn time, it calculates the mass of the rocket by using the fuel consumed - initial mass, then calculates acceleration
    # if the velocity is negative and altitude is less than 5000 m, it simulates a parachute deploying by setting a high drag coefficient and area
    # updates the velocity and altitude using the calculated acceleration and time interval
    # if the altitude drops below 0, it sets the altitude and velocity to 0, simulating the rocket hitting the ground, and breaks the loop if this is after burn time
    for i in 2:total_steps
        rho = calc_density(altitude[i-1])
        v = velocity[i-1]
        
        Cd = calc_drag_coefficient(v, altitude[i-1])
        q = 0.5 * rho * v^2
        max_q = max(max_q, q)

        if q > q_limit
            return time[1:i-1], altitude[1:i-1], velocity[1:i-1], altitude[i-1], max_q
        end

        current_thrust = (q > q_safety_margin && i <= time_steps) ? (full_thrust * 0.6) : full_thrust

        drag_force = q * Cd * Area

        if altitude[i-1] > 60000
            Area = π * (payload_diameter / 2)^2 
        end

        if i <= time_steps
            m = (dry_mass + fuel_mass) - (mdot * (i-1) * ti)
            acc = (current_thrust - drag_force) / m - g0
        else
            acc = -g0 - (drag_force / dry_mass)
            if velocity[i-1] < 0 && altitude[i-1] < 5000
                Cd = 1.5 
                Area = 20.0 
            end
        end

        velocity[i] = velocity[i-1] + acc * ti
        altitude[i] = altitude[i-1] + velocity[i] * ti

        if altitude[i] < 0
            altitude[i] = 0; velocity[i] = 0
            if i > time_steps break end
        end
    end
    
    # returns the time, altitude, velocity, maximum altitude, and maximum dynamic pressure
    return time, altitude, velocity, maximum(altitude), max_q
end

# Altitude Graph:
# Plots a Graph of Altitude vs Time for all fuels
function altitude_graph(fuels, burn_time=50)
    plt = plot(title="Altitude vs Time for Different Rocket Fuels", 
                xlabel="Time (s)", ylabel="Altitude (km)", 
                legend=:topright, size=(1000, 600))
    for fuel in fuels
        time, altitude, _, _, _ = rocket_simulation(fuel.Isp, fuel.dry_mass, fuel.fuel_mass, fuel.Diameter, burn_time)
        plot!(plt, time, altitude/1000, label="$(fuel.name) (Isp: $(fuel.Isp))", color=fuel.color)
    end
    display(plt)
end

# Sensitivity Analysis Graph:
# Plots a Graph of Altitude vs Time for all fuels, but with the Isp now varied by ±10% 
# This shows how sensitive the fuel is to changes in Isp
function sensitivity_graph(fuels, burn_time=50)
    plt = plot(title="Sensitivity Analysis: Altitude with Isp ±10%", 
               xlabel="Time (s)", ylabel="Altitude (km)", 
               legend=:topright, size=(1000, 600))
    
    factors = [0.9, 1.0, 1.1]
    styles = [:dash, :solid, :dot]

    for fuel in fuels
        for (factor, style) in zip(factors, styles)
            isp_varied = fuel.Isp * factor
            time, altitude, _, _, _ = rocket_simulation(isp_varied, fuel.dry_mass, fuel.fuel_mass, fuel.Diameter, burn_time)
            
            plot!(plt, time, altitude ./ 1000, 
                  label="$(fuel.name) Isp*$(round(factor, digits=1))", 
                  color=fuel.color, linestyle=style, lw=1.5)
        end
    end
    display(plt)
end

# Monte Carlo Simulation:
# Runs a Monte Carlo simulation for each fuel, where the Isp and fuel mass are randomly varied by ±5% 
# This simulates realistic variability in fuel performance
function monte_carlo_simulation(fuels, burn_time=50, runs=100000)
    plt = plot(title="Monte Carlo Apogee Distribution ($runs runs)", 
               xlabel="Apogee Altitude (km)", ylabel="Frequency", 
               size=(1000, 600))
    
    for fuel in fuels
        apogees = Float64[]
        for _ in 1:runs
            isp_mc = fuel.Isp + (fuel.Isp * 0.05) * randn()
            fm_mc = fuel.fuel_mass + (fuel.fuel_mass * 0.05) * randn()
            fm_mc = clamp(fm_mc, 1.0, 2 * fuel.fuel_mass)
            
            _, _, _, apogee, _ = rocket_simulation(isp_mc, fuel.dry_mass, fm_mc, fuel.Diameter, burn_time)
            push!(apogees, apogee / 1000)
        end
        
        histogram!(plt, apogees, bins=30, alpha=0.5, label=fuel.name, color=fuel.color)
    end
    display(plt)
end

# Velocity Graph:
# Plots a Graph of Velocity vs Time for all fuels
function velocity_graph(fuels, burn_time=50)
    plt = plot(title="Velocity Profiles for All Fuels", 
               xlabel="Time (s)", ylabel="Velocity (m/s)", 
               legend=:topright, size=(1000, 600), lw=2)
    
    for fuel in fuels
        time, _, velocity, _, _ = rocket_simulation(fuel.Isp, fuel.dry_mass, fuel.fuel_mass, fuel.Diameter, burn_time)
        plot!(plt, time, velocity, label=fuel.name, color=fuel.color)
    end
    display(plt)
end

# Acceleration Graph:
# Plots a Graph of Acceleration vs Time for all fuels
function acceleration_graph(fuels, burn_time=50)
    plt = plot(title="Real-Time Acceleration (With Drag & Throttle)", 
               xlabel="Time (s)", ylabel="Acceleration (m/s²)", 
               size=(1200, 800), lw=2)
    
    for fuel in fuels
        time, _, velocity, _, _ = rocket_simulation(fuel.Isp, fuel.dry_mass, fuel.fuel_mass, fuel.Diameter, burn_time)
        accel = diff(velocity) ./ 0.1
        
        plot!(plt, time[1:end-1], accel, label=fuel.name, color=fuel.color)
    end
    display(plt)
end

# Runs all the graphs
altitude_graph(fuels)
sensitivity_graph(fuels)
monte_carlo_simulation(fuels)
velocity_graph(fuels)
acceleration_graph(fuels)
