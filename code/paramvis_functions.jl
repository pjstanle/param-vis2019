

"""
Reads an optimization output file and returns an array `M` where `M[i, j]`
is the value associated to the i-th turbine in the j-th iteration.
"""
function read_opt_output(filename::String, path::String, nturbines::Int)

    # Parse file into values
    vals = Float64[]
    open(joinpath(path, filename)) do file
        for line in eachline(file)
            for elem in split(line, ' ')
                if length(elem) != 0
                    push!(vals, parse(elem))
                end
            end
        end
    end

    # Check that the number of values is a multiple of the number of turbines
    if length(vals)%nturbines != 0
        error("Number of entries ($(length(valus))) is not divisable by the"*
                " number of turbines ($nturbines).")
    end

    # Reshape into a matrix
    return reshape(vals, (nturbines, Int(length(vals)/nturbines)))
end

"""
Generate a Gaussian wake function to feed `generate_windfarm(...)`.
"""
function generate_wake_gaussian(magVinf::Real, wind_direction::Real,
                                turbine_x, turbine_y, turbine_z,
                                hub_height, rotor_diameter;
                                    alpha::Real=0.15, refH::Real=-1,
                                    shear::Bool=true,
                                    k::Real=0.0325, CT::Real=8.0/9.0)


    # Guassian wake function
    # magVinf = 12.0                        # (m/s) wind speed
    Vinf_Oaxis = gt.rotation_matrix2(0, 0, wind_direction) # Free-stream coordinate system
    invVinf_Oaxis = Vinf_Oaxis'             # Inverse transformation
    # refH = mean(hub_height)               # (m) reference height for shear layer

    if shear && refH==-1
        error("Shear layer requested without a reference height (refH).")
    end

    function wake_fun(X)

        loss = 0

        # Calculate loss in freestream direction
        for i = 1:length(turbine_x)     # Iterate over turbines

            # Distance from this turbine in global coordinates
            DX = X - [turbine_x[i], turbine_y[i], turbine_z[i]+hub_height[i]]

            # Distance in free-stream coordinates
            dX = invVinf_Oaxis*DX
            dx = dX[1]

            if dx > 0.
                dy = dX[2]
                dz = dX[3]
                sigma = k*dx + rotor_diameter[i]/sqrt(8.0)
                this_loss = (1.0-sqrt(1.0-CT/(8*sigma^2/(rotor_diameter[i]^2))))*
                               exp(-1.0/2.0*(dy/sigma)^2)*exp(-1.0/2.0*(dz/sigma)^2)
                loss += this_loss^2
            end

        end
        loss = sqrt(loss)

        windspeed= magVinf

        if shear
            windspeed *= (X[3]/refH)^alpha
        end

        # Return velocity vector in global coordinate system
        return Vinf_Oaxis*[ windspeed*(1-loss), 0, 0 ]
    end

    return wake_fun
end
