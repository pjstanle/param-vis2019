
# ------------------------------------------------------------------------------
function generate_vtks(; data_path="data/direct4",
                            save_path="temps/direct4_00",
                            filepref="windfarm",
                            prompt=true, verbose=true, v_lvl=0,
                            previs=false, paraview=false)

    # data_path = "data/direct4"              # Optimization data path
    # save_path = "temps/direct4_00"          # Path where to save vtks
    previs = previs && paraview

    # Create save path
    gt.create_path(save_path, prompt)

    # Visualization options
    shear = false                           # Add shear layer


    # --------------------- OPTIMIZATION DATA ----------------------------------

    # Optimization constant parameters
    nturbines = 100                         # Number of turbines in optimization
    hhub = 110.0                            # (m) hub height
    D = 130.0                               # (m) rotor diameter
    nB = 3                                  # Number of blades
    magVinf = 10.0                          # (m/s) wind speed
    wind_direction = 0.0                    # (deg) wind direction
    refH = hhub                             # (m) reference height for shear

    turbine_z = zeros(nturbines)            # (m) z-position of tower base
    hub_height = hhub*ones(nturbines)       # (m) hub height
    rotor_diameter = D*ones(nturbines)      # (m) rotor diameter
    yaw = zeros(nturbines) - wind_direction # (deg) yaw
    nBlades = nB*ones(Int, nturbines)       # Number of blades

    # Optimizationp path
    if verbose; println("\t"^(v_lvl)*"Reading optimization outputs..."); end;
    turbine_xs = read_opt_output("turbineX_history.txt", data_path, nturbines)
    turbine_ys = read_opt_output("turbineY_history.txt", data_path, nturbines)

    nites = size(turbine_xs, 2)             # Number of optimization iterations

    # Generate starting layout
    if verbose; println("\t"^(v_lvl)*"Generating initial layout..."); end;
    windfarm = wfv.generate_layout(rotor_diameter, hub_height, nBlades,
                              turbine_xs[:, 1], turbine_ys[:, 2], turbine_z,
                              yaw;
                              save_path= previs ? save_path : nothing,
                              file_name=filepref, paraview=false
                             )



    # --------------------- PERIMETER AND FLUID DOMAIN -------------------------
    if verbose; println("\t"^(v_lvl)*"Generating farm perimeter..."); end;
    NDIVSx = 70              # Cells in the parametric x-direction
    NDIVSy = 70              # Cells in the parametric y-direction
    NDIVSz = 70              # Cells in the geometric z-direction
    spl_s = 0.1            # Spline smoothness
    spl_k = 1                # Spline degree

    # Read perimeter
    aux = CSV.read(joinpath(data_path, "perimeter.txt"); delim=',',
                                                header=["x", "y"], datarow=1)

    # Physical perimeter of the wind farm
    physical_perimeter = [[aux[i, 1], aux[i, 2], 0.0] for i in 1:size(aux, 1)]

    # Perimeter of fluid domain
    Dadd = 5.0                                    # Number of diameters to add
    C = mean(physical_perimeter)                # Windfarm centroid
    # fdom_perimeter = [C + (X-C)*(1+D*Dadd/norm(X-C)) for X in physical_perimeter]
    zmin = 0.0                                  # Minimum z to mesh
    zmax = hhub + 1.25*D/2                      # Maximum z to mesh

    # Here I stretch the fluid domain perimeter
    fdom_perimeter = [zeros(3) for i in physical_perimeter]
    np = length(physical_perimeter)
    for i in 1:np
        Xi = physical_perimeter[i]
        Xip1 = physical_perimeter[i != np ? i+1 : 2]
        Xim1 = physical_perimeter[i != 1 ? i-1 : np-1]

        # Stretch in the normal direction
        nip1 = cross(Xip1-Xi, [0,0,1])
        nim1 = cross(Xi-Xim1, [0,0,1])
        ni = (nip1 + nim1)/2
        ni /= norm(ni)

        # Extra stretching to make it injective
        extraX = [sign(Xi[1])*0.01*exp(abs(Xi[1]-C[1])/240), 0, 0]

        fdom_perimeter[i][:] = Xi + D*Dadd*ni + extraX
    end


    perimeter = wfv.generate_perimetergrid(physical_perimeter, 25, 25, 0;
                                            z_min=0.0, z_max=0.0,
                                            verify_spline=!true, spl_s=spl_s,
                                            spl_k=spl_k, save_path=save_path,
                                            paraview=false,
                                            file_name=filepref*"_perimeter")

    if verbose; println("\t"^(v_lvl)*"Generating fluid domain..."); end;
    fdom = wfv.generate_perimetergrid(fdom_perimeter,
                                        NDIVSx, NDIVSy, NDIVSz;
                                        z_min=zmin, z_max=zmax,
                                        verify_spline=!true,
                                        spl_s=spl_s, spl_k=spl_k,
                                        save_path=previs ? save_path : nothing,
                                        file_name=filepref*"_fdom",
                                        paraview=false
                                      )

    # --------------------- PREVISUALIZATION -----------------------------------
    # Previsualize farm
    if previs
        if verbose
            println("\t"^(v_lvl)*"Calling Paraview for previsualization...")
        end

        str = save_path*"/"
        for ti in 1:nturbines
            for bi in 1:nB
                str *= filepref*"_turbine$(ti)_rotor_blade$(bi).vtk;"
            end
            str *= filepref*"_turbine$(ti)_rotor_hub.vtk;"
            str *= filepref*"_turbine$(ti)_tower.vtk;"
        end
        str *= filepref*"_perimeter.vtk;"
        str *= filepref*"_fdom.vtk;"

        run(`paraview --data=$str`)
    end

    # --------------------- GENERATE VISUALIZATION DATA ------------------------
    if verbose; println("\t"^(v_lvl)*"Generating visualization data..."); end;

    for ite in 1:nites # Iterate over optimization steps
        if verbose # && (ite-1)%ceil(Int, nites/20)==0
            println("\t"^(v_lvl+1)*"Step $ite out of $nites")
        end

        if ite != 1
            # Translates every turbine
            for ti in 1:nturbines
                # Translation
                dX = [turbine_xs[ti, ite] - turbine_xs[ti, ite-1],
                      turbine_ys[ti, ite] - turbine_ys[ti, ite-1],
                      0]
                gt.lintransform!(gt.get_grid(windfarm, ti), eye(3), dX)
            end
        end

        # Generate wake function based in new positions
        wake_fun = generate_wake_gaussian(magVinf, wind_direction,
                                            turbine_xs[:, ite],
                                            turbine_ys[:, ite], turbine_z,
                                            hub_height, rotor_diameter;
                                            shear=shear, refH=refH)

        # Evaluate wake
        gt.calculate_field(fdom, wake_fun, "wake", "vector", "node";
                                                            raise_warn=false)

        # Save turbine and fluid domain vtks
        gt.save(windfarm, filepref; path=save_path, num=ite)
        gt.save(fdom, filepref*"_fdom"; path=save_path, num=ite)

    end

    # Previsualize farm
    if paraview
        if verbose
            println("\t"^(v_lvl)*"Calling Paraview...")
        end

        str = save_path*"/"
        for ti in 1:nturbines
            for bi in 1:nB
                str *= filepref*"_turbine$(ti)_rotor_blade$(bi)...vtk;"
            end
            str *= filepref*"_turbine$(ti)_rotor_hub...vtk;"
            str *= filepref*"_turbine$(ti)_tower...vtk;"
        end
        str *= filepref*"_fdom...vtk;";
        str *= filepref*"_perimeter.vtk;";

        run(`paraview --data=$str`)
    end



end
