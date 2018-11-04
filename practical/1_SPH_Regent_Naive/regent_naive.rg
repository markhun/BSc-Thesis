import "regent"

local Config = require("sph_config")

local c = regentlib.c
local papi = terralib.includec("papi.h")
terralib.linklibrary("/usr/local/lib/libpapi.so")



fspace fs_point3d {
    { x, y, z } : double,
}

fspace fs_vel {
    { x, y, z } : double,
}

fspace fs_accel {
    { x, y, z } : double,
}

fspace fs_mass {
    mass : double,
}

fspace fs_density {
    density : double,
}

terra int_to_range(i : int, cfg : SPHConfig) : double
    var v : double = [double](i and ((1 << 15) - 1))
    v = v / (1 << 15)
    v = v - 0.5
    v = v * 2
    v = v * cfg.max_coord
    return v
end

function stall(x)
    for i=0,1000 do
        x = x*0.1
    end
    
    return x
end

function max2(x1, x2)
   return ( (x1 > x2) and x1 or x2)()
end

function max3(x1, x2, x3)
    return max2(max2(x1, x2), x3)
end

function min(x1, x2)
    return ( (x1 < x2 and x1 or x2) )()
end

function abs(x)
    -- testing if lua lib functions are available
    return math.abs(x)
end

-- //L_inf norm - the maximum of all distances
function L_inf(x1,y1,z1,x2,y2,z2)
	return math.max(abs(x1 - x2), abs(y1 - y2), abs(z1 - z2))
end
terraL_inf = terralib.cast( {double,double,double,double,double,double} -> double, L_inf)

-- //L_2 norm - the usual Euclidean norm
function L_2(x1, y1, z1, x2, y2, z2)
    return math.sqrt( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2) )
end
terraL_2 = terralib.cast( {double,double,double,double,double,double} -> double, L_2)

function monaghan_cubic_spline(h, distance)
	local distance = distance + stall(distance);
	local x = distance / h;
	if (x > 1) then return 0 end
	local res = 1 -- //8/(PI*h*h*h);
	if (x >= 0.5) then return res * 2 * (1 - x)*(1 - x)*(1 - x) end
	return res*(1 - 6 * x*x + 6 * x*x*x);
end
terraMonaghan_cubic_spline = terralib.cast( {double,double} -> double, monaghan_cubic_spline)

task update_density(point_region: region(ispace(int1d), fs_point3d),
                        density_region: region(ispace(int1d), fs_density),
                        accel_region: region(ispace(int1d), fs_accel),
                        h: double
                     )
where
    reads (point_region.x), reads (point_region.y), reads (point_region.z),
    reads writes(density_region.density),
    reads writes(accel_region.x), reads writes(accel_region.y), reads writes(accel_region.z),
    reads writes(density_region.density)
do
    var limits = point_region.bounds

    for i = [int](limits.lo), [int](limits.hi) + 1 do
        density_region[i].density = 0

        for j = [int](limits.lo), [int](limits.hi) + 1 do
            var L_inf = terraL_inf(point_region[i].x,point_region[i].y,point_region[i].z,point_region[j].x,point_region[j].y,point_region[j].z)
            if( (i ~= j) and ( L_inf < h ) ) then 
                var distance = terraL_2(point_region[i].x,point_region[i].y,point_region[i].z,point_region[j].x,point_region[j].y,point_region[j].z)
                if( distance < h) then
                    density_region[i].density = terraMonaghan_cubic_spline(h, distance)
                    
                end
            end

            accel_region[i].x = -0.1 * point_region[i].x
            accel_region[i].y = -0.1 * point_region[i].y
            accel_region[i].z = -0.1 * point_region[i].z

        end

    end

end

task update_by_pressure_force(point_region: region(ispace(int1d), fs_point3d),
                        density_region: region(ispace(int1d), fs_density),
                        mass_region: region(ispace(int1d), fs_mass),
                        accel_region: region(ispace(int1d), fs_accel),
                        h: double
                     )
where
    reads (point_region.x), reads (point_region.y), reads (point_region.z),
    reads (density_region.density), reads (mass_region.mass),
    reads writes(accel_region.x), reads writes(accel_region.y), reads writes(accel_region.z)
do
    var limits = accel_region.bounds

    for i = [int](limits.lo), [int](limits.hi) + 1 do

        for j = [int](limits.lo), [int](limits.hi) + 1 do
            var L_inf = terraL_inf(point_region[i].x,point_region[i].y,point_region[i].z,point_region[j].x,point_region[j].y,point_region[j].z)
            if( (i ~= j) and ( L_inf < h ) ) then 
                var distance = terraL_2(point_region[i].x,point_region[i].y,point_region[i].z,point_region[j].x,point_region[j].y,point_region[j].z)
                if( distance < h) then
                    
                    accel_region[i].x = accel_region[i].x + (point_region[i].x - point_region[j].x)*(1e-8*(density_region[i].density + density_region[j].density)*terraMonaghan_cubic_spline(h, distance)) / mass_region[i].mass;
                    accel_region[i].y = accel_region[i].y +  (point_region[i].y - point_region[j].y)*(1e-8*(density_region[i].density + density_region[j].density)*terraMonaghan_cubic_spline(h, distance)) / mass_region[i].mass;
                    accel_region[i].z = accel_region[i].z + (point_region[i].z - point_region[j].z)*(1e-8*(density_region[i].density + density_region[j].density)*terraMonaghan_cubic_spline(h, distance)) / mass_region[i].mass;

                end
            end

        end
    
    end

end

task integrate(point_region: region(ispace(int1d), fs_point3d),
                accel_region: region(ispace(int1d), fs_accel),
                vel_region: region(ispace(int1d), fs_vel),
                dt: double)
where 
    reads writes(point_region.x), reads writes(point_region.y), reads writes(point_region.z),
    reads writes(accel_region.x), reads writes(accel_region.y), reads writes(accel_region.z),
    reads writes(vel_region.x), reads writes(vel_region.y), reads writes(vel_region.z)
do
    var limits = point_region.bounds

    for i = [int](limits.lo), [int](limits.hi) + 1 do
        vel_region[i].x = vel_region[i].x + accel_region[i].x * dt;
        vel_region[i].y = vel_region[i].y + accel_region[i].y * dt;
        vel_region[i].z = vel_region[i].z + accel_region[i].z * dt;
        point_region[i].x = point_region[i].x + vel_region[i].x * dt;
        point_region[i].y = point_region[i].y + vel_region[i].y * dt;
        point_region[i].z = point_region[i].z + vel_region[i].z * dt;
        accel_region[i].x = 0;
        accel_region[i].y = 0;
        accel_region[i].z = 0;
    end

end

task init_points(point_region: region(ispace(int1d), fs_point3d), cfg: SPHConfig )
where
    reads writes(point_region.x), reads writes(point_region.y), reads writes(point_region.z)
do   
    for p in point_region do
        p.x = int_to_range(c.rand(), cfg)
        p.y = int_to_range(c.rand(), cfg)
        p.z = int_to_range(c.rand(), cfg)
    end
end


task printer_points( point_region: region(ispace(int1d), fs_point3d) )
where
    reads(point_region.x), reads(point_region.y), reads(point_region.z)
do
    var limits = point_region.bounds

    for i = [int](limits.lo), [int](limits.hi) + 1 do
        c.printf("X: %f \t Y: %f \t Z: %f \n", point_region[i].x, point_region[i].y, point_region[i].z)
    end
    c.printf("  --------------------------------  \n\n")
end

task printer_accel( accel_region: region(ispace(int1d), fs_accel) )
where
    reads(accel_region.x), reads(accel_region.y), reads(accel_region.z)
do
    var limits = accel_region.bounds

    for i = [int](limits.lo), [int](limits.hi) + 1 do
        c.printf("X: %f \t Y: %f \t Z: %f \n", accel_region[i].x, accel_region[i].y, accel_region[i].z)
    end
end

-- insert artificial dependencies to ensure that this executed last
task time_measurement(point_region: region(ispace(int1d), fs_point3d)) : int64
where
    reads writes(point_region)
do
    return papi.PAPI_get_real_usec();
end

task main()

    var times : int64
    var timee : int64
    var cfg : SPHConfig
    cfg:initialize_from_command()

    c.printf("-- Started Similuation with values:\t--\n")
    c.printf("-- particle_count:\t %u \t\t--\n", cfg.particle_count)
    c.printf("-- iteration_count:\t %u \t\t--\n", cfg.iteration_count)
    c.printf("-- max_coord:\t\t %f \t--\n", cfg.max_coord)
    c.printf("-- h:\t\t\t %f \t--\n", cfg.h)
    c.printf("-- dt:\t\t\t %f \t--\n", cfg.dt)
    c.printf("------------------------------------------\n\n")


    var structured_particle_is = ispace(int1d, cfg.particle_count, 0)

    var point_region = region(structured_particle_is, fs_point3d)
    var vel_region = region(structured_particle_is, fs_vel)
    var accel_region = region(structured_particle_is, fs_accel)
    var mass_region = region(structured_particle_is, fs_mass)
    var density_region = region(structured_particle_is, fs_density)


    init_points(point_region, cfg)

    fill(vel_region.x, 0)
    fill(vel_region.y, 0)
    fill(vel_region.z, 0)
    fill(accel_region.x, 0)
    fill(accel_region.y, 0)
    fill(accel_region.z, 0)
    fill(mass_region.mass, 1)
    fill(density_region.density, 0)

    -- Init finished
    -- printer_points(point_region)

    -- Starttime measured by PAPI
	times = time_measurement(point_region);	

    for t = 0, cfg.iteration_count do
        update_density(point_region,
                    density_region,
                    accel_region,
                    cfg.h)
        
        update_by_pressure_force(point_region,
                    density_region,
                    mass_region,
                    accel_region,
                    cfg.h)

        integrate(point_region,
                    accel_region,
                    vel_region,
                    cfg.dt)
    end
    
    -- BEWARE !!! Calling legion runtime routines like this barrier breaks Legion Spy !!!
    -- c.legion_runtime_issue_execution_fence(__runtime(), __context())
    -- c.printf("-- Finished Similuation --\n\n")
   
    timee = time_measurement(point_region);	
	c.printf("Time in microseconds: %lld\n",timee-times);

    -- printer_points(point_region)
    -- printer_accel(accel_region)

end

regentlib.start(main)