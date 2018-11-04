import "regent"

local Config = require("sph_config")

local c = regentlib.c
local sqrt = regentlib.sqrt(double)
local cmath = terralib.includecstring [[
   #include <math.h>
]]

fspace fs_vel {
    { x, y, z } : double,
}

fspace fs_accel {
    { x, y, z } : double,
}

fspace fs_point3d {
    { x, y, z } : double,
    vel : fs_vel,
    accel : fs_accel,
    mass : double,
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

terra stall(x : double) : double
    for i=0,1000 do
        x = x*0.1
    end
    
    return x
end

terra max2(x1: double, x2: double) : double
    return cmath.fmax(x1,x2)
end

terra max3(x1: double, x2: double, x3: double) : double  
    var local_max : double = cmath.fmax(x1,x2)
    if x3 > local_max then
        return x3 
    else
        return local_max
    end
end

terra min(x1: double, x2: double)
    if x1 < x2 then
        return x1 
    else
        return x2
    end 
end

terra abs(x:double) : double
    -- testing if c lib functions are available
    return cmath.fabs(x)
end

-- //L_inf norm - the maximum of all distances
terra L_inf(x1: double,y1: double,z1: double,x2: double, y2: double, z2: double) : double
	return max3(abs(x1 - x2), abs(y1 - y2), abs(z1 - z2))
end

-- //L_2 norm - the usual Euclidean norm
terra L_2(x1 : double, y1 : double, z1 : double, x2 : double, y2 : double, z2 : double) : double
    return sqrt( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2) )
end

terra monaghan_cubic_spline(h: double, distance: double) : double
	var distance : double = distance + stall(distance);
	var x : double = distance / h;
	if (x > 1) then return 0 end
	var res : double = 1 -- //8/(PI*h*h*h);
	if (x >= 0.5) then return res * 2 * (1 - x)*(1 - x)*(1 - x) end
	return res*(1 - 6 * x*x + 6 * x*x*x);
end

task update_density(point_region: region(ispace(int1d), fs_point3d), h: double)
where
    reads (point_region.x), reads (point_region.y), reads (point_region.z),
    reads writes(point_region.density),
    reads writes(point_region.accel)
do
    var limits = point_region.bounds

    for i = [int](limits.lo), [int](limits.hi) + 1 do
        point_region[i].density = 0

        for j = [int](limits.lo), [int](limits.hi) + 1 do
            var L_inf = L_inf(point_region[i].x,point_region[i].y,point_region[i].z,point_region[j].x,point_region[j].y,point_region[j].z)
            if( (i ~= j) and ( L_inf < h ) ) then 
                var distance = L_2(point_region[i].x,point_region[i].y,point_region[i].z,point_region[j].x,point_region[j].y,point_region[j].z)
                if( distance < h) then
                    point_region[i].density = monaghan_cubic_spline(h, distance)
                end
            end

        end

        point_region[i].accel.x = -0.1 * point_region[i].x
        point_region[i].accel.y = -0.1 * point_region[i].y
        point_region[i].accel.z = -0.1 * point_region[i].z

    end

end

task update_by_pressure_force(point_region: region(ispace(int1d), fs_point3d), h: double)
where
    reads (point_region.x), reads (point_region.y), reads (point_region.z),
    reads (point_region.density), reads (point_region.mass),
    reads writes(point_region.accel)
do
    var limits = point_region.bounds

    for i = [int](limits.lo), [int](limits.hi) + 1 do

        for j = [int](limits.lo), [int](limits.hi) + 1 do
            var L_inf = L_inf(point_region[i].x,point_region[i].y,point_region[i].z,point_region[j].x,point_region[j].y,point_region[j].z)
            if( (i ~= j) and ( L_inf < h ) ) then 
                var distance = L_2(point_region[i].x,point_region[i].y,point_region[i].z,point_region[j].x,point_region[j].y,point_region[j].z)
                if( distance < h) then
                    
                    point_region[i].accel.x = point_region[i].accel.x + (point_region[i].x - point_region[j].x)*(1e-8*(point_region[i].density + point_region[j].density)*monaghan_cubic_spline(h, distance)) / point_region[i].mass;
                    point_region[i].accel.y = point_region[i].accel.y + (point_region[i].y - point_region[j].y)*(1e-8*(point_region[i].density + point_region[j].density)*monaghan_cubic_spline(h, distance)) / point_region[i].mass;
                    point_region[i].accel.z = point_region[i].accel.z + (point_region[i].z - point_region[j].z)*(1e-8*(point_region[i].density + point_region[j].density)*monaghan_cubic_spline(h, distance)) / point_region[i].mass;

                end
            end

        end
    
    end

end

task integrate(point_region: region(ispace(int1d), fs_point3d), dt: double)
where 
    reads writes(point_region.x), reads writes(point_region.y), reads writes(point_region.z),
    reads writes(point_region.accel),
    reads writes(point_region.vel)
do
    var limits = point_region.bounds

    for i = [int](limits.lo), [int](limits.hi) + 1 do
        point_region[i].vel.x =  point_region[i].vel.x + point_region[i].accel.x * dt;
        point_region[i].vel.y =  point_region[i].vel.y + point_region[i].accel.y * dt;
        point_region[i].vel.z =  point_region[i].vel.z + point_region[i].accel.z * dt;
        point_region[i].x = point_region[i].x + point_region[i].vel.x * dt;
        point_region[i].y = point_region[i].y + point_region[i].vel.y * dt;
        point_region[i].z = point_region[i].z + point_region[i].vel.z * dt;
        point_region[i].accel.x = 0;
        point_region[i].accel.y = 0;
        point_region[i].accel.z = 0;
    end

end

task init_points(point_region: region(ispace(int1d), fs_point3d), cfg: SPHConfig )
where
    reads writes(point_region.x), reads writes(point_region.y), reads writes(point_region.z)
do
    
    -- for i : int1d(structured_particle_is) in structured_particle_is do
    --    point_region[i].x = int_to_range(c.rand(), cfg)
    --    point_region[i].y = int_to_range(c.rand(), cfg)
    --    point_region[i].z = int_to_range(c.rand(), cfg)
    
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


task main()

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

    -- build a region of all points, this one will be partitioned and distributed
    var point_region = region(structured_particle_is, fs_point3d)
    -- partition the region
    var point_region_partition_small = partition(equal, point_region, ispace(int1d, 2))
    

    for color in point_region_partition_small.colors do
        init_points(point_region_partition_small[color], cfg)
    end
    fill(point_region.vel, {0,0,0})
    fill(point_region.accel, {0,0,0})
    fill(point_region.mass, 1)
    fill(point_region.density, 0)
    -- Init finished


    for t = 0, cfg.iteration_count do
        -- The following loops produce incorrect results!
        -- Partitioning can lead to partial iterations over datasets,
        -- where a iteration over all the data would be necessary
        for color in point_region_partition_small.colors do
            update_density(point_region_partition_small[color], cfg.h)
            -- c.printf(" Updated densities --\n\n")
        end
        
        for color in point_region_partition_small.colors do
            update_by_pressure_force(point_region_partition_small[color], cfg.h)
            -- c.printf(" Updated pressure forces --\n\n")
        end
        
        for color in point_region_partition_small.colors do
            integrate(point_region_partition_small[color], cfg.dt)
            -- c.printf(" Integrated --\n\n")
        end
        -- c.printf("-- Finished 1 Iteration --\n\n")
    end
    
    -- BEWARE !!! Calling legion runtime routines like this barrier breaks Legion Spy !!!
    -- c.legion_runtime_issue_execution_fence(__runtime(), __context())
    -- c.printf("-- Finished Similuation --\n\n")

    -- printer_points(point_region)
    -- printer_accel(accel_region)

end

regentlib.start(main)