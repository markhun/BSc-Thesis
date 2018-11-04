import "regent"

local c = regentlib.c
local cstr = terralib.includecstring [[
   #include <string.h>
]]

struct SPHConfig {
    particle_count : uint32,
	iteration_count : uint32,
	iteration : uint32,
	chunk_size : uint32,
	max_coord : double,
	h : double,
	dt : double,
  }



terra print_usage_and_abort()
    c.printf("Usage: .../regent.py sph.rg [OPTIONS]\n")
    c.printf("OPTIONS\n")
    c.printf("  -h  \t\t\t\t: Print the usage and exit.\n")
    c.printf("  -c {particle count}\n")
    c.printf("  -i {iteration count}\n")
    c.printf("  -m {maximum coordinates}\t: Sets the space boundaries to [-max;max]\n")
    c.abort()
end

terra SPHConfig:initialize_from_command()

    self.particle_count = 1000
    self.iteration_count = 2
    self.iteration = 0
    -- chunk_size = particle_count / (4 * 240),
    self.chunk_size = 60 / (4 * 240)
    self.max_coord = 15
    self.h = 10
    self.dt = 0.1

    var args = c.legion_runtime_get_input_args()
    var i = 1
    while i < args.argc do
    if cstr.strcmp(args.argv[i], "-h") == 0 then
        print_usage_and_abort()
    elseif cstr.strcmp(args.argv[i], "-c") == 0 then
        i = i + 1
        self.particle_count = c.strtol(args.argv[i], nil, 10)
    elseif cstr.strcmp(args.argv[i], "-i") == 0 then
        i = i + 1
        self.iteration_count = c.strtol(args.argv[i], nil, 10)
    elseif cstr.strcmp(args.argv[i], "-m") == 0 then
        i = i + 1
        self.iteration_count = c.strtol(args.argv[i], nil, 10)
    end
    i = i + 1
    end

end

return SPHConfig