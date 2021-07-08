--[[
 Salim PERCHY
 
 17/06/2021
 
 Copyright (c) 2021, INRIA
 All rights reserved.
 MFX Team
--]]

-- Disclaimer
print(
"\nCopyright (c) 2021, INRIA\n"..
"All rights reserved.\n"..
"\nSupports via Wave fuction collapse.\n\n"..
"MFX Team\n"
)

-- PARAMETERS
models           = "models.txt"
model_file       = "model.stl"
voxel_size       = 0.5
overhang_angle   = 45
rotation_z       = 0
no_model_lift    = false
models_folder    = true
batch_processing = true
stats_file       = "batch_stats.txt"
n_models         = 0
n_successes      = 0

output_voxelizer = "out"
wave_method_output= "segments.csv"
exec_time = 0

-- Create list of models
files_models = {}
for model in io.lines(Path..models) do
  if (models_folder) then
    files_models[#files_models+1] = 'models/'..model
  else
  files_models[#files_models+1] = model
  end
end

-- Tweaks
idx = ui_combo("Model", files_models)
model_file = files_models[idx+1]
voxel_size = ui_scalar("Voxel size (mm)", voxel_size, 0.1, 2.0)
rotation_z = ui_scalarBox("Rot Z (°)", rotation_z, 5)
overhang_angle = ui_number("Overhang angle (°)", overhang_angle, 0, 90)
no_model_lift = ui_bool("No bed gap", no_model_lift)

function pipeline()

  -- Statistics
  if (batch_processing) then
    io.write('Model: '..model_file..'\n')
    n_models = n_models + 1
  end

  -- Voxelizer
  print("\nCalling voxelizer...\n")
  t0 = os.time()
  lift_arg = ""
  if (no_model_lift) then
    lift_arg = " -nolifting"
  end
  exec_code = os.execute("VoxSurf.exe".." -model "..model_file.." -mm "..voxel_size.." -angle "..overhang_angle.." -rotz "..rotation_z..lift_arg)
  t1 = os.time()
  exec_time = os.difftime(t1-t0)
  if (exec_code ~= 0) then
    error("Voxelizer failed!\n")
  end
  print(
  "\n\ndone!"..
  "\nVoxelizer executed in "..exec_time.."s\n"
  )
  
  -- Statistics
  if (batch_processing) then
    if (exec_code ~= 0) then
      io.write("Voxelixer: failed\n")
      model_file = model_file - 1
    else
      io.write("Voxelixer: passed\n")
    end
    io.write("Voxelizer time: "..exec_time.."s\n")
  end

  if (exec_code == 0) then
  -- Wave function collapse
    print("\nCalling wave function collapse method...\n")
    t0 = os.time()
    exec_code = os.execute("ProcSupGen.exe".." "..output_voxelizer)
    t1 = os.time()
    exec_time = os.difftime(t1-t0)
    if (exec_code ~= 0) then
      error("Supports failed!\n")
    end
    print(
    "\n\ndone!"..
    "\nWave function collapse method executed in "..exec_time.."s\n"
    )
    exec_time_wave_func = exec_time
    
    -- Statistics
    if (batch_processing) then
      if (exec_code ~= 0) then
        io.write("Supports: failed\n")
      else
        io.write("Supports: passed\n")
        n_successes = n_successes + 1
      end
      io.write("Supports time: "..exec_time.."s\n")
    end
  end

  if (not batch_processing) then
    -- IceSL script
    print("\Calling IceSL script...\n")
    pipeline_call = true
    modelfile = model_file
    segmentsfile = wave_method_output
    rotZ = rotation_z
    t0 = os.time()
    dofile("script.lua")
    t1 = os.time()
    exec_time = os.difftime(t1-t0)
    print(
    "\n\ndone!"..
    "\nIceSL script executed in "..exec_time.."s\n"
    )
  end
  
  
  if (batch_processing) then
    io.write('---------------------------\n')
  end
  
end


if (batch_processing) then
  -- open statistics file
  batch_stats = io.open(Path..stats_file, "a")
  io.output(batch_stats)
  -- call pipeline (without IceSL rendering) for every model
  for _,filename in ipairs(files_models) do
    model_file = filename..'.stl'
    pipeline()
    os.execute('del 2p-init-out.slab.vox 2p-synth-out.slab.vox out.slab.info out.slab.points out.slab.vox segments.csv')
  end
  io.write('\n\n')
  io.write('N. of models processed: '..n_models..'\n')
  io.write('Supports succeeded: '..n_successes..'\n')
  io.write('Supports failed: '..(n_models - n_successes)..'\n')
  io.write('Success rate: '..math.floor(n_successes / n_models * 100.0)..'%\n')
  io.close(batch_stats)
else
 model_file = model_file..'.stl'
 pipeline()
end

emit(Void)