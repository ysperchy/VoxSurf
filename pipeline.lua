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
models          = "models.txt"
model_file      = "model.stl"
voxel_size      = 0.5
overhang_angle  = 45
rotation_z      = 0
models_folder   = false

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
idx = ui_combo("Model", files_models)
model_file = files_models[idx+1]

-- Voxelizer
print("\nCalling voxelizer...\n")
t0 = os.time()
exec_code = os.execute("VoxSurf.exe".." -model "..model_file.." -mm "..voxel_size.." -angle "..overhang_angle.." -rotz "..rotation_z)
t1 = os.time()
exec_time = os.difftime(t1-t0)
if (exec_code ~= 0) then
  error("Voxelizer failed!\n")
end
print(
"\n\ndone!"..
"\nVoxelizer executed in "..exec_time.."s\n"
)

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

emit(Void)