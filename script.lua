--[[
 Salim PERCHY
 
 17/06/2021
 
 Copyright (c) 2021, INRIA
 All rights reserved.
 MFX Team
--]]

-- Global variables
total_length = 0
radius_supports = 0.3
anchor_extension = 0

-- helper function
function find_closest_point(point, points, h)
  local distance = math.huge
  for _,p in ipairs(points) do
    local d = length(p - point)
	  if (d < distance) then
		distance = d
		result = p
	  end
  end
  return result
end

-- Support shape based on type
-- A --> anchors: between an anchor voxel and the object
-- B --> bars: horizontal bars
-- D --> diagonals: offsets and anchor merges
-- F --> for the anchor pillars that touch the object, they signal the connection between the pillar and the object
-- P --> pillars: vertical bars
function support_shapes(x0,y0,z0,x1,y1,z1,mtx,t)
  v0 = mtx * v(x0,y0,z0)
  v1 = mtx * v(x1,y1,z1)
  vp = v1 - v0
  l = length(vp)
  total_length = total_length + l
  r = radius_supports
  if (t == 'A') then
    local v2 = find_closest_point(v1, pointsArray, 0, t)
	local v3 = normalize(vp) * anchor_extension
    return {translate(v0) * sphere(r), cone(r, 0, v0, v2 + v3)}
  elseif (t == 'B') then
	return {translate(v0) * frame(vp) * translate(0, 0, l/2) * ccube(r*2, 0.4, l)}
  elseif (t == 'D') then
    return {translate(v0) * sphere(r), translate(v0) * frame(vp) * cylinder(r, l), translate(v0) * translate(vp) * sphere(r)}
  elseif (t == 'F') then
    local v2 = find_closest_point(v1, pointsArray, 0, t)
	local v3 = normalize(vp) * anchor_extension
	return {translate(v0) * sphere(r), cone(r, 0, v0, v2 + v3)}
  elseif (t == 'P') then
	return {translate(v0) * frame(vp) * translate(0, 0, l/2) * ccube(r*2, r*2, l)}
  end
end

-- Load model
if (not pipeline_call) then
  modelfile = 'model.stl'
end
model = load(Path..modelfile)

-- Load model transformation from voxelizer output
transfile = 'out.slab.info'
trans = io.open(Path..transfile)
line = trans:read()
m00,m10,m20,m30,m01,m11,m21,m31,m02,m12,m22,m32,m03,m13,m23,m33 = string.match(line, "(-?%d+.?%d*),(-?%d+.?%d*),(-?%d+.?%d*),(-?%d+.?%d*),(-?%d+.?%d*),(-?%d+.?%d*),(-?%d+.?%d*),(-?%d+.?%d*),(-?%d+.?%d*),(-?%d+.?%d*),(-?%d+.?%d*),(-?%d+.?%d*),(-?%d+.?%d*),(-?%d+.?%d*),(-?%d+.?%d*),(-?%d+.?%d*)") -- bounding box extent
m00 = tonumber(m00) m10 = tonumber(m10) m20 = tonumber(m20) m30 = tonumber(m30)
m01 = tonumber(m01) m11 = tonumber(m11) m21 = tonumber(m21) m31 = tonumber(m31)
m02 = tonumber(m02) m12 = tonumber(m12) m22 = tonumber(m22) m32 = tonumber(m32)
m03 = tonumber(m03) m13 = tonumber(m13) m23 = tonumber(m23) m33 = tonumber(m33)
obj2box = m{m00,m10,m20,m30,m01,m11,m21,m31,m02,m12,m22,m32,m03,m13,m23,m33}
box2obj = inverse(obj2box)

-- Load voxel poins from voxelizer output
pointsfile = 'out.slab.points'
pointsArray = {}
for line in io.lines(Path..pointsfile) do
  local x,y,z = string.match(line, "(-?%d+.?%d*),(-?%d+.?%d*),(-?%d+.?%d*)")
  pointsArray[#pointsArray + 1] = v(tonumber(x),tonumber(y),tonumber(z))
end

-- Read segments
if (not pipeline_call) then
  segmentsfile = 'segments.csv'
end
segmentsA = {}
segmentsB = {}
segmentsD = {}
segmentsF = {}
segmentsP = {}
supports  = {}

for line in io.lines(Path..segmentsfile) do

  x0,y0,z0,x1,y1,z1,t = string.match(line, "(%d+.?%d*),(%d+.?%d*),(%d+.?%d*),(%d+.?%d*),(%d+.?%d*),(%d+.?%d*),(%u)") -- read segment from file
  --print(x0..','..y0..','..z0..','..x1..','..y1..','..z1..','..t)
  x0 = tonumber(x0) y0 = tonumber(y0) z0 = tonumber(z0)
  x1 = tonumber(x1) y1 = tonumber(y1) z1 = tonumber(z1)
  x0 = 1.0 - x0 x1 = 1.0 - x1 -- flip X axis (due to MagicaVoxel weirdness)
  y0 = 1.0 - y0 y1 = 1.0 - y1 -- flip Y axis (due to MagicaVoxel weirdness)

  local shapes =
    support_shapes(
      -- p0
      x0, y0, z0,
      -- p1
      x1, y1, z1,
      -- transformation matrix
      box2obj,
      -- type
      t)

  for _,s in ipairs(shapes) do
    supports[#supports + 1] = s
  end

  if (t == 'A') then
    segmentsA[#segmentsA + 1] = shape
  elseif (t == 'B') then
	  segmentsB[#segmentsB + 1] = shape
  elseif (t == 'D') then
    segmentsD[#segmentsD + 1] = shape
  elseif (t == 'F') then
	  segmentsF[#segmentsF + 1] = shape
  elseif (t == 'P') then
	  segmentsP[#segmentsP + 1] = shape
  else
    error("Unknown type!")
  end
end

-- Output shape and supports
scaleFactor = 1
if (not pipeline_call) then
  rotZ = 0
end
emit(scale(scaleFactor) * rotate(rotZ, Z) * model, 0)
--[[
for i = 1,#pointsArray,10 do
  emit(translate(pointsArray[i].x,pointsArray[i].y,pointsArray[i].z)*cube(0.2), 2)
end
]]--
emit(scale(scaleFactor) * union(supports), 1)

print('Total lenght of supports: '..math.floor(total_length)..'mm')

-- Print settings
set_setting_value('num_shells_1', 0)
set_setting_value('add_brim', true)
set_setting_value('brim_distance_to_print_mm', 0)
set_setting_value('brim_num_contours', 8)
set_setting_value('enable_z_lift', true)
set_setting_value('is_support_1', true)