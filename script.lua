-- Support shape based on type
-- A --> anchors: between an anchor voxel (white) and the object
-- B --> bars: horizontal bars
-- D --> diagonals: offsets and anchor merges
-- F --> for the anchor pillars that touch the object, they signal the connection between the pillar and the object (I have created a new voxel ID for this)
-- P --> pillars: vertical bars
function support_shape(x0,y0,z0,x1,y1,z1,mtx,t)
  v0 = mtx * v(x0,y0,z0)
  v1 = mtx * v(x1,y1,z1)
  vp = v1 - v0
  l = length(vp)
  r = 0.1
  if (t == 'A') then
    return translate(v0) * frame(vp) * cylinder(r,l)
  elseif (t == 'B') then
	return translate(v0) * frame(vp) * cylinder(r,l)
  elseif (t == 'D') then
    return translate(v0) * frame(vp) * cylinder(r,l)
  elseif (t == 'F') then
	return translate(v0) * frame(vp) * cylinder(r,l)
  elseif (t == 'P') then
	return translate(v0) * frame(vp) * cylinder(r,l)
  else
    return Void
  end
end

-- Load model
modelfile = 'model.stl'
model = load(Path..modelfile)
bx = bbox(model) -- bounding box

-- Load model bounding box extent (output from voxelizer)
transfile = 'out.slab.info'
trans = io.open(Path..transfile)
line = trans:read()
print(line)
m00,m10,m20,m30,m01,m11,m21,m31,m02,m12,m22,m32,m03,m13,m23,m33 = string.match(line, "(-?%d+.?%d*),(-?%d+.?%d*),(-?%d+.?%d*),(-?%d+.?%d*),(-?%d+.?%d*),(-?%d+.?%d*),(-?%d+.?%d*),(-?%d+.?%d*),(-?%d+.?%d*),(-?%d+.?%d*),(-?%d+.?%d*),(-?%d+.?%d*),(-?%d+.?%d*),(-?%d+.?%d*),(-?%d+.?%d*),(-?%d+.?%d*)") -- bounding box extent
m00 = tonumber(m00) m10 = tonumber(m10) m20 = tonumber(m20) m30 = tonumber(m30)
m01 = tonumber(m01) m11 = tonumber(m11) m21 = tonumber(m21) m31 = tonumber(m31)
m02 = tonumber(m02) m12 = tonumber(m12) m22 = tonumber(m22) m32 = tonumber(m32)
m03 = tonumber(m03) m13 = tonumber(m13) m23 = tonumber(m23) m33 = tonumber(m33)
obj2box = m{m00,m10,m20,m30,m01,m11,m21,m31,m02,m12,m22,m32,m03,m13,m23,m33}
box2obj = inverse(obj2box)

-- Read segments
segmentsfile = 'segments.csv'
segmentsA = {}
segmentsB = {}
segmentsD = {}
segmentsF = {}
segmentsP = {}
for line in io.lines(Path..segmentsfile) do
  x0,y0,z0,x1,y1,z1,t = string.match(line, "(%d+.?%d*),(%d+.?%d*),(%d+.?%d*),(%d+.?%d*),(%d+.?%d*),(%d+.?%d*),(%u)") -- read segment from file
  --print(x0..','..y0..','..z0..','..x1..','..y1..','..z1..','..t)
  x0 = tonumber(x0) y0 = tonumber(y0) z0 = tonumber(z0)
  x1 = tonumber(x1) y1 = tonumber(y1) z1 = tonumber(z1)
  x0 = 1.0 - x0 x1 = 1.0 - x1 -- flip X axis (due to MagicaVoxel weirdness)
  y0 = 1.0 - y0 y1 = 1.0 - y1 -- flip Y axis (due to MagicaVoxel weirdness)
  local shape = 
    support_shape(
	-- p0
	x0, y0, z0,
	-- p1
	x1, y1, z1,
	-- transformation matrix
	box2obj,
	-- type
	t)
	
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
emit(model,0)
emit(union(segmentsA),1)
emit(union(segmentsB),2)
emit(union(segmentsD),3)
emit(union(segmentsF),4)
emit(union(segmentsP),5)

-- Debug shapes for purposes of checking aligment with normalized voxel segments and model extent from voxelizer
--boundingboxV = translate(bx:min_corner()) * ocube(ex,ey,ez * 0.999)
--boundingboxM = translate(bx:min_corner()) * ocube(bx:extent().x,bx:extent().y,bx:extent().z * 0.999)
--barX = translate(bx:min_corner()) * ocube(ex,1,1)
--barY = translate(bx:min_corner()) * ocube(1,ey,1)
--barZ = translate(bx:min_corner()) * ocube(1,1,ez)
--emit(union{barX, barY, barZ})
--emit(boundingboxV)
--emit(boundingboxM)