function support_shape(x0,y0,z0,x1,y1,z1,t)
  v0 = v(x0,y0,z0)
  v1 = v(x1,y1,z1)
  vp = v1 - v0
  l = length(vp)
  r = 0.5
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
model = load(Path..'model.stl')
bx = bbox(model)

-- Load model bounding box extent
trans = io.open(Path..'out.slab.info')
line = trans:read()
ex,ey,ez = string.match(line, "(%d+.?%d*),(%d+.?%d*),(%d+.?%d*)")
trans:close()
ex = tonumber(ex)
ey = tonumber(ey)
ez = tonumber(ez)

-- Read segments
segmentsA = {}
segmentsB = {}
segmentsD = {}
segmentsF = {}
segmentsP = {}
for line in io.lines(Path..'segments.csv') do
  x0,y0,z0,x1,y1,z1,t = string.match(line, "(%d+.?%d*),(%d+.?%d*),(%d+.?%d*),(%d+.?%d*),(%d+.?%d*),(%d+.?%d*),(%u)")
  --print(x0..','..y0..','..z0..','..x1..','..y1..','..z1..','..t)
  x0 = tonumber(x0) y0 = tonumber(y0) z0 = tonumber(z0)
  x1 = tonumber(x1) y1 = tonumber(y1) z1 = tonumber(z1)
  x0 = 1.0 - x0 y0 = 1.0 - y0
  x1 = 1.0 - x1 y1 = 1.0 - y1
  local shape = translate(bx:min_corner()) * 
    support_shape(ex * x0, ey * y0, ez * z0, ex * x1, ey * y1, ez * z1, t)
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

-- Output shape
emit(model,0)
emit(union(segmentsA),1)
emit(union(segmentsB),2)
emit(union(segmentsD),3)
emit(union(segmentsF),4)
emit(union(segmentsP),5)

--[[
boundingbox = translate(bx:min_corner()) * ocube(ex,ey,ez)
barX = translate(bx:min_corner()) * ocube(ex,1,1)
barY = translate(bx:min_corner()) * ocube(1,ey,1)
barZ = translate(bx:min_corner()) * ocube(1,1,ez)
emit(union{barX, barY, barZ})
emit(boundingbox)
]]--