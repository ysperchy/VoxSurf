// --------------------------------------------------------------
// SL 2018-01-02
// A simple, easily hackable CPU surface voxelizer
// MIT-license
// (c) Sylvain Lefebvre, https://github.com/sylefeb
// --------------------------------------------------------------

/*

Takes as input a file 'model.stl' from the source directory.
Outputs a voxel file named 'out.slab.vox' that can be imported 
by 'MagicaVoxel' https://ephtracy.github.io/

Change VOXEL_RESOLUTION  to fit your needs.
Set    VOXEL_FILL_INSIDE to 1 to fill in the interior
Set    VOXEL_ROBUST_FILL to 1 to fill in the interior using 
       a voting scheme (more robust, slower)

The basic principle is to rasterize triangles using three 2D axis
aligned grids, using integer arithmetic (fixed floating point)
for robust triangle interior checks.

Very simple and quite efficient despite a straightforward implementation.
Higher resolutions could easily be reached by not storing the
voxels as a 3D array of booleans (e.g. use blocking or an octree).

For the inside fill to work properly, the mesh has to be perfectly
watertight, with exactly matching vertices between neighboring 
vertices.

*/

#include "path.h"

#include <LibSL/LibSL.h>

#include <iostream>
#include <algorithm>
#include <queue>
#include <fstream>


// --------------------------------------------------------------

#define VOXEL_RESOLUTION  128
#define VOXEL_FILL_INSIDE 1
#define VOXEL_ROBUST_FILL 0

// --------------------------------------------------------------

#define FP_POW    16
#define FP_SCALE  (1<<FP_POW)
#define BOX_SCALE v3f(VOXEL_RESOLUTION*FP_SCALE)

#define ALONG_X  1
#define ALONG_Y  2
#define ALONG_Z  4
#define INSIDE   8
#define INSIDE_X 16
#define INSIDE_Y 32
#define INSIDE_Z 64
#define OVERHANG 128

#define OVERHANG_ANGLE 45
#define GRAVITY_VECTOR v3f(0,0,1)
#define BED_LEVEL 2
#define FILTER_PARITY 1
#define PARITY_RULE 1 // 0: odd, 1: even
#define FILTER_SPHERE 1
#define SPHERE_RADIUS 3.0f
#define FILTER_COLLISION 1
#define COLLISION_RETRACT 0

// --------------------------------------------------------------

// saves a voxel file (.slab.vox format, can be imported by MagicaVoxel)
void saveAsVox(const char *fname, const Array3D<uchar>& voxs)
{
  Array<v3b> palette(256); // RGB palette
  palette.fill(0);
  palette[0  ] = v3b(  0,   0, 255); // pillar
  palette[252] = v3b(255, 255, 255); // anchor
  palette[253] = v3b(  0,   0,   0); // solid
  std::ofstream f(fname, std::ios::binary);
  sl_assert(f.is_open())
  long sx = voxs.xsize(), sy = voxs.ysize(), sz = voxs.zsize();
  f.write(reinterpret_cast<char*>(&sx), 4 * 1);
  f.write(reinterpret_cast<char*>(&sy), 4 * 1);
  f.write(reinterpret_cast<char*>(&sz), 4 * 1);
  ForRangeReverse(i, sx - 1, 0) {
    ForIndex(j, sy) {
      ForRangeReverse(k, sz - 1, 0) {
        uchar v   = voxs.at(i, j, k);
        uchar pal = v != 0 ? 253 : 255;
        if (v == INSIDE) {
          pal = 253;
        }
        if (v & OVERHANG) {
          pal = 252;
        }
        f.write(reinterpret_cast<char*>(&pal), sizeof(uchar) * 1);
      }
    }
  }
  f.write(reinterpret_cast<char*>(palette.raw()), sizeof(v3b) * 256);
}

// --------------------------------------------------------------

inline bool isInTriangle(int i, int j, const v3i& p0, const v3i& p1, const v3i& p2, int& _depth)
{
  v2i delta_p0 = v2i(i, j) - v2i(p0);
  v2i delta_p1 = v2i(i, j) - v2i(p1);
  v2i delta_p2 = v2i(i, j) - v2i(p2);
  v2i delta10 = v2i(p1) - v2i(p0);
  v2i delta21 = v2i(p2) - v2i(p1);
  v2i delta02 = v2i(p0) - v2i(p2);

  int64_t c0 = (int64_t)delta_p0[0] * (int64_t)delta10[1] - (int64_t)delta_p0[1] * (int64_t)delta10[0];
  int64_t c1 = (int64_t)delta_p1[0] * (int64_t)delta21[1] - (int64_t)delta_p1[1] * (int64_t)delta21[0];
  int64_t c2 = (int64_t)delta_p2[0] * (int64_t)delta02[1] - (int64_t)delta_p2[1] * (int64_t)delta02[0];
  bool inside = (c0 <= 0 && c1 <= 0 && c2 <= 0) || (c0 >= 0 && c1 >= 0 && c2 >= 0);

  if (inside) {
    int64_t area = c0 + c1 + c2;
    int64_t b0 = (c1 << 10) / area;
    int64_t b1 = (c2 << 10) / area;
    int64_t b2 = (1 << 10) - b0 - b1;
    _depth = (int)((b0 * p0[2] + b1 * p1[2] + b2 * p2[2]) >> 10);
  }
  return inside;
}

// --------------------------------------------------------------

class swizzle_xyz
{
public:
  inline v3i forward(const v3i& v)  const { return v; }
  inline v3i backward(const v3i& v) const { return v; }
  inline int along() const { return ALONG_Z; }
};

class swizzle_zxy
{
public:
  inline v3i   forward(const v3i& v)  const { return v3i(v[2], v[0], v[1]); }
  inline v3i   backward(const v3i& v) const { return v3i(v[1], v[2], v[0]); }
  inline uchar along() const { return ALONG_Y; }
};

class swizzle_yzx
{
public:
  inline v3i   forward(const v3i& v)  const { return v3i(v[1], v[2], v[0]); }
  inline v3i   backward(const v3i& v) const { return v3i(v[2], v[0], v[1]); }
  inline uchar along() const { return ALONG_X; }
};

// --------------------------------------------------------------

template <class S>
void rasterize(
  const v3u&                  tri,
  const std::vector<v3i>&     pts,
  Array3D<uchar>&             _voxs)
{
  const S swizzler;
  v3i tripts[3] = {
    swizzler.forward(pts[tri[0]]),
    swizzler.forward(pts[tri[1]]),
    swizzler.forward(pts[tri[2]])
  };
  // check if triangle is valid
  v2i delta10 = v2i(tripts[1]) - v2i(tripts[0]);
  v2i delta21 = v2i(tripts[2]) - v2i(tripts[1]);
  v2i delta02 = v2i(tripts[0]) - v2i(tripts[2]);
  if (delta10 == v2i(0)) return;
  if (delta21 == v2i(0)) return;
  if (delta02 == v2i(0)) return;
  if (delta02[0] * delta10[1] - delta02[1] * delta10[0] == 0) return;
  // check if triangle is overhanging
  v3f v10 = v3f(pts[tri[1]] - pts[tri[0]]);
  v3f v20 = v3f(pts[tri[2]] - pts[tri[0]]);
  v3f nrml = cross(v20, v10);
  nrml = normalize_safe(nrml);
  float cosAngle = dot(nrml, GRAVITY_VECTOR);
  float angleRad = std::acosf(cosAngle); // return value on interval [0 - pi]
  float angleDeg = angleRad * 180.0f / static_cast<float>(M_PI);
  bool overhang = angleDeg < static_cast<float>(OVERHANG_ANGLE);
  // proceed
  AAB<2, int> pixbx;
  pixbx.addPoint(v2i(tripts[0]) / FP_SCALE);
  pixbx.addPoint(v2i(tripts[1]) / FP_SCALE);
  pixbx.addPoint(v2i(tripts[2]) / FP_SCALE);
  for (int j = pixbx.minCorner()[1]; j <= pixbx.maxCorner()[1]; j++) {
    for (int i = pixbx.minCorner()[0]; i <= pixbx.maxCorner()[0]; i++) {
      int depth;
      if (isInTriangle(
        (i << FP_POW) + (1 << (FP_POW - 1)), // centered
        (j << FP_POW) + (1 << (FP_POW - 1)), // centered
        tripts[0], tripts[1], tripts[2], depth)) {
        v3i vx = swizzler.backward(v3i(i, j, depth >> FP_POW));
        // tag the voxel as occupied
        // NOTE: voxels are likely to be hit multiple times (e.g. thin features)
        //       we flip the bit every time a hit occurs in a voxel
        _voxs.at(vx[0], vx[1], vx[2]) = ( _voxs.at(vx[0], vx[1], vx[2]) ^ swizzler.along() );
        // overhang mark
        if (overhang && vx[2] > BED_LEVEL) {
          _voxs.at(vx[0], vx[1], vx[2]) |= OVERHANG;
        }
      }
    }
  }
}

// --------------------------------------------------------------

// This version is more robust by using all three direction
// and voting among them to decide what is filled or not
void fillInsideVoting(Array3D<uchar>& _voxs)
{
  // along x
  ForIndex(k, _voxs.zsize()) {
    ForIndex(j, _voxs.ysize()) {
      bool inside = false;
      ForIndex(i, _voxs.xsize()) {
        if (_voxs.at(i, j, k) & ALONG_X) {
          inside = !inside;
        }
        if (inside) {
          _voxs.at(i, j, k) |= INSIDE_X;
        }
      }
    }
  }
  // along y
  ForIndex(k, _voxs.zsize()) {
    ForIndex(j, _voxs.xsize()) {
      bool inside = false;
      ForIndex(i, _voxs.ysize()) {
        if (_voxs.at(j, i, k) & ALONG_Y) {
          inside = !inside;
        }
        if (inside) {
          _voxs.at(j, i, k) |= INSIDE_Y;
        }
      }
    }
  }
  // along z
  ForIndex(k, _voxs.ysize()) {
    ForIndex(j, _voxs.xsize()) {
      bool inside = false;
      ForIndex(i, _voxs.zsize()) {
        if (_voxs.at(j, k, i) & ALONG_Z) {
          inside = !inside;
        }
        if (inside) {
          _voxs.at(j, k, i) |= INSIDE_Z;
        }
      }
    }
  }
  // voting
  ForArray3D(_voxs, i, j, k) {
    uchar v = _voxs.at(i, j, k);
    int votes =
      (  (v & INSIDE_X) ? 1 : 0)
      + ((v & INSIDE_Y) ? 1 : 0)
      + ((v & INSIDE_Z) ? 1 : 0);
    // clean
    _voxs.at(i, j, k) &= ~(INSIDE_X | INSIDE_Y | INSIDE_Z);
    if (votes > 1) {
      // tag as inside
      _voxs.at(i, j, k) |= INSIDE;
    }
  }
}

// --------------------------------------------------------------

void fillInside(Array3D<uchar>& _voxs)
{
  ForIndex(k, _voxs.zsize()) {
    ForIndex(j, _voxs.ysize()) {
      bool inside = false;
      ForIndex(i, _voxs.xsize()) {
        if (_voxs.at(i, j, k) & ALONG_X) {
          inside = !inside;
        }
        if (inside) {
          _voxs.at(i, j, k) |= INSIDE;
        }
      }
    }
  }
}

// --------------------------------------------------------------

void insideFilter(Array3D<uchar>& _voxs)
{
  ForIndex(k, _voxs.zsize()) {
    ForIndex(j, _voxs.ysize()) {
      ForIndex(i, _voxs.xsize()) {
        if ((_voxs.at(i, j, k) & OVERHANG) && (_voxs.at(i, j, k - 1) & (ALONG_X | ALONG_Y | ALONG_Z))) {
          _voxs.at(i, j, k) &= ~OVERHANG;
        }
      }
    }
  }
}

// --------------------------------------------------------------

void parityFilter(Array3D<uchar>& _voxs)
{
  ForIndex(k, _voxs.zsize()) {
    ForIndex(j, _voxs.ysize()) {
      ForIndex(i, _voxs.xsize()) {
        if (!(i % 2 == PARITY_RULE && j % 2 == PARITY_RULE)) {
          _voxs.at(i, j, k) &= ~OVERHANG;
        }
      }
    }
  }
}

// --------------------------------------------------------------

void sphereFilter(Array3D<uchar>& _voxs)
{
  const int r = static_cast<int>(std::ceilf(SPHERE_RADIUS));
  ForIndex(k, _voxs.zsize()) {
    ForIndex(j, _voxs.ysize()) {
      ForIndex(i, _voxs.xsize()) {
        if (_voxs.at(i, j, k) & OVERHANG) { // sphere center
          ForRange(c, std::max(0, k - r), std::min(k + r, static_cast<int>(_voxs.zsize()) - 1)) {
            ForRange(b, std::max(0, j - r), std::min(j + r, static_cast<int>(_voxs.ysize()) - 1)) {
              ForRange(a, std::max(0, i - r), std::min(i + r, static_cast<int>(_voxs.xsize()) - 1)) {
                if (static_cast<float>(sqLength(v3i(i-a,j-b,k-c))) < SPHERE_RADIUS * SPHERE_RADIUS) { // inside radius
                  _voxs.at(a, b, c) &= ~OVERHANG;
                }
              }
            }
          }
          _voxs.at(i, j, k) |= OVERHANG;
        }
      }
    }
  }
}

// --------------------------------------------------------------

template <bool sign> // true: positive, false: negative
class ProjX
{
public:
  inline v2u  projSize (const Array3D<uchar>& voxs) const { return v2u(voxs.ysize(), voxs.zsize());               }
  inline v3u  projLoop (const Array3D<uchar>& voxs) const { return v3u(voxs.ysize(), voxs.zsize(), voxs.xsize()); }
  inline v3u  projCoord(const v3u coord)            const { return v3u(coord[2], coord[0], coord[1]);             }
  inline uint initValue(const Array3D<uchar>& voxs) const { return (sign ? voxs.xsize() - 1 : 0);                 }
  inline bool getSign  ()                           const { return sign;                                          }
};

template <bool sign> // true: positive, false: negative
class ProjY
{
public:
  inline v2u  projSize (const Array3D<uchar>& voxs) const { return v2u(voxs.xsize(), voxs.zsize());               }
  inline v3u  projLoop (const Array3D<uchar>& voxs) const { return v3u(voxs.xsize(), voxs.zsize(), voxs.ysize()); }
  inline v3u  projCoord(const v3u coord           ) const { return v3u(coord[0], coord[2], coord[1]);             }
  inline uint initValue(const Array3D<uchar>& voxs) const { return (sign ? voxs.ysize() - 1 : 0);                 }
  inline bool getSign  ()                           const { return sign;                                          }
};

template <bool sign> // true: positive, false: negative
class ProjZ
{
public:
  inline v2u  projSize (const Array3D<uchar>& voxs) const { return v2u(voxs.xsize(), voxs.ysize());               }
  inline v3u  projLoop (const Array3D<uchar>& voxs) const { return v3u(voxs.xsize(), voxs.ysize(), voxs.zsize()); }
  inline v3u  projCoord(const v3u coord)            const { return v3u(coord[0], coord[1], coord[2]);             }
  inline uint initValue(const Array3D<uchar>& voxs) const { return (sign ? voxs.zsize() - 1 : 0);                 }
  inline bool getSign  ()                           const { return sign;                                          }
};

template<class P>
Array2D<uint> orthogonalProjection(const Array3D<uchar>& voxs, const v3i axis)
{
  const P projAxis;
  Array2D<uint> proj;
  v2u projSize = projAxis.projSize(voxs);
  proj.allocate(projSize[0], projSize[1]);
  proj.fill(projAxis.initValue(voxs));
  v3u loopSize = projAxis.projLoop(voxs);
  ForIndex(a, loopSize[0]) {
    ForIndex(b, loopSize[1]) {
      bool sign = projAxis.getSign();
      for (int c = (sign ? 0 : static_cast<int>(loopSize[2] - 1)); (sign ? c < static_cast<int>(loopSize[2]) : c >= 0); (sign ? c++ : c--)) {
        v3u coord = projAxis.projCoord(v3u(a, b, c));
        if (voxs.at(coord[0], coord[1], coord[2]) & (ALONG_X | ALONG_Y | ALONG_Z)) {
          proj.at(a, b) = c;
          break;
        }
      }
    }
  }
  return proj;
}

// --------------------------------------------------------------

void collisionFilter(Array3D<uchar>& _voxs)
{
  Array2D<uint> projX = orthogonalProjection<ProjX<true>>(_voxs, v3i(1, 0, 0));
  Array2D<uint> projY = orthogonalProjection<ProjY<true>>(_voxs, v3i(0, 1, 0));
  Array2D<uint> projZ = orthogonalProjection<ProjZ<true>>(_voxs, v3i(0, 0, 1));

  ForIndex(k, _voxs.zsize()) {
    ForIndex(j, _voxs.ysize()) {
      ForIndex(i, _voxs.xsize()) {
        if (_voxs.at(i, j, k) & OVERHANG) {
          if (projZ.at(i, j) < static_cast<uint>(k)) { // collides with shape
            ForRangeReverse(c, k - 1, 0) {             // all the way down to the bed or...
              if (_voxs.at(i, j, c) & (ALONG_X | ALONG_Y | ALONG_Z)) { // ...stop when coliding with shape
                break;
              }
              _voxs.at(i, j, c) |= OVERHANG;           // move this after the if below if anchors can connect to bridges (bars)
              if (projX.at(j, c) == static_cast<uint>(_voxs.xsize() - 1) || projY.at(i, c) == static_cast<uint>(_voxs.ysize() - 1)) { // visible from X+- or Y+-
                for (int z = c; z < c + COLLISION_RETRACT && z < k && z < static_cast<int>(_voxs.zsize()); z++) {
                  _voxs.at(i, j, z) &= ~OVERHANG;      // unmark retract length
                }
                break;
              }
            }
          }
        }
      }
    }
  }
}

// --------------------------------------------------------------


int main(int argc, char **argv)
{

  try {

    // load triangle mesh
    TriangleMesh_Ptr mesh(loadTriangleMesh(SRC_PATH "/model.stl"));
    // produce (fixed fp) integer vertices and triangles
    std::vector<v3i> pts;
    std::vector<v3u> tris;
    {
      float factor = 0.95f;
      m4x4f boxtrsf = scaleMatrix(BOX_SCALE)
        * scaleMatrix(v3f(1.f) / tupleMax(mesh->bbox().extent()))
        * translationMatrix((1 - factor) * 0.5f * mesh->bbox().extent())
        * scaleMatrix(v3f(factor))
        * translationMatrix(-mesh->bbox().minCorner());

      // transform vertices
      pts.resize(mesh->numVertices());
      ForIndex(p, mesh->numVertices()) {
        v3f pt   = mesh->posAt(p);
        v3f bxpt = boxtrsf.mulPoint(pt);
        v3i ipt  = v3i(clamp(round(bxpt), v3f(0.0f), BOX_SCALE - v3f(1.0f)));
        pts[p]   = ipt;
      }
      // prepare triangles
      tris.reserve(mesh->numTriangles());
      ForIndex(t, mesh->numTriangles()) {
        v3u tri = mesh->triangleAt(t);
        tris.push_back(tri);
      }
    }

    // rasterize into voxels
    v3u resolution(mesh->bbox().extent() / tupleMax(mesh->bbox().extent()) * static_cast<float>(VOXEL_RESOLUTION));
    Array3D<uchar> voxs(resolution);
    voxs.fill(0);
    {
      Timer tm("rasterization");
      Console::progressTextInit((int)tris.size());
      ForIndex(t, tris.size()) {
        Console::progressTextUpdate(t);
        rasterize<swizzle_xyz>(tris[t], pts, voxs); // xy view
        rasterize<swizzle_yzx>(tris[t], pts, voxs); // yz view
        rasterize<swizzle_zxy>(tris[t], pts, voxs); // zx view
      }
      Console::progressTextEnd();
      std::cerr << std::endl;
    }

    {
      Timer tm("filtering");
      std::cerr << "filtering parity/sphere/collision ... ";

      insideFilter(voxs);

#if FILTER_PARITY
      parityFilter(voxs);
#endif

#if FILTER_SPHERE
      sphereFilter(voxs);
#endif

#if FILTER_COLLISION
      collisionFilter(voxs);
#endif

      std::cerr << " done." << std::endl;
    }

    // add inner voxels
#if VOXEL_FILL_INSIDE
    {
      Timer tm("fill");
      std::cerr << "filling in/out ... ";
#if VOXEL_ROBUST_FILL
      fillInsideVoting(voxs);
#else
      fillInside(voxs); // winding order; determines if a voxel is solid or empty (detects holes inside stl)
#endif
      std::cerr << " done." << std::endl;
    }
#endif

    // save the result
    saveAsVox(SRC_PATH "/out.slab.vox", voxs);

    // report some stats
    int num_in_vox = 0;
    ForArray3D(voxs, i, j, k) {
      if (voxs.at(i, j, k) > 0) {
        num_in_vox++;
      }
    }
    std::cerr << "number of set voxels: " << num_in_vox << std::endl;

  } catch (Fatal& e) {
    std::cerr << "[ERROR] " << e.message() << std::endl;
  }

}

/* -------------------------------------------------------- */
