
OpenFOAM 2D Pillar Case (WSL / OF 2315)

Steps:

1. Place your Fusion-exported STL here:
   constant/triSurface/pillars.stl

2. Run:
   blockMesh
   snappyHexMesh -overwrite
   checkMesh
   pimpleFoam

Notes:
- 2D is enforced via frontAndBack = empty
- Ensure STL is watertight
- locationInMesh MUST be inside fluid region
