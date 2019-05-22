//G4
#include "G4Material.hh"
#include "G4VSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Point3D.hh"
#include "G4TwoVector.hh"
#include "globals.hh"

//My
#include "fTOFConst.hh"
#include "crtConst.hh"

struct VolumeStruct {
  G4Material*        material;
  G4VSolid*          solid;
  G4LogicalVolume*   logical;
  G4VPhysicalVolume* physical;
  VolumeStruct() :
    material(0),
    solid(0),
    logical(0),
    physical(0)
  {;}
  ~VolumeStruct() {;}
};

struct WorldStruct : VolumeStruct {
  // Defined as a G4Box
  const G4double sizeX;
  const G4double sizeY;
  const G4double sizeZ;
  WorldStruct() :
    sizeX(200.0*cm),
    sizeY(200.0*cm),
    sizeZ(600.0*cm)      
  {;}
};

struct SecStruct : VolumeStruct {
  // Defined as a G4Box 
  const G4double sizeX;
  const G4double sizeY;
  const G4double sizeZ;
  SecStruct() :
    sizeX(fTOFConst::barSizeX),
    sizeY(fTOFConst::barSizeY),
    sizeZ(fTOFConst::barSizeZ)
  {;}
};

struct PmtWindow1Struct : VolumeStruct {
  // Defined as a G4Box
  const G4double sizeX;
  const G4double sizeY;
  const G4double sizeZ;
  PmtWindow1Struct() :
    sizeX(fTOFConst::pmtWin1SizeX),
    sizeY(fTOFConst::pmtWin1SizeY),
    sizeZ(fTOFConst::pmtWin1SizeZ)
  {;}
};

struct PmtWindow2Struct : VolumeStruct {
  // Defined as a G4Box
  const G4double sizeX;
  const G4double sizeY;
  const G4double sizeZ;
  PmtWindow2Struct() :
    sizeX(fTOFConst::pmtWin2SizeX),
    sizeY(fTOFConst::pmtWin2SizeY),
    sizeZ(fTOFConst::pmtWin2SizeZ)
  {;}
};

struct SensitiveStruct : VolumeStruct {
  // Defined as a G4Box
  const G4double sizeX;
  const G4double sizeY;
  const G4double sizeZ;
  SensitiveStruct() :
    sizeX(fTOFConst::pmtChSizeX),
    sizeY(fTOFConst::pmtChSizeY),
    sizeZ(fTOFConst::pmtChSizeZ)
  {;}
};

struct PmtBoxStruct : VolumeStruct {
  // Defined as a G4Box
  const G4double sizeX;
  const G4double sizeY;
  const G4double sizeZ;
  const G4double sizeXgapp;
  PmtBoxStruct() :
    sizeX(fTOFConst::pmtBoxSizeX),
    sizeY(fTOFConst::pmtBoxSizeY),
    sizeZ(fTOFConst::pmtBoxSizeZ),
    sizeXgapp(fTOFConst::pmtBoxGapp)
  {;}
};

struct PmtAbs1Struct : VolumeStruct {
  // Defined as a G4Box
  const G4double sizeX;
  const G4double sizeY;
  const G4double sizeZ;
  PmtAbs1Struct() :
    sizeX(fTOFConst::pmtAbs1SizeX),
    sizeY(fTOFConst::pmtAbs1SizeY),
    sizeZ(fTOFConst::pmtAbs1SizeZ)
  {;}
};

struct Abs1Struct : VolumeStruct {
  // Defined as a G4Box
  const G4double sizeX;
  const G4double sizeY;
  const G4double sizeZ;
  Abs1Struct() :
    sizeX(fTOFConst::abs1SizeX),
    sizeY(fTOFConst::abs1SizeY),
    sizeZ(fTOFConst::abs1SizeZ)
  {;}
};

struct Abs2Struct : VolumeStruct {
  // Defined as a G4Box
  const G4double sizeX;
  const G4double sizeY;
  const G4double sizeZ;
  Abs2Struct() :
    sizeX(fTOFConst::abs2SizeX),
    sizeY(fTOFConst::abs2SizeY),
    sizeZ(fTOFConst::abs2SizeZ)
  {;}
};





///////////////////// fTOF bar 10.04.18 ///////////////

struct BarBoxStruct : VolumeStruct {
  // Defined as a G4Box
  const G4double sizeX;
  const G4double sizeY;
  const G4double sizeZ;
  BarBoxStruct() :
    sizeX(fTOFConst::barBoxSizeX),
    sizeY(fTOFConst::barBoxSizeY),
    sizeZ(fTOFConst::barBoxSizeZ)
  {;}
};

struct BarTrdXYStruct : VolumeStruct {
  // Defined as a G4Trd
  const G4double dX1;
  const G4double dX2;
  const G4double dY1;
  const G4double dY2;
  const G4double dZ;
  BarTrdXYStruct() :
    dX1(fTOFConst::barBoxSizeX/2.),
    dX2(fTOFConst::barBoxSizeX/2. - fTOFConst::barChamfer/2.),
    dY1(fTOFConst::barBoxSizeY/2.),
    dY2(fTOFConst::barBoxSizeY/2. - fTOFConst::barChamfer/2.),
    dZ(fTOFConst::barChamfer/4.)
  {;}
};

struct BarTrdYZStruct : VolumeStruct {
  // Defined as a G4Trd
  const G4double dX1;
  const G4double dX2;
  const G4double dY1;
  const G4double dY2;
  const G4double dZ;
  BarTrdYZStruct() :
    dX1(fTOFConst::barBoxSizeZ/2.),
    dX2(fTOFConst::barBoxSizeZ/2. - fTOFConst::barChamfer/2.),
    dY1(fTOFConst::barBoxSizeY/2.),
    dY2(fTOFConst::barBoxSizeY/2. - fTOFConst::barChamfer/2.),
    dZ(fTOFConst::barChamfer/4.)
  {;}
};

struct BarTrdZXStruct : VolumeStruct {
  // Defined as a G4Trd
  const G4double dX1;
  const G4double dX2;
  const G4double dY1;
  const G4double dY2;
  const G4double dZ;
  BarTrdZXStruct() :
    dX1(fTOFConst::barBoxSizeX/2.),
    dX2(fTOFConst::barBoxSizeX/2. - fTOFConst::barChamfer/2.),
    dY1(fTOFConst::barBoxSizeZ/2.),
    dY2(fTOFConst::barBoxSizeZ/2. - fTOFConst::barChamfer/2.),
    dZ(fTOFConst::barChamfer/4.)
  {;}
};


struct HamPmtWindowStruct : VolumeStruct { //////////////////// 11.04.18
  // Defined as a G4Box
  const G4double sizeX;
  const G4double sizeY;
  const G4double sizeZ;
  HamPmtWindowStruct() :
    sizeX(fTOFConst::hamWindowThick),
    sizeY(fTOFConst::hamSensitiveSize),
    sizeZ(fTOFConst::hamSensitiveSize)
  {;}
};

struct PlanPmtWindowStruct : VolumeStruct { //////////////////// 11.04.18
  // Defined as a G4Box
  const G4double sizeX;
  const G4double sizeY;
  const G4double sizeZ;
  PlanPmtWindowStruct() :
    sizeX(fTOFConst::planWindowThick),
    sizeY(fTOFConst::planSensitiveSize),
    sizeZ(fTOFConst::planSensitiveSize)
  {;}
};

struct PlanPmtChannelStruct : VolumeStruct { //////////////////// 11.04.18
  // Defined as a G4Box
  const G4double sizeX;
  const G4double sizeY;
  const G4double sizeZ;
  PlanPmtChannelStruct() :
    sizeX(fTOFConst::planChanThick),
    sizeY(fTOFConst::planChanSize),
    sizeZ(fTOFConst::planChanSize)
  {;}
};

struct HamPmtChannelStruct : VolumeStruct { //////////////////// 11.04.18
  // Defined as a G4Box
  const G4double sizeX;
  const G4double sizeY;
  const G4double sizeZ;
  HamPmtChannelStruct() :
    sizeX(fTOFConst::hamChanThick),
    sizeY(fTOFConst::barBoxSizeY),
    sizeZ(fTOFConst::barBoxSizeZ)  //hamChanSize
  {;}
};

struct HamPmtBoxStruct : VolumeStruct { //////////////////// 11.04.18
  // Defined as a G4Box
  const G4double sizeX;
  const G4double sizeY;
  const G4double sizeZ;
  HamPmtBoxStruct() :
    sizeX(fTOFConst::hamBoxDepth),
    sizeY(fTOFConst::hamBoxSize),
    sizeZ(fTOFConst::hamBoxSize)
  {;}
};

struct PlanPmtBoxStruct : VolumeStruct { //////////////////// 11.04.18
  // Defined as a G4Box
  const G4double sizeX;
  const G4double sizeY;
  const G4double sizeZ;
  PlanPmtBoxStruct() :
    sizeX(fTOFConst::planBoxDepth),
    sizeY(fTOFConst::planBoxSize),
    sizeZ(fTOFConst::planBoxSize)
  {;}
};


//////////////////////////////////////////////////////////////////
