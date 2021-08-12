//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// 
/// \file He3DetectorConstruction.cc
/// \brief Implementation of the He3DetectorConstruction class

#include "He3DetectorConstruction.hh"

#include "G4Material.hh"
#include "G4RunManager.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Polycone.hh"
#include "G4Tubs.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"

#include "G4PVReplica.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "G4SubtractionSolid.hh"
#include "G4BooleanSolid.hh"
#include "G4VSolid.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"
#include "G4Transform3D.hh"
#include "G4AffineTransform.hh"

 #include "G4MultiUnion.hh"



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreadLocal 
G4GlobalMagFieldMessenger* He3DetectorConstruction::fMagFieldMessenger = nullptr; 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

He3DetectorConstruction::He3DetectorConstruction()
 : G4VUserDetectorConstruction(),
   fAbsorberPV(nullptr),
   fGapPV(nullptr),
   fAbsorberPV1(nullptr),
   fAbsorberPV2(nullptr),
   fCheckOverlaps(true)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

He3DetectorConstruction::~He3DetectorConstruction()
{ 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* He3DetectorConstruction::Construct()
{
  // Define materials 
  DefineMaterials();
  
  // Define volumes
  return DefineVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void He3DetectorConstruction::DefineMaterials()
{ 
  // Define materials

  G4double A;  // atomic mass (mass of a mole)
  G4double Z;  // atomic number (mean number of protons)
  G4double d;  // density
  G4double fractionMass; // Fraction by Mass (Weight %)
  G4double abundance;
  // Lead material defined using NIST Manager
  auto nistManager = G4NistManager::Instance();
  nistManager->FindOrBuildMaterial("G4_Pb");
  //nistManager->FindOrBuildMaterial("G4_TEFLON");
  nistManager->FindOrBuildMaterial("G4_POLYETHYLENE");
  nistManager->FindOrBuildMaterial("G4_POLYPROPYLENE"); 
  nistManager->FindOrBuildMaterial("G4_AIR");
  nistManager->FindOrBuildMaterial("G4_Al");
  nistManager->FindOrBuildMaterial("G4_W");

//*******************************************************************

  //General Material defination
  G4Element* elH = new G4Element ("Hydrogen","H",Z = 1.,A = 1.01*g/mole);
  G4Element* elBe = new G4Element("Beryllium","Be", Z=4, A=9.1218*g/mole);
  G4Element* elC = new G4Element("Carbon","C",Z = 6.,A = 12.011*g/mole);
  G4Element* elO = new G4Element("Oxygen","O",Z = 8.,A = 16.00*g/mole);
  G4Element* elCr  = new G4Element("Chromium","Cr",Z = 24.,A = 52.00*g/mole);
  G4Element* elFe = new G4Element("Iron","Fe", Z=26., A = 55.845*g/mole);
  G4Element* elNi = new G4Element("Nickel","Ni", Z=28., A = 58.69*g/mole);
  G4Element* elW = new G4Element("Tungsten","W",  Z=74., A = 183.84*g/mole);
  G4Element* elAl  = new G4Element("Aluminum","Al", Z = 13.0, A = 26.982*g/mole);
  G4Element* elMg  = new G4Element("Magnesium","Mg", Z = 12.0, A = 24.305*g/mole);
  G4Element* elS  = new G4Element("Sulfur","S", Z = 16.0, A = 32.06*g/mole);

 //Molybdenum definition
    G4Isotope* Mo92 = new G4Isotope("Mo92", 42, 92, 91.906809 * g / mole);
    G4Isotope* Mo94 = new G4Isotope("Mo94", 42, 94, 93.9050853 * g / mole);
    G4Isotope* Mo95 = new G4Isotope("Mo95", 42, 95, 94.9058411 * g / mole);
    G4Isotope* Mo96 = new G4Isotope("Mo96", 42, 96, 95.9046785 * g / mole);
    G4Isotope* Mo97 = new G4Isotope("Mo97", 42, 97, 96.9060205 * g / mole);
    G4Isotope* Mo98 = new G4Isotope("Mo98", 42, 98, 97.9054073 * g / mole);
    G4Isotope* Mo100 = new G4Isotope("Mo100", 42, 100, 99.907477 * g / mole);

  G4Element* elMo = new G4Element("Natural Mo", "elMo", 7);
    elMo->AddIsotope(Mo92, abundance=14.84*perCent);
    elMo->AddIsotope(Mo94, abundance=9.25*perCent);
    elMo->AddIsotope(Mo95, abundance=15.92*perCent);
    elMo->AddIsotope(Mo96, abundance=16.68*perCent);
    elMo->AddIsotope(Mo97, abundance=9.55*perCent);
    elMo->AddIsotope(Mo98, abundance=24.13*perCent);
    elMo->AddIsotope(Mo100, abundance=9.63*perCent);

  //YBe Source Assembly Materials

  //MT-185 alloys of Tungsten used to make YBe Pig (LZYBePig)
  // Midwest Tungsten Service -> 97.0% (W), 2.1% (Ni), 0.9% (Fe)   
  // Mi Tech Tungsten Metals -> 97.05% (W), 2.04% (Ni), 0.91% (Fe)  
  d = 18.52*g/cm3;
  G4Material* YBeTungsten_MT185 = new G4Material("Tungsten_MT185",d,3);
  YBeTungsten_MT185->AddElement(elW, fractionMass=97.05*perCent);
  YBeTungsten_MT185->AddElement(elNi, fractionMass=2.04*perCent);
  YBeTungsten_MT185->AddElement(elFe, fractionMass=0.91*perCent);
   //YBeTungsten_MT185->AddElement(elW,fractionMass=0.97);
   //YBeTungsten_MT185->AddElement(elNi,fractionMass=0.021);
   //YBeTungsten_MT185->AddElement(elFe,fractionMass=0.009);

/*
  // LZYBeSource (Beryllium metal)
  d = 1.85*g/cm3;
  G4Material* Beryllium = new G4Material("Beryllium", d, 1);
  Beryllium->AddElement(elBe, fractionMass=100.0*perCent);
*/
  // Beryllium metal
  G4Material* Beryllium = new G4Material("Beryllium", d = 1.85*g/cm3, 7);
  Beryllium->AddElement(elBe, fractionMass=99.0892*perCent);
  Beryllium->AddElement(elO,  fractionMass= 0.6378*perCent);
  Beryllium->AddElement(elC,  fractionMass= 0.06*perCent);
  Beryllium->AddElement(elFe, fractionMass= 0.11*perCent);
  Beryllium->AddElement(elAl, fractionMass= 0.05*perCent);
  Beryllium->AddElement(elMg, fractionMass= 0.02*perCent);
  Beryllium->AddElement(elS,  fractionMass= 0.03*perCent);


  // LZYBeDisk (Made up of SS-316)
  d = 7.99*g/cm3;
  G4Material* SS316 = new G4Material("SS316",d,4);
  SS316->AddElement(elFe, fractionMass=68.5*perCent);
  SS316->AddElement(elCr, fractionMass=17.0*perCent);
  SS316->AddElement(elNi, fractionMass=12.0*perCent);
  SS316->AddElement(elMo, fractionMass=2.5*perCent);
  

//************************************************************************
  // plexiglass, lucite 
  d = 1.19*g/cm3;
  G4Material* matplexiglass = new G4Material("Plexiglass",d,3);
  matplexiglass->AddElement(elH, fractionMass=8.0*perCent);
  matplexiglass->AddElement(elC, fractionMass=60.0*perCent);
  matplexiglass->AddElement(elO, fractionMass=32.0*perCent);
//*****************************************************************************

// He3 gas for UCB/LBNL He3 tube

  //Volume of the He3 tubes
  G4double volume_He3 = (8.0*2.54*pi*1.18872*1.18872);

  // Build He3 gas
  //(From Junsong's email) LND replies 
  // The gas mix for the LND 252 is 91.2 Torr CO2, 2948.8 Torr He-3. 
  // Total fill pressure is 3040 Torr
  G4int protons=2, neutrons=1, nucleons=protons+neutrons;
  G4double elements;
  G4double atomicMass_He3 = 3.016*g/mole; //molar mass
  G4Isotope* isoHe3 = new G4Isotope("He3", protons, nucleons, atomicMass_He3);
  G4Element* elHe3 = new G4Element("Helium3", "He3", 1);
  elHe3->AddIsotope(isoHe3, 100*perCent);
  G4double pressure_He3 = 2948.8/760.0*atmosphere; // 2948.8 Torr He-3
  G4double temperature = 293.15*kelvin;
  G4double molar_constant = Avogadro*k_Boltzmann;
  //G4double density = (atomicMass*pressure)/(temperature*molar_constant);
  // PV = nRT --> PV= (m/M)RT --> P/RT = D/M --> D = MP/RT
  G4double density_He3 = (atomicMass_He3*pressure_He3)/(temperature*molar_constant);

  G4Material* Helium3 = new G4Material("Helium3", density_He3, elements=1, kStateGas, temperature, pressure_He3);
  Helium3->AddElement(elHe3, fractionMass=100*perCent);

  // carbon dioxide (Co2) gas in non STP conditions
  G4int ncomponents, natoms;
  G4double atomicMass_Co2 = 44.01*g/mole; //molar mass (12.01 + 2*16.0)
  G4double pressure_Co2 = (91.2/760.0)*atmosphere; //91.2 Torr CO2
  G4double temperature_Co2 = 293.15*kelvin;
  // PV = nRT --> PV= (m/M)RT --> P/RT = D/M --> D = MP/RT
  G4double density_Co2 = (atomicMass_Co2*pressure_Co2)/(temperature_Co2*molar_constant); 
  G4Material* CO2 = new G4Material("carbon dioxide", density_Co2, ncomponents=2, kStateGas, temperature_Co2, pressure_Co2);
  CO2->AddElement(elC, natoms=1);
  CO2->AddElement(elO, natoms=2);

  // Now mixture He3 + Co2
  //Amount of He3 with respect to total  = 2948.8/3040*100 =  97.0%
  //Amount of Co2 with respect to total  = 91.2/3040*100 =  3.0%
  //G4double atomicMass_He3Co2 = ((2948.8/3040.0)*3.016 + (91.2/3040.0)*44.01)*g/mole;
  G4double atomicMass_He3Co2 = ((0.97*3.016) + (0.03*44.01))*g/mole;
  G4double pressure_He3Co2 = (3040.0/760.0)*atmosphere;
  G4double density_He3Co2 = (atomicMass_He3Co2*pressure_He3Co2)/(temperature*molar_constant);
  G4Material* He3Co2 = new G4Material("He3Co2"  , density_He3Co2, 2,  kStateGas, temperature, pressure_He3Co2);
  He3Co2->AddMaterial( Helium3, (0.97 *3.016) / ((0.97 *3.016) + (0.03* 44.01)) );
  He3Co2->AddMaterial( CO2, (0.03* 44.01) / ((0.97 *3.016) + (0.03* 44.01)) );

/*
  //Volume of the He3 tubes
  G4double volume_He3 = (13*2.54*pi*1.18872*1.18872);

  // Build He3 gas
  G4int protons=2, neutrons=1, nucleons=protons+neutrons;
  G4double elements;
  G4double atomicMass_He3 = 3.016*g/mole; //molar mass
  G4Isotope* isoHe3 = new G4Isotope("He3", protons, nucleons, atomicMass_He3);
  G4Element* elHe3 = new G4Element("Helium3", "He3", 1);
  elHe3->AddIsotope(isoHe3, 100*perCent);
  G4double pressure_He3 = 11.02*atmosphere;
  G4double temperature = 293.15*kelvin;
  G4double molar_constant = Avogadro*k_Boltzmann;
  //G4double density = (atomicMass*pressure)/(temperature*molar_constant);
  G4double density_He3 = (atomicMass_He3*pressure_He3)/(temperature*molar_constant);

  G4Material* Helium3 = new G4Material("Helium3", density_He3, elements=1, kStateGas, temperature, pressure_He3);
  Helium3->AddElement(elHe3, fractionMass=100*perCent);

  // Argon 
  G4double atomicMass_Ar = 39.948*g/mole;
  G4double pressure_Ar = 0.58*atmosphere;
  G4double density_Ar = (atomicMass_Ar*pressure_Ar)/(temperature*molar_constant);
  G4Element* elAr = new G4Element("Argon", "Ar", Z=18., atomicMass_Ar);
  G4Material* Argon = new G4Material("Argon"  , density_Ar, 1,  kStateGas, temperature, pressure_Ar);
  Argon->AddElement(elAr, 1);

  // 95% He3 + 4.95% Ar + 0.05% CH4,
  G4double atomicMass_He3Ar = ((0.95 *3.016) + (0.05* 39.948))*g/mole;
  G4double pressure_He3Ar = 11.60*atmosphere;
  G4double density_He3Ar = (atomicMass_He3Ar*pressure_He3Ar)/(temperature*molar_constant);

  G4Material* He3Ar = new G4Material("He3Ar"  , density_He3Ar, 2,  kStateGas, temperature, pressure_He3Ar);
  He3Ar->AddMaterial( Helium3, (0.95 *3.016) / ((0.95 *3.016) + (0.05* 39.948)) );
  He3Ar->AddMaterial( Argon, (0.05* 39.948) / ((0.95 *3.016) + (0.05* 39.948)) );

*/
  // UHMW (Ultra High Molecular Weight Polyethylene)
  d = 0.93*g/cm3;
  nistManager->BuildMaterialWithNewDensity("UHMWPE","G4_POLYETHYLENE",d);


  // UHMW for Thermal Scattering of neutron
  // For containt (http://www.apc.univ-paris7.fr/~franco/g4doxy/html/G4NistMaterialBuilder_8cc-source.html)
  G4Element *H = new G4Element("TS_H_of_Polyethylene", "H", 1., 1.0079*g/mole);
  G4Material *TS_of_Polyethylene = new G4Material("TS_of_Polyethylene", 0.94*g/cm3, 2, kStateSolid, 293.15*kelvin);
  TS_of_Polyethylene->AddElement(H, fractionMass=14.3711*perCent);
  TS_of_Polyethylene->AddElement(elC, fractionMass=85.6289*perCent);


  // POLYPROPYLENE
  d = 913.43685*mg/cm3; // Geant4 density: 900.000 mg/cm3
  nistManager->BuildMaterialWithNewDensity("POLYPROPYLENE","G4_POLYPROPYLENE",d);

  // POLYPROPYLENE for Thermal Scattering of neutron
  // For containt (http://www.apc.univ-paris7.fr/~franco/g4doxy/html/G4NistMaterialBuilder_8cc-source.html)
  G4Material *TS_of_Polypropylene = new G4Material("TS_of_Polypropylene", 913.43685*mg/cm3, 2, kStateSolid, 293.15*kelvin);
  TS_of_Polypropylene->AddElement(H, fractionMass=14.3711*perCent);
  TS_of_Polypropylene->AddElement(elC, fractionMass=85.6289*perCent);


  G4cout<< G4endl;
  G4cout << "**************************************************************************"<< G4endl;
  G4cout << "**************************************************************************"<< G4endl;
  G4cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<G4endl; 
  G4cout << "MY CHECK....!!!!  " << "pressure_He3  =  " << pressure_He3/atmosphere << " atmosphere"<< G4endl;
  G4cout << "MY CHECK....!!!!  " << "pressure_Co2  =  " << pressure_Co2/atmosphere<< " atmosphere" << G4endl;
  G4cout << "MY CHECK....!!!!  " << "pressure_He3Co2  =  " << pressure_He3Co2/atmosphere<< " atmosphere" << G4endl;
  G4cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<G4endl;
  G4cout << "MY CHECK....!!!!  " << "density_He3  =  " << density_He3/(mg/cm3)<< " kg/m3" << G4endl;
  G4cout << "MY CHECK....!!!!  " << "density_Co2  =  " << density_Co2/(mg/cm3)<< " kg/m3" << G4endl;
  G4cout << "MY CHECK....!!!!  " << "density_He3Co2  =  " << density_He3Co2/(mg/cm3)<< " kg/m3" << G4endl;
  G4cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<G4endl; 
  G4cout << "MY CHECK....!!!!  " << "volume of the tube  =  " << volume_He3 << " cm3"<< G4endl;
  G4cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<G4endl;
  G4cout << "MY CHECK....!!!!  " << "Mass_He3  =  " << (density_He3/(mg/cm3))*volume_He3 << " mg" << G4endl;
  G4cout << "MY CHECK....!!!!  " << "Mass_Co2  =  " << (density_Co2/(mg/cm3))*volume_He3 << " mg"<<G4endl;
  G4cout << "MY CHECK....!!!!  " << "mass_He3Co2  =  " << (density_He3Co2/(mg/cm3))*volume_He3 << " mg"<< G4endl;
  G4cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<G4endl;
  G4cout << "**************************************************************************"<< G4endl;
  G4cout << "**************************************************************************"<< G4endl;

  // Print materials
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* He3DetectorConstruction::DefineVolumes()
{
  // Geometry parameters

  G4double startAngle = 0.*deg;
  G4double spanningAngle = 360.*deg;

  // Outer layer of He3 tube
  G4double outerRadius_He3Tubes = 25.4/2*mm;
  G4double thickness_AlTube = (0.032*25.4)*mm;


  // Parameters: Inner layer of He3
  // Tubes filled with He3 gas + ...
  G4double outerRadius_He3Gass = outerRadius_He3Tubes - thickness_AlTube;
  G4double inerrRadius_He3Gass = 0.*mm;

  // This is the effective length (13") of He-3 tubes
  // divide by 2 because of to match G4 simulation style (-33.02/2 to +33.02/2)
  G4double zHalfHeight_He3Tube = (11.20*2.54)/2*cm; 
  G4double zHalfHeight_He3Gass = (8.0*2.54)/2*cm; 
  
  // Rotation of the Tubes 
  // (both outer and inner filled with He3 gas)
  G4RotationMatrix* rotD1 = new G4RotationMatrix();
  G4RotationMatrix* rotD2 = new G4RotationMatrix();
  G4RotationMatrix* rotD3 = new G4RotationMatrix();
  rotD1->rotateZ(90.*deg);
  rotD2->rotateX(90.*deg);
  rotD3->rotateY(90.*deg);

  //G4double moderatorThickness = 2.5*2.540*cm;  // Thickness of moderator (Varies: 1", 2", 2.5", 5" etc)
  
  // (From Junsong email) The round-tube-shape UHMW PE are all 12 inch long. 
  //I attach the Xometry order for the UHMW  PE tubes. 
  //There is a 1-inch dimaeter, 1 inch long UHMW PE plug on the far side of the He-3 tube.
  G4double innerRadius_PE = (1.02/2.0)*25.4*mm; // 1.02" internal diameter
  G4double outerRadius_PE = 1.50*25.4*mm; // (UCB measurement for 0.5", 1.0", 1.5", 2.0", 2.5", 3" thickness (radius))
  G4double zHalfHeight_PE = 6.0*25.4*mm;
  // 1" PE plug
  G4double innerRadius_PE2 = 0.0*mm;
  G4double outerRadius_PE2 = 0.5*25.4*mm;
  G4double zHalfHeight_PE2 = 0.5*25.4*mm;

  //G4double boxSizeX  = 30.48*cm;  // 12"
  G4double boxSizeX  = 31.115*cm;  // 12.25"
  //G4double boxSizeY  = 30.48*cm;   // 12" 
  G4double boxSizeXY = 100.0*cm;
  auto boxSizeZ = 100.0*cm;
  auto worldSizeXY = 5.0 * boxSizeX; 
  auto worldSizeZ  = 5.0 * boxSizeZ;
  
  // Parameters: source holder YBe source
  G4double innerRadius_sourceHolder = 0.*cm;
  G4double outerRadius_sourceHolder = 10.*cm;
  G4double zHalfHeight_YBeSourceHolder = (20./2.)*cm;
  G4double outerRadius_YBeSourceHolderCap = (10.52/2.0)*cm;
  G4double zHalfHeight_YBeSourceHolderCap = (4./2.)*cm;
  //G4double Position_YBeSourceHolder = moderatorThickness/2 + zHalfHeight_YBeSourceHolder; // just add thickness here for table height etc.
  // (From Junsong email) The He-3 tube was horizontally orientated, and its center was 4.875" (12.3825*cm) below the bottom of the tungsten. 
  //The middle of the 12-inch long UHMW PE tube is center to the Tungsten Shielding. 
  //We keep the He-3 tube positions relative to the tungsten shielding the same in the different scenarios of PE thicknesses.
  G4double Position_YBeSourceHolder = 12.3825*cm + zHalfHeight_YBeSourceHolder; // just add thickness here for table height etc.
  G4double Position_YBeSourceHolderCap = Position_YBeSourceHolder+ zHalfHeight_YBeSourceHolder + zHalfHeight_YBeSourceHolderCap;
  G4double outerRadius_YBeBeO = (2.54/2.)*cm;
  //G4double zHalfHeight_YBeBeO = (6.46/2.)*cm; //My old value
  G4double zHalfHeight_YBeBeO = (6.858/2.)*cm;  // new value from Andy
  //Bottom of BeO start at 12.31646 (4.894 in CAd Drawing) cm from base of YBe Pig
  //Centre of Be at total 15.74546 (6.199") from bottom of the W pig
  G4double Position_YBeBeO =  Position_YBeSourceHolder + 2.31646*cm + (6.858/2.)*cm; //Centre of Be at total 15.74546 (6.199") from bottom of the W pig
  G4double outerRadius_YDisk = (2.35)*mm;
  G4double zHalfHeight_YDisk = (4.6/2.)*mm;
  //From CAD drawing distance from bottom of the pig to bottom of Y88 disk 6.349" 
  G4double Position_YDisk = Position_YBeSourceHolder + 6.12646*cm; //directly from CAD drawing
  //G4double Position_YDisk = Position_YBeBeO + 0.635*cm;

  // Get materials
  auto defaultMaterial = G4Material::GetMaterial("G4_AIR");
  //auto absorberMaterial = G4Material::GetMaterial("Helium3");
  auto absorberMaterial = G4Material::GetMaterial("He3Co2");
  //auto absorberMaterial = G4Material::GetMaterial("He3Ar");
  //auto moderatorMaterial = G4Material::GetMaterial("G4_POLYETHYLENE");
  //auto moderatorMaterial = G4Material::GetMaterial("UHMWPE");
  auto moderatorMaterial = G4Material::GetMaterial("TS_of_Polyethylene");
  //auto moderatorMaterial2 = G4Material::GetMaterial("POLYPROPYLENE");
  //auto moderatorMaterial2 = G4Material::GetMaterial("TS_of_Polypropylene");
  //auto holdMaterial = G4Material::GetMaterial("Plexiglass");
  auto moderatorMaterialT = G4Material::GetMaterial("G4_Al");
  //auto holdMaterial_Tungsten = G4Material::GetMaterial("G4_W");
  auto holdMaterial_Tungsten_MT185 = G4Material::GetMaterial("Tungsten_MT185");
  auto BeMaterial = G4Material::GetMaterial("Beryllium");
  auto YMaterial = G4Material::GetMaterial("SS316");

  
  if ( ! defaultMaterial || ! absorberMaterial || ! moderatorMaterial ) {
    G4ExceptionDescription msg;
    msg << "Cannot retrieve materials already defined."; 
    G4Exception("He3DetectorConstruction::DefineVolumes()",
      "MyCode0001", FatalException, msg);
  }  
  // Volume defination:
  //     
  // World
  //
  auto worldS 
    = new G4Box("World",           // its name
                 worldSizeXY/2, worldSizeXY/2, worldSizeZ/2); // its size
                         
  auto worldLV
    = new G4LogicalVolume(
                 worldS,           // its solid
                 defaultMaterial,  // its material
                 "World");         // its name
                                   
  auto worldPV
    = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(),  // at (0,0,0)
                 worldLV,          // its logical volume                         
                 "World",          // its name
                 0,                // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 
   
   //                               
  // Layer (This is a Box to put every little geometry on it !!!)
  //
  auto layerS 
    = new G4Box("Layer",           // its name
                 boxSizeXY/2, boxSizeXY/2, boxSizeZ/2.); // its size
                         
  auto layerLV
    = new G4LogicalVolume(
                 layerS,           // its solid
                 defaultMaterial,  // its material
                 "Layer");         // its name

/*  new G4PVReplica(
                 "Layer",          // its name
                 layerLV,          // its logical volume
                 worldLV,          // its mother
                 kZAxis,           // axis of replication
                 1,        // number of replica
                 boxSizeZ);  // witdth of replica*/

    new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(),  // at (0,0,0)
                 layerLV,          // its logical volume                         
                 "Layer",    // its name
                 worldLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 

  //***************************************************************************************************//                               
  //
  // Source holder for YBe neutron source
  //YBe source Pig Cap 
  //             
  G4Tubs* sourceTubeYBeCap = 
    new G4Tubs("YBeSourceCap", innerRadius_sourceHolder, outerRadius_YBeSourceHolderCap, zHalfHeight_YBeSourceHolderCap, startAngle, spanningAngle);

  G4LogicalVolume* YBeCap =                         
    new G4LogicalVolume(sourceTubeYBeCap,            //its solid
                        holdMaterial_Tungsten_MT185, //its material
                        "YBeSourceCap"); 

  fHoldPV
   = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0., -Position_YBeSourceHolderCap), // its position (moderatorThickness/2 + zHalfHeight_YBeSourceHolder)
                 YBeCap,            // its logical volume                         
                 "YBeSourceCap",            // its name
                 layerLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 
  //
  //YBe Pig
  //               
  G4Tubs* sourceTube = 
    new G4Tubs("YBeSourcePig", innerRadius_sourceHolder, outerRadius_sourceHolder, zHalfHeight_YBeSourceHolder, startAngle, spanningAngle);
      /*   
     // substracting the voulume if we do not on different logic volume!!!
     // outerRadius_YBeBeO+0.3937*mm, 0.7874/2 added because of 1.03 inch of source inner can diameter 
     // zHalfHeight_YBeBeO+2.63*mm, 5.25/2 added because of 2.75 inch of source inner can height 

      G4Tubs* sourceTube_BeOVol = new G4Tubs("Source_BeoVol", innerRadius_sourceHolder, outerRadius_YBeBeO+0.3937*mm, zHalfHeight_YBeBeO+2.63*mm, startAngle, spanningAngle); 
      // substraction of Beo overlap volume from Tungsten Pig
      G4VSolid* subtract_BeoVol = new G4SubtractionSolid("SourcePig-Source_BeoVol", sourceTube, sourceTube_BeOVol, 0, G4ThreeVector(0., 0., -5.18*cm));
      G4LogicalVolume* YBePig = new G4LogicalVolume(subtract_BeoVol,holdMaterial_Tungsten_MT185,"SourcePig"); 
      G4LogicalVolume* YBePig = new G4LogicalVolume(sourceTube, holdMaterial_Tungsten_MT185,"SourcePig"); */
  
   G4LogicalVolume* YBePig = new G4LogicalVolume(sourceTube,holdMaterial_Tungsten_MT185,"YBeSourcePig"); 
  
  fHoldPV
   = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0., -Position_YBeSourceHolder), // its position (moderatorThickness/2 + zHalfHeight_YBeSourceHolder)
                 YBePig,            // its logical volume                         
                 "YBeSourcePig",            // its name
                 layerLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 

  //
  // Source BeO
 /*  // YBe source position is based on the Andy's BACCARAT simulation
   // https://luxzeplin.gitlab.io/docs/softwaredocs/simulations/baccarat/event_generator/generators_details/YBe.html
   //3.23*cm --> z-offset of the center of the BeO volume relative to the center of the total tungsten volume (main cylinder + cap) 
   // 5.23*cm --> (3.23 + 2.0)*cm, I added 2.0*cmm because Andy's values based on the (main cylinder + cap) 
   // which is opposite then I assumed but in both cases base of the Ybe source start at 12 cm from bottom of the YBe Pig
   //So distance from the bottom of the tungston pig to center of Beo volume is 12 + 3.23 = 15.23 cm
  */
  // New math based on the CAD drawing
  //(and also match with Andy new simulation https://gitlab.com/biekerta/photoneutron-safety-sims/-/blob/main/YBe_pig/src/YBeSafeDetectorConstruction.cc)
  //Bottom of BeO start at 12.31646 (4.894 in CAd Drawing) cm from base of YBe Pig
  // length of Be metal 2.7" --> 6.858 cm 
  //Centre of Be at total 12.31646 + (6.858/2) = 15.74546 (6.199") from bottom of the W pig
  // base of Be start at 12.31646 (4.894 in CAd Drawing) cm from base of YBe Pig 
  // so z -position for z with respect to ventre of the 20 cm height w pig is 2.31646 cm +  6.858/2 cm (length of Be metal 2.7") = 5.74546
  //  
  G4Tubs* sourceTubeBeO = 
    new G4Tubs("YBeSource_Be", innerRadius_sourceHolder, outerRadius_YBeBeO, zHalfHeight_YBeBeO, startAngle, spanningAngle);

/*  
  // substracting the voulume if we do not on different logic volume!!!
  G4Tubs* sourceTube_YDiskVol = new G4Tubs("Source_YDiskVol", innerRadius_sourceHolder, outerRadius_YDisk+0.01*mm, zHalfHeight_YDisk+0.04*mm, startAngle, spanningAngle);
  // substraction of overlap Y-88 Disk volume from BeO
  G4VSolid* subtract_YDiskVol = new G4SubtractionSolid("SourceBeo-Source_YDiskVol", sourceTubeBeO, sourceTube_YDiskVol, 0, G4ThreeVector(0., 0., -1.23357*cm));
  G4LogicalVolume* BeVolume = new G4LogicalVolume(subtract_YDiskVol,BeMaterial,"SourceBeo"); 
    */
  G4LogicalVolume* BeVolume = new G4LogicalVolume(sourceTubeBeO,BeMaterial,"YBeSource_Be"); 
  
  fHoldPV
   = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0., -(5.74546*cm)), // its position (moderatorThickness/2 + zHalfHeight_YBeSourceHolder)
                 BeVolume,            // its logical volume                         
                 "YBeSource_Be",            // its name
                 YBePig,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 


  //
  // source YDisk
  // New math based on the CAD drawing 
  //(and also match with Andy new simulation https://gitlab.com/biekerta/photoneutron-safety-sims/-/blob/main/YBe_pig/src/YBeSafeDetectorConstruction.cc)
  //Bottom of Be start at 12.31646 (4.894 in CAd Drawing) cm from base of W Pig
  //Centre of Be at total 15.74546 (6.199") from bottom of the W pig
  // bottom of the Y88 start at 16.12646 cm (6.349") from the bottom of the W pig
  // and centre at 16.12646 cm + zHalfHeight_YDisk (4.6 mm/2.) = 16.35646
  // so the off set is 16.35646-15.74546 = 0.611 cm (check with Andy code and he has same value)
     G4Tubs* sourceTubeYDisk = 
    new G4Tubs("YBeSource_YDisk", innerRadius_sourceHolder, outerRadius_YDisk, zHalfHeight_YDisk, startAngle, spanningAngle);

  G4LogicalVolume* YDisk =                         
    new G4LogicalVolume(sourceTubeYDisk,           //its solid
                        YMaterial,                 //its material
                        "YBeSource_YDisk"); 

  fHoldPV
   = new G4PVPlacement(
                 0,                // no rotation
                 //G4ThreeVector(0., 0., -(1.23*cm)),
                 G4ThreeVector(0., 0., (-0.611*cm)),
                 YDisk,            // its logical volume                         
                 "YBeSource_YDisk",            // its name
                 BeVolume,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 
/*
    // if we substract the voulume for the different logic volume!!! here e.g layerLV
     fHoldPV
   = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0., -(Position_YDisk)), // its position (moderatorThickness/2 + zHalfHeight_YBeSourceHolder)
                 YDisk,            // its logical volume                         
                 "SourceYDisk",            // its name
                 layerLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps */

//***************************************************************************************************//                               
  // Neutron Moderator 
  // Here neutron moderator material is UHMWPE (Ultra High Molecular Weight Polyethylene)
  //
  //G4double innerRadius_PE = 25.6*mm;
  //G4double outerRadius_PE = 3*25.4*mm;
  //G4double zHalfHeight_PE = 6*25.4*mm;
  auto gapS = new G4Tubs("Moderator", innerRadius_PE, outerRadius_PE, zHalfHeight_PE, startAngle, spanningAngle);
  //auto gapS = new G4Box("Moderator", boxSizeX/2, boxSizeY/2, moderatorThickness/2);                     
  auto gapLV = new G4LogicalVolume(gapS, moderatorMaterial,"Moderator");
                                   
  fGapPV
    = new G4PVPlacement(
                 rotD3,                // no rotation
                 G4ThreeVector(0., 0., 0.), // its position
                 gapLV,            // its logical volume                         
                 "Moderator",            // its name
                 layerLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 
  
  // For 1 inch PE cap at far end
  auto gapS2 = new G4Tubs("Moderator2", innerRadius_PE2, outerRadius_PE2, zHalfHeight_PE2, startAngle, spanningAngle);
  //auto gapS = new G4Box("Moderator", boxSizeX/2, boxSizeY/2, moderatorThickness/2);                     
  auto gapLV2 = new G4LogicalVolume(gapS2, moderatorMaterial,"Moderator2");
                                   
  fGapPV
    = new G4PVPlacement(
                 rotD3,                // no rotation
                 G4ThreeVector(-5.5*25.4, 0., 0.), // its position
                 gapLV2,            // its logical volume                         
                 "PE_Plug_1in",            // its name
                 layerLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 
  //
  // For He3 tube outerlayer
  // Tube 1 (center tube) outerlayer
  // centre tube is off by 0.635cm(0.25")
  // Length of the tube (not symmetric)
  // 1.6"-1.25" = 0.35" active length extending out of moderator at NEAR END
  // 1.6"-0.95" = 0.65" active length extending out of moderator at FAR END
  // So the X-center is pushed by the 0.15" 
  //

  //G4double zPosition_He3Tubes = (moderatorThickness/2.) + outerRadius_He3Tubes + 3.6*mm; 
  //G4double zPosition_He3Gass = (moderatorThickness/2.) + outerRadius_He3Tubes + 3.6*mm;
 
  // Almunium Outer layer 
  auto gapsT1 
    = new G4Tubs("AlmuniumOuterLayer_centralTube", outerRadius_He3Gass, outerRadius_He3Tubes, zHalfHeight_He3Tube, startAngle, spanningAngle); 
                         
  auto gapLVT1
    = new G4LogicalVolume(
                 gapsT1,        // its solid
                 moderatorMaterialT, // its material
                 "AlmuniumOuterLayer_centralTube");          // its name
                                   
  fGapPV
    = new G4PVPlacement(
                 rotD3,                // no rotation
                 G4ThreeVector(0.6*25.4, 0., 0.),
                 //G4ThreeVector((0.15 * 2.54)*cm, 0, zPosition_He3Tubes),
                 gapLVT1,       // its logical volume                         
                 "AlmuniumOuterLayer_centralTube",           // its name
                 layerLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 

  //                               
  // Absorber
  // Absorber is 95% He3 + 4.95% Ar + 0.05% CH4,
  // For He3 gas Tube 1 (Center tube) 
  // Centre tube is off by 0.635cm (0.25")
  // Length of the tube (not symmetric)
  // 1.6"-1.25" = 0.35" active length extending out of moderator at NEAR END
  // 1.6"-0.95" = 0.65" active length extending out of moderator at FAR END
  // So the X-center is pushed by the 0.15"  
  auto absorberS 
    = new G4Tubs("Absorber_centralTube", inerrRadius_He3Gass, outerRadius_He3Gass, zHalfHeight_He3Gass, startAngle, spanningAngle); 
                         
  auto absorberLV
    = new G4LogicalVolume(
                 absorberS,        // its solid
                 absorberMaterial, // its material
                 "Absorber_centralTube");          // its name
                                   
  fAbsorberPV
    = new G4PVPlacement(
                 rotD3,                // no rotation
                 G4ThreeVector(0., 0., 0.),
                 //G4ThreeVector((0.15 * 2.54)*cm, 0, zPosition_He3Gass), // its position Centre tube is off by 0.635cm (0.25")
                 absorberLV,       // its logical volume                         
                 "Absorber_centralTube",           // its name
                 layerLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 

  // print parameters
  ///////
  G4cout<< G4endl;
  G4cout << "**************************************************************************"<< G4endl;
  G4cout << "**************************************************************************"<< G4endl;
  G4cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<G4endl;
  G4cout << "MY CHECK....!!!!  " << "Gap Thickness (Moderator Radius)  =  " << outerRadius_PE << " mm"<< G4endl;
  G4cout << "MY CHECK....!!!!  " << "Position of Y-88 source disk  =  " << Position_YDisk << " mm"<< G4endl;
  G4cout << "MY CHECK....!!!!  " << "Position of center of BeO volume  =  " << Position_YBeBeO << " mm"<< G4endl;
  G4cout << "MY CHECK....!!!!  " << "Position of YBe Source Holder  =  " <<  Position_YBeSourceHolder << " mm"<< G4endl;
  G4cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<G4endl; 
  //G4cout << "MY CHECK....!!!!  " << "X-position of all He3 tubes (const. /w M.T.)  =  " << 0.15*2.54 << " mm "<< G4endl;
  //G4cout << "MY CHECK....!!!!  " << "Y-position of 2 side He3 tubes (const. /w M.T.)  =  +/-" << 5.0*25.4 << " mm "<< G4endl;
  //G4cout << "MY CHECK....!!!!  " << "Y-position of central He3 tube (const. /w M.T.)  =  " << 0.0*25.4 << " mm "<< G4endl;
  //G4cout << "MY CHECK....!!!!  " << "Z-position of all He3 tubes (change /w M.T.)  =  " << zPosition_He3Tubes << " mm "<< G4endl;
  //G4cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<G4endl;
  G4cout << "**************************************************************************"<< G4endl;
  G4cout << "**************************************************************************"<< G4endl;
  //                                        
  // Visualization attributes
  //
  G4VisAttributes* blue_clear = new G4VisAttributes(G4Colour(0.0, 0.0, 1.0, 0.4));
  G4VisAttributes* aqua_clear = new G4VisAttributes(G4Colour(0.0, 0.5, 1.0, 0.7));
  G4VisAttributes* brown_clear = new G4VisAttributes(G4Colour(0.45, 0.25, 0.0, 0.7));
  G4VisAttributes* yellow_clear = new G4VisAttributes(G4Colour(1.0, 1.0, 0.0, 0.5));
  G4VisAttributes* red_clear = new G4VisAttributes(G4Colour(1.0, 0.0, 0.0, 0.5));
  G4VisAttributes* red_clear2 = new G4VisAttributes(G4Colour(1.0, 0.0, 0.0, 0.2));
  
  worldLV->SetVisAttributes (G4VisAttributes::GetInvisible());
  layerLV->SetVisAttributes (G4VisAttributes::GetInvisible());
  YBePig->SetVisAttributes(brown_clear);
  BeVolume->SetVisAttributes(yellow_clear);
  YDisk->SetVisAttributes(red_clear2);
  //YDisk->SetVisAttributes(G4VisAttributes(G4Colour::Red()));
  YBeCap->SetVisAttributes(brown_clear);

  //absorberLV->SetVisAttributes(G4VisAttributes(G4Colour::Red()));
  absorberLV->SetVisAttributes(red_clear);
  //absorberLV2->SetVisAttributes(G4VisAttributes(G4Colour::Red()));
  //absorberLV3->SetVisAttributes(G4VisAttributes(G4Colour::Red()));

  gapLVT1->SetVisAttributes(blue_clear);
  //gapLVT2->SetVisAttributes(blue_clear);
  //gapLVT3->SetVisAttributes(blue_clear);

  gapLV->SetVisAttributes(aqua_clear);
  //gapLV->SetVisAttributes(G4VisAttributes(G4Colour::Green()));

  //auto simpleBoxVisAtt= new G4VisAttributes(G4Colour(1.0,2.0,3.0));
  //simpleBoxVisAtt->SetVisibility(true);
  //calorLV->SetVisAttributes(simpleBoxVisAtt);

  //
  // Always return the physical World
  //
  return worldPV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void He3DetectorConstruction::ConstructSDandField()
{ 
  // Create global magnetic field messenger.
  // Uniform magnetic field is then created automatically if
  // the field value is not zero.
  G4ThreeVector fieldValue;
  fMagFieldMessenger = new G4GlobalMagFieldMessenger(fieldValue);
  fMagFieldMessenger->SetVerboseLevel(1);
  
  // Register the field messenger for deleting
  G4AutoDelete::Register(fMagFieldMessenger);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
