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
// $Id: YBe_PrimaryGeneratorAction.cc $
//
/// \file YBe_PrimaryGeneratorAction.cc
/// \brief Implementation of the YBe_PrimaryGeneratorAction class

#include "YBe_PrimaryGeneratorAction.hh"
#include "YBe_PrimaryGeneratorMessenger.hh"

#include "G4RunManager.hh"
#include "G4Navigator.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4VPhysicalVolume.hh"
#include "G4SolidStore.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4ThreeVector.hh"
#include "G4RandomDirection.hh"
#include "G4Neutron.hh"
#include "G4TransportationManager.hh"
#include "G4Navigator.hh"
#include "G4GenericIon.hh"
#include "G4IonTable.hh"
#include "G4VisExtent.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

YBe_PrimaryGeneratorAction::YBe_PrimaryGeneratorAction()
 : G4VUserPrimaryGeneratorAction(),
   fParticleSource(0),
   gamma_energy(1836.1*keV),
   generator_mode("undefined"),
   first_event(true),
   disk_global_position(0., 0., 0.),
   disk_PV(0),
   source_PV(0),
   disk_xmin(0),
   disk_xmax(0),
   disk_ymin(0),
   disk_ymax(0),
   disk_zmin(0),
   disk_zmax(0),
   boundingSphereRadius(0),
   process_threshold(-1664.54*keV),
   target_mass(7.456894 * GeV),
   particleDefinition(0),
   event_position(0., 0., 0.),
   accepted(false),
   gamma_direction(0., 0., 0.),
   gamma_position(0., 0., 0.),
   neutron_direction(0., 0., 0.),
   neutron_angle(0),
   neutron_energy(0),
   gamma_count(0),
   primaryGeneratorMessenger(0)

{
  fParticleSource = new G4ParticleGun();

  gNavigator = G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();

  primaryGeneratorMessenger = new YBe_PrimaryGeneratorMessenger(this);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

YBe_PrimaryGeneratorAction::~YBe_PrimaryGeneratorAction()
{
  delete fParticleSource;
  delete primaryGeneratorMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void YBe_PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  // Do the setup stuff here because detector construction initialized after
  // primary generator action
  if(first_event){
    G4VPhysicalVolume* source_disk_PV = G4PhysicalVolumeStore::GetInstance()->GetVolume("YBeSource_YDisk");
    G4VPhysicalVolume* current_PV = source_disk_PV;

    // here we get the global coordinates of the disk by looping up through mother volumes
    while( current_PV->GetName() != "World" ){
      disk_global_position += current_PV->GetTranslation();
      // This is really ugly... only one name per PV and LV and they must match!!
      current_PV = G4PhysicalVolumeStore::GetInstance()->GetVolume(current_PV->GetMotherLogical()->GetName());
    }

    gNavigator->LocateGlobalPointAndSetup(disk_global_position);
    G4TouchableHistoryHandle gTouchable = gNavigator->CreateTouchableHistoryHandle();

    // This gets the parent volume of the source disk (BeO)
    disk_PV = gTouchable->GetSolid(0);
    source_PV = gTouchable->GetSolid(1);

    // The extent is a bounding rectangular solid of the source volume
    // Use this to generate gamma points
    G4VisExtent disk_extent = disk_PV->GetExtent();
    disk_xmin = disk_extent.GetXmin();
    disk_xmax = disk_extent.GetXmax();
    disk_ymin = disk_extent.GetYmin();
    disk_ymax = disk_extent.GetYmax();
    disk_zmin = disk_extent.GetZmin();
    disk_zmax = disk_extent.GetZmax();

    // Use this to generate neutron points
    G4VisExtent source_extent = source_PV->GetExtent();
    G4double xmin = source_extent.GetXmin();
    G4double xmax = source_extent.GetXmax();
    G4double ymin = source_extent.GetYmin();
    G4double ymax = source_extent.GetYmax();
    G4double zmin = source_extent.GetZmin();
    G4double zmax = source_extent.GetZmax();

    // Use maximum ray in source as radius of bounding sphere to ensure total envelopment
    boundingSphereRadius = sqrt((xmax-xmin)*(xmax-xmin) + (ymax-ymin)*(ymax-ymin) + (zmax-zmin)*(zmax-zmin));

    first_event = false;
  }

  if(generator_mode == "gammas"){
    particleDefinition = G4ParticleTable::GetParticleTable()->GetIonTable()->GetIon(39, 88, 0);

    event_position = disk_global_position; 

    fParticleSource->SetParticlePosition(event_position);
    fParticleSource->SetParticleDefinition(particleDefinition);
    fParticleSource->SetParticleEnergy(0.*eV);

    fParticleSource->SetParticleMomentumDirection( G4RandomDirection() );
  } else if (generator_mode == "neutrons" ){
    particleDefinition = G4Neutron::Definition();

    accepted = false;
    gamma_direction = G4ThreeVector(0., 0., 0.);

    while( !accepted ){
      gamma_position = disk_global_position; 

      accepted = GenerateNeutronPoint(boundingSphereRadius, gamma_position,
                                      gNavigator, source_PV, gamma_direction, event_position);
    }

    neutron_direction = G4RandomDirection();
    neutron_angle = gamma_direction.angle(neutron_direction);

    neutron_energy = PhotoneutronEnergy(process_threshold, neutron_angle,
                                        gamma_energy, target_mass);

    fParticleSource->SetParticleDefinition(particleDefinition);

    fParticleSource->SetParticlePosition(event_position);
    fParticleSource->SetParticleEnergy(neutron_energy);
    fParticleSource->SetParticleMomentumDirection(neutron_direction);
  } else {
    G4cout << "ERROR! Bad generator mode definition: " << generator_mode << G4endl;
    return;
  }

  fParticleSource->GeneratePrimaryVertex(anEvent);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void YBe_PrimaryGeneratorAction::SetGeneratorMode(G4String mode)
{
  generator_mode = mode;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreeVector YBe_PrimaryGeneratorAction::GenerateGammaPoint(G4double xmax,
                G4double xmin, G4double ymax, G4double ymin, G4double zmax, G4double zmin,
                G4ThreeVector sourcePosition, G4Navigator* navigator, G4VSolid* diskVolume)
{
  G4double rand_x = (xmax - xmin) * G4UniformRand() + xmin;
  G4double rand_y = (ymax - ymin) * G4UniformRand() + ymin;
  G4double rand_z = (zmax - zmin) * G4UniformRand() + zmin;

  G4ThreeVector decayPosition = G4ThreeVector(rand_x, rand_y, rand_z) + sourcePosition;

  navigator->LocateGlobalPointAndSetup(decayPosition);

  while( navigator->CreateTouchableHistoryHandle()->GetSolid() != diskVolume){
    
    rand_x = (xmax - xmin) * G4UniformRand() + xmin;
    rand_y = (ymax - ymin) * G4UniformRand() + ymin;
    rand_z = (zmax - zmin) * G4UniformRand() + zmin;

    decayPosition = G4ThreeVector(rand_x, rand_y, rand_z) + sourcePosition;

    navigator->LocateGlobalPointAndSetup(decayPosition);

  }

  return decayPosition;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool YBe_PrimaryGeneratorAction::GenerateNeutronPoint(G4double sphereRadius,
                G4ThreeVector sourcePosition, G4Navigator* navigator, G4VSolid* parentVolume,
                G4ThreeVector& gammaDirection, G4ThreeVector& interactionPoint)
{
  G4double interactionLength = G4UniformRand() * sphereRadius;

  gammaDirection = G4RandomDirection();

  interactionPoint = sourcePosition + gammaDirection * interactionLength;

  navigator->LocateGlobalPointAndSetup(interactionPoint);

  return navigator->CreateTouchableHistoryHandle()->GetSolid() == parentVolume;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double YBe_PrimaryGeneratorAction::PhotoneutronEnergy(G4double threshold,
                G4double scatteringAngle, G4double gammaEnergy, G4double nucleusMass)
{
  G4double neutronMass = G4Neutron::Definition()->GetPDGMass();
  G4double firstTerm = nucleusMass * (gammaEnergy + threshold) / (neutronMass + nucleusMass);
  G4double secondTerm = gammaEnergy * sqrt((2 * neutronMass * nucleusMass) *
       (neutronMass + nucleusMass) * (gammaEnergy + threshold)) * cos(scatteringAngle) /
      ((neutronMass + nucleusMass) * (neutronMass + nucleusMass));

  return firstTerm + secondTerm;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
