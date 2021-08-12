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
/// \file He3aEventAction.cc
/// \brief Implementation of the He3EventAction class

#include "He3EventAction.hh"
#include "He3RunAction.hh"
#include "He3Analysis.hh"

#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4UnitsTable.hh"

#include "Randomize.hh"
#include <iomanip>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

He3EventAction::He3EventAction()
 : G4UserEventAction(),
   fEnergyAbs(0.),
   fEnergyGap(0.),
   fTrackLAbs(0.),
   fTrackLGap(0.),
   fEnergyAbs1(0.),
   fEnergyAbs2(0.),
   fTrackLAbs1(0.),
   fTrackLAbs2(0.),
   fpositionHe3{0.},
   fpositionHe3c{0.},
   fpositionHe3s1{0.},
   fpositionHe3s2{0.}
   //fpositionX(0.),
   //fpositionY(0.),
   //fpositionZ(0.)

{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

He3EventAction::~He3EventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void He3EventAction::BeginOfEventAction(const G4Event* /*event*/)
{  
  // initialisation per event
  fEnergyAbs = 0.;
  fEnergyGap = 0.;
  fTrackLAbs = 0.;
  fTrackLGap = 0.;
  fEnergyAbs1 = 0.;
  fEnergyAbs2 = 0.;
  fTrackLAbs1 = 0.;
  fTrackLAbs2 = 0.;
  //fpositionX = 0.;
  //fpositionY = 0.;
  //fpositionZ = 0.;
  fpositionHe3[0] = 0.0;
  fpositionHe3[1] = 0.0;
  fpositionHe3[2] = 0.0;

  fpositionHe3c[0] = 0.0;
  fpositionHe3c[1] = 0.0;
  fpositionHe3c[2] = 0.0;
  //fpositionHe3[3] = 0.0;

  fpositionHe3s1[0] = 0.0;
  fpositionHe3s1[1] = 0.0;
  fpositionHe3s1[2] = 0.0;

  fpositionHe3s2[0] = 0.0;
  fpositionHe3s2[1] = 0.0;
  fpositionHe3s2[2] = 0.0;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void He3EventAction::EndOfEventAction(const G4Event* event)
{
  // Accumulate statistics
  //

  // get analysis manager
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  // fill histograms
  if (fEnergyAbs > 0) analysisManager->FillH1(0, fEnergyAbs);
  if (fEnergyGap > 0) analysisManager->FillH1(1, fEnergyGap);

  // Central Tube with cuts
  // x-axis +/- 4", y & z center

  if ( fpositionHe3[0] > (-4.0*25.4)  && fpositionHe3[0] < (4.0*25.4) &&
       //std::sqrt( ((fpositionHe3[1] - 0.0) * (fpositionHe3[1] - 0.0)) + (((fpositionHe3[2] - (29.16)) * (fpositionHe3[2] - (29.16)) ))) <= 11.8872
       std::sqrt( ((fpositionHe3[1] - 0.0) * (fpositionHe3[1] - 0.0)) + ((fpositionHe3[2] - 0.0) * (fpositionHe3[2] - 0.0)) ) <= 11.8872
     )
      {  //if (fEnergyAbs > 0)  analysisManager->FillH1(7, fEnergyAbs);
         if (fEnergyAbs > 0)  analysisManager->FillH1(2, fEnergyAbs); } 

     //filling 2D histogram --> Central tube

      if ( fpositionHe3[0] > (-4.0*25.4)  && fpositionHe3[0] < (4.0*25.4) &&
           std::sqrt( ((fpositionHe3[1] - 0.0) * (fpositionHe3[1] - 0.0)) + ((fpositionHe3[2] - 0.0) * (fpositionHe3[2] - 0.0)) ) <= 11.8872
         )
      {
        analysisManager->FillH2(0, fpositionHe3[0], fpositionHe3[1]); 
        analysisManager->FillH2(1, fpositionHe3[1], fpositionHe3[2]); 
        analysisManager->FillH2(2, fpositionHe3[2], fpositionHe3[0]);
      }  

//filling 2D histogram --> Side tube 1 (left)

  // fill ntuple
  analysisManager->FillNtupleDColumn(0, fEnergyAbs);
  analysisManager->FillNtupleDColumn(1, fEnergyGap);
  analysisManager->FillNtupleDColumn(2, fEnergyAbs1);

  // Central tubes with cuts 
  if ( fpositionHe3[0] > (-4.0*25.4)  && fpositionHe3[0] < (4.0*25.4) &&
       std::sqrt( ((fpositionHe3[1] - 0.0) * (fpositionHe3[1] - 0.0)) + ((fpositionHe3[2] - 0.0) * (fpositionHe3[2] - 0.0)) ) <= 11.8872
     )
      {  if (fEnergyAbs > 0)  analysisManager->FillNtupleDColumn(3, fEnergyAbs); } 

  analysisManager->FillNtupleDColumn(4, fpositionHe3[0]);
  analysisManager->FillNtupleDColumn(5, fpositionHe3[1]);
  analysisManager->FillNtupleDColumn(6, fpositionHe3[2]);

  analysisManager->AddNtupleRow();  
  
  // Print per event (modulo n)
  //
  auto eventID = event->GetEventID();
  auto printModulo = G4RunManager::GetRunManager()->GetPrintProgress();
  if ( ( printModulo > 0 ) && ( eventID % printModulo == 0 ) ) {
    //G4cout << "---> End of event: " << eventID << G4endl <<G4endl;     
  }
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
