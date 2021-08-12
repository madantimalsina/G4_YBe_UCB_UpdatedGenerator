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
// $Id: YBe_PrimaryGeneratorMessenger.cc $
// 
/// \file YBe_PrimaryGeneratorMessenger.cc
/// \brief Definition of the YBe_PrimaryGeneratorMessenger class

#include "YBe_PrimaryGeneratorMessenger.hh"

#include "YBe_PrimaryGeneratorAction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"

YBe_PrimaryGeneratorMessenger::YBe_PrimaryGeneratorMessenger(YBe_PrimaryGeneratorAction* generator)
  : G4UImessenger(),
    primaryGenerator(generator),
    primaryGeneratorDir(0),
    generatorModeCmd(0)
{
	primaryGeneratorDir = new G4UIdirectory("/generator/");
	primaryGeneratorDir->SetGuidance("Generator control.");

	generatorModeCmd = new G4UIcmdWithAString("/generator/setMode", this);
	generatorModeCmd->SetGuidance("Select gammas, neutrons");
	generatorModeCmd->SetParameterName("mode", true);
	generatorModeCmd->SetDefaultValue("neutrons");
	generatorModeCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

YBe_PrimaryGeneratorMessenger::~YBe_PrimaryGeneratorMessenger()
{
	delete generatorModeCmd;
	delete primaryGeneratorDir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void YBe_PrimaryGeneratorMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
	if( command == generatorModeCmd )
	{	
		if(newValue != "gammas" && newValue != "neutrons")
		{
			G4cout << "\n ERROR! Bad generator mode definition: " << newValue << G4endl;
			return;
		}
		G4cout << "Setting generator mode to " << newValue << G4endl;
		primaryGenerator->SetGeneratorMode(newValue);
	}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
