#ifndef STEPPINGACTION_H
#define STEPPINGACTION_H

#include "G4UserSteppingAction.hh"
#include "globals.hh"

G4ThreadLocal extern int cont_PEN;

class SteppingAction : public G4UserSteppingAction 
{
public:
    SteppingAction();
    virtual ~SteppingAction();
    virtual void UserSteppingAction(const G4Step* step) override;
    
};

#endif
