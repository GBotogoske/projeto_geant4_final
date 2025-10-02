#ifndef RunAction_h
#define RunAction_h 1

#include "G4UserRunAction.hh"
#include "G4Accumulable.hh"
#include "G4AccumulableManager.hh"
#include "globals.hh"

class G4Run;
class RunActionMessenger;


class RunAction : public G4UserRunAction
{
  public:
    RunAction();
    virtual ~RunAction();

    // virtual G4Run* GenerateRun();
    virtual void BeginOfRunAction(const G4Run*);
    virtual void   EndOfRunAction(const G4Run*);
   
    void SetNameOfOutputFile(G4String name){fOutputFileName=name;};
    G4String GetNameOfOutputFile() const { return fOutputFileName;};

    void  GetInstance();
    
    void AddPenCount(G4int n)
    {
        fPenCount += n;
    };

  private:
    RunActionMessenger *fRunMessenger;
    G4String fOutputFileName;

    G4Accumulable<G4int> fPenCount;
};

#endif
