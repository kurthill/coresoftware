#ifndef PHG4TrackUserInfoV1_H__
#define PHG4TrackUserInfoV1_H__

#include <Geant4/G4VUserTrackInformation.hh>

class PHG4Shower;

// Made this with "V1" in the name in case we ever want to inherit from
// it with other versions...

// Use the UserTrackInformation to attach a flag telling the framework 
// to save the track in the truth output.  Other uses might include keeping
// track of the 

class PHG4TrackUserInfoV1 : public G4VUserTrackInformation
{
public:
  PHG4TrackUserInfoV1() : G4VUserTrackInformation("TrackUserInfoV1"),
			  usertrackid(0), userparentid(0), userprimaryid(0),
			  wanted(0), keep(0) , shower(NULL) {}
  virtual ~PHG4TrackUserInfoV1() {}

  void Print() const {
    G4cout << "PHG4TrackUserInfoV1: " << std::endl;
    G4cout << "   UserTrackId = " << usertrackid << std::endl;
    G4cout << "   UserParentId = " << userparentid << std::endl;
    G4cout << "   UserPrimaryId = " << userprimaryid << std::endl;
    G4cout << "   Wanted = " << wanted << std::endl;
    G4cout << "   Keep = " << keep << std::endl;
  }

  void SetUserTrackId(const int val) {usertrackid = val;}
  int GetUserTrackId() const {return usertrackid;}
  
  void SetUserParentId(const int val) {userparentid = val;}
  int GetUserParentId() const {return userparentid;}

  void SetUserPrimaryId(const int val) {userprimaryid = val;}
  int GetUserPrimaryId() const {return userprimaryid;}

  void SetWanted(const int val) {wanted = val;}
  int GetWanted() const {return wanted;}

  void SetKeep(const int val) {keep = val;}
  int GetKeep() const {return keep;}

  void SetShower(PHG4Shower *ptr) {shower = ptr;}
  PHG4Shower* GetShower() const {return shower;}

private:
  int usertrackid;
  int userparentid;
  int userprimaryid;
  int wanted;
  int keep;
  PHG4Shower* shower;
};

// Utility function to wrap up the operations involved with adding user info
// to the track.
class G4Track;

namespace PHG4TrackUserInfo 
{
  void SetUserTrackId(G4Track* track, const int usertrackid);
  void SetUserParentId(G4Track* track, const int userparentid);
  void SetUserPrimaryId(G4Track* track, const int userprimaryid);
  void SetWanted(G4Track* track, const int wanted);
  void SetKeep(G4Track* track, const int keep);
  void SetShower(G4Track* track, PHG4Shower* ptr);
};

#endif
