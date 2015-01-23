#include <bits/hhc-phasespace.h>
#include <bits/hhc-jetfunc.h>

// --- fastNLO v2.2: define user class to be used with NLOJet++
class FastNLOUserHHC : public nlo::basic_user_set<nlo::user0d_hhc, nlo::user1h_hhc, nlo::user2h_hhc> {
public:
   FastNLOUserHHC();
   // --- fastNLO user: evaluate steering file and define physics output (called once before first event)
   virtual void phys_output(const std::basic_string<char>& fname, unsigned long nsave = 10000UL, bool txt = false);
   // --- fastNLO user: analyze parton event (called once for each event)
   virtual void userfunc(const event_type&, const nlo::amplitude_hhc&);
   // --- fastNLO v2.2: count events and store table (called after each event)
   virtual void end_of_event();
   // --- fastNLO v2.2: read settings from steering file
   virtual void read_steering() = 0;

   // --- fastNLO v2.2: define global pointer to fastNLO steering file
   struct fastNLOCreate *ftable;

   static FastNLOUserHHC *instance;

private:
   // --- fastNLO definitions (not for user)
   double nevents;            // No. of events calculated so far
   unsigned long nwrite;      // Actual no. of events after which to write out the table
   unsigned long nwritemax;   // Maximal no. of events after which to write out the table
};
