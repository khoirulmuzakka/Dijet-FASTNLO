#ifndef __fastNLOUserBlock__
#define __fastNLOUserBlock__

#include <fstream>
#include <iostream>
#include <istream>
#include <string>
#include "speaker.h"

using namespace std;

class fastNLOUserBlock : public PrimalScream {

public:
   fastNLOUserBlock();
   ~fastNLOUserBlock();

   // i/o
   int Read(istream *table);					//! read UserBlock
   virtual void Write(ostream *table);				//! write to disk
   virtual void Print() const;

   int GetUserFlag() const {return fUserFlag;}			//! Get UserFlag of UserBlock	
   void SetUserFlag(int id) {fId = id;}				//! Set UserFlag of UserBlock
   
   int GetNLines() const { return fNLines;}			//! Get number of lines in table on disk
   void SetNLines(int nl) { fNLines = nl; fContent.resize(nl);}	//! Set number of lines. Should be consistent with content of user block. Must always be greater than 3

   vector<string> GetDescription() const { return fContent; }	//! Get description of user block
   void SetDescription(vector<string> description) { fDescr = Description;} //! Set description for user block

protected:
   int fNLines;							//! number of lines in user block (here: identical to fContent.size() when read)
   int fUserFlag;						//! UserFlag according to table definition
   vector<string> fDescr;					//! Description of user block

private:
   vector<string> fContent;					//! Plain-text content if read from disk

};

#endif
