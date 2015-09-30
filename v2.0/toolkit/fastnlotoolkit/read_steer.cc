// Author: Daniel Britzger
// DESY, 17/07/2012
//
//**********************************************************************************
//
//     D. Britzger
//     daniel.britzger@desy.de
//
//
//**********************************************************************************

#include "fastnlotk/read_steer.h"
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
#include <vector>

using namespace std;

map<string,read_steer*>* read_steer::instances = NULL;
//const string& read_steer::stdID = *(new string("SingleFileMode"));
string read_steer::stdID = "SingleFileMode";
int read_steer::fVerbosity = 2;

read_steer::read_steer() :
   str_sep(" \t")     , str_cmt("#"),
   str_arrbeg("{")    , str_arrend("}"),
   str_tabbeg("{{")   , str_tabend("}}"),
   str_nmspcbeg("{{{"), str_nmspcend("}}}"),
   //   str_inc("#include:"), fParseIncMode(0),
   str_inc(">>"), fParseIncMode(0),
   oW(" # WARNING! [read_steer] "), oI(" # INFO.    [read_steer] "), oE(" # ERROR!   [read_steer] ") {
}

read_steer* read_steer::Steering(string steerID) {
   if (instances==NULL) instances = new std::map<string,read_steer*>();
   // get singleton class
   if (!(*instances)[steerID]) { //new instance
      (*instances)[steerID] = new read_steer();
      if (steerID.compare(read_steer::stdID)!=0) {
         if (fVerbosity > 2) {
            cout << " # INFO.    [read_steer] " << "Initialized new read_steer namespace with steerID = '" << steerID << "'." << endl;
         }
      }
   }
   return (read_steer*)((*instances)[steerID]);
}


void read_steer::inits(string filename) {
   //    if ( ffilename != "" )
   //       cout << oW << "Filename already set (old=" << ffilename << ", new=" << filename << "). Is the used steerID unique?" << endl;
   if (filename == "")
      if (fVerbosity > 1) {cout << oW << "No filename specified." << endl;}
   if (ffilename !="") ffilename+=", ";
   ffilename += filename;
   fcurrentfilename = filename;
   read_stdin(fcurrentfilename);
}


void read_steer::destroy() {
   for (map<string, read_steer*>::iterator ii=(*instances).begin(); ii!=(*instances).end(); ++ii) {
      if ((*ii).second) {
         delete(*ii).second;
      }
      (*instances).erase((*ii).first);
   }
}

void read_steer::printall() {
   const string linesep = " +----------------------------------------------------------------------------+\n";
   const string l = " | ";
   cout << linesep;
   cout << l << "    read_steer. Printing all steering information.                         |" << endl;
   cout << linesep;
   for (map<string, read_steer*>::iterator ii=(*instances).begin(); ii!=(*instances).end(); ++ii) {
      cout << l << endl;
      cout << l << "steerID = '" << (*ii).first << "'" << endl;
      cout << linesep;
      (*ii).second->prt();
      cout << linesep;
   }
}


void read_steer::print(string steerID) {
   const string linesep = " +----------------------------------------------------------------------------+\n";
   const string l = " | ";
   cout << linesep;
   cout << l << "    read_steer. Printing steering information of steerID = '" << steerID << "'" << endl;
   cout << linesep;
   read_steer::Steering(steerID)->prt();
   cout << linesep;
}

void read_steer::prt() {
   const string l = " | ";
   // filename
   printf("%s%-30s%s\n",l.c_str(),"Filename(s)",ffilename.c_str());
   cout << l << endl;;
   //single values
   cout << l << "Single values" << endl;
   for (map<string,string>::iterator ii=fstrings.begin(); ii!=fstrings.end(); ++ii)
      if ((*ii).first!="")
         printf("%s   %-27s\t\t%s\n",l.c_str(),(*ii).first.c_str(),(*ii).second.c_str());
   // arrays
   if (!ffields.empty()) {
      cout << l << endl;
      cout << l << "Arrays" << endl;
      for (map<string,vector<string> >::iterator ii=ffields.begin(); ii!=ffields.end(); ++ii) {
         cout << l << "   " << (*ii).first <<  " {" << endl;
         for (unsigned int j = 0 ; j<(*ii).second.size() ; j++)
            cout  << l << "     [" << j << "]\t" << (*ii).second[j] << endl;
         cout  << l << "   }" << endl;
      }
   }
   // tables
   if (!ftables.empty()) {
      cout << l << endl;
      cout << l << "Tables/Matrices" << endl;
      for (map<string,vector<string> >::iterator ii=ftableheaders.begin(); ii!=ftableheaders.end(); ++ii) {
         cout << l << "   " << (*ii).first <<  " {{" << endl;
         cout << l << "     [H]\t";
         for (unsigned int j = 0 ; j<(*ii).second.size() ; j++)
            printf("%-10s",(*ii).second[j].c_str());
         cout <<  endl;
         vector<vector<string> > tab = ftables[(*ii).first];
         for (unsigned int ll = 0 ; ll<tab.size() ; ll++) {
            cout << l << "     [" << ll << "]\t";
            for (unsigned int j = 0 ; j<tab[ll].size() ; j++)
               printf("%-10s",tab[ll][j].c_str());
            cout << endl;
         }
         cout  << l << "   }}" << endl;
      }
   }
}

int read_steer::initnmspc(ifstream& strm, string filename) {
   if (ffilename !="") ffilename+=", ";
   ffilename += filename;
   //ffilename = filename;
   return readstrm(strm);
}

int read_steer::readstrm(ifstream& strm,unsigned int lstart, unsigned int lend, bool incfile) {
   if (!strm) {
      if (fVerbosity > 0) {cerr << oE << "This is not a valid stream." << endl;}
      return EXIT_FAILURE;
   }
   string lineread;
   if (!incfile) {
      fParseFieldMode = false;
      fParseTableMode = 0;
      ffieldlabel = "";
   }
   fParseIncMode++;
   unsigned int nlines=0;
   unsigned int rlines=0;
   while (std::getline(strm, lineread)) {
      nlines++;
      if (nlines < lstart) continue;
      if (lend!=0 && nlines > lend) break;
      rlines++;
      bool goon = ParseString(lineread);
      if (!goon) break;
   }
   fParseIncMode--;
   return rlines; // return nlines including comments
}

int read_steer::read_stdin(const string& filename) {
   //If the steering has alread been read -> do nothing
   ffile.open(filename.c_str());
   if (!ffile) {
      if (fVerbosity > 0) {cerr << oE << " Could not open file ('" << filename << "')." << endl;}
      return EXIT_FAILURE;
   }
   int n = readstrm(ffile);
   ffile.close();
   return n;
};


vector<bool> read_steer::getbf(const string& label) {
   vector<string> sf = ffields[label];
   vector<bool> ret(sf.size());
   for (unsigned int i = 0 ; i<sf.size() ; i++)
      ret[i] =  StringToBool(sf[i],label);
   return ret;
}

vector<int> read_steer::getif(const string& label) {
   vector<int> ret;
   vector<string> sf = ffields[label];
   for (unsigned int i = 0 ; i<sf.size() ; i++) {
      string val = sf[i];
      bool isnan = CheckInt(val.c_str());
      if (!isnan)
         if (fVerbosity > 1) {cout << oW << "Value number " <<  i << " of label='" << label << "' does not seem to be an integer number. value=" << val << endl;}
      ret.push_back(atoi(val.c_str()));
   }
   return ret;
}

vector<double> read_steer::getdf(const string& label) {
   vector<double> ret;
   vector<string> sf = ffields[label];
   for (unsigned int i = 0 ; i<sf.size() ; i++) {
      string val = sf[i];
      bool isnan = CheckNumber(val.c_str());
      if (!isnan)
         if (fVerbosity > 1) {cout << oW << "Value number " <<  i << " of label='" << label << "' does not seem to be a numeric number. value=" << val << endl;}
      ret.push_back(atof(val.c_str()));
   }
   return ret;
}

vector<string> read_steer::getsf(const string& label) {
   vector<string> ret = ffields[label];
   if (ret.empty())
      if (fVerbosity > 1) {cout  <<  oW << "Label '" <<   label  << "' was not found in list or has no values." <<  endl;}
   return ret;
}

vector<string> read_steer::getstcol(const string& label,const string& col) {
   // get column of a table with header 'col' as string values
   vector<string> ret;
   vector<string> head = getsthead(label);
   vector<vector<string> > tab = getst(label);
   for (vector<string>::size_type i = 0; i != head.size(); i++) {
      if (col.compare(head[i])==0) {
         bool IsErr = false;
         for (vector<string>::size_type j = 0; j != tab.size(); j++) {
            if ( i >= tab[j].size() ) {
               if (fVerbosity > 1) {cout << oW << "Row "<<j<<" has too little columns ("<<tab[j].size()<<"). Column '"<<col<<"' cannot be read correctly."<<endl;}
               IsErr = true;
            }
            else
               ret.push_back(tab[j][i]);
         }
         if ( IsErr )  ret.resize(tab.size(),"read_steer::warning");
         return ret;
      }
   }
   if (fVerbosity > 1) {cout << oW << "Column '" << col << "' was not found in table '" << label << "'." << endl;}
   return ret;
}


vector<bool> read_steer::getbtcol(const string& label,const string& col) {
   // get column of a table with header 'col' as string values
   vector<string> scol = getstcol(label,col);
   vector<bool> ret(scol.size());
   for (vector<string>::size_type i = 0; i != scol.size(); i++)
      ret[i] = StringToBool(scol[i]);
   return ret;
}


vector<int> read_steer::getitcol(const string& label,const string& col) {
   // get column of a table with header 'col' as string values
   vector<int> ret;
   vector<string> scol = getstcol(label,col);
   for (vector<string>::size_type i = 0; i != scol.size(); i++) {
      string val = scol[i];
      if (!CheckInt(val.c_str()))
         if (fVerbosity > 1) {
            cout << oW << "Value number " << i << " of table='" << label
                 << "' in column '" << col << "' does not seem to be an integer number. value=" << val << endl;
         }
      ret.push_back(atoi(val.c_str()));
   }
   return ret;
}

vector<double> read_steer::getdtcol(const string& label,const string& col) {
   // get column of a table with header 'col' as string values
   vector<double> ret;
   vector<string> scol = getstcol(label,col);
   for (vector<string>::size_type i = 0; i != scol.size(); i++) {
      string val = scol[i];
      if (!CheckNumber(val.c_str()))
         if (fVerbosity > 1) {
            cout << oW << "Value number " << i << " of table='" << label
                 << "' in column '" << col << "' does not seem to be a numeric number. value=" << val << endl;
         }
      ret.push_back(atof(val.c_str()));
   }
   return ret;
}

vector<string> read_steer::getsthead(const string& label) {
   // get table header
   vector<string> ret = ftableheaders[label];
   if (ret.empty())
      if (fVerbosity > 1) {cout  <<  oW << "Label '" <<   label  << "' was not found in list or has no values." <<  endl;}
   return ret;
}

vector<vector<string> > read_steer::getst(const string& label) {
   // get table values as strings
   vector<vector<string> > ret = ftables[label];
   if (ret.empty() && ftableheaders[label].empty()) {
      if (fVerbosity > 1) {cout  <<  oW << "Label '" <<   label  << "' was not found as a table." <<  endl;}
   } else if (ret.empty()) {
      if (fVerbosity > 2) {cout  <<  oI << "Table '" <<   label  << "' is empty." <<  endl;}
   }
   return ret;
}

vector<vector<double> > read_steer::getdt(const string& label) {
   // get table values as doubles
   vector<vector<double> > ret;
   vector<vector<string> > sf = getst(label);
   for (unsigned int i = 0 ; i<sf.size() ; i++) {
      ret.push_back(vector<double>());
      for (unsigned int j = 0 ; j<sf[i].size() ; j++) {
         string val = sf[i][j];
         if (!CheckNumber(val.c_str()))
            if (fVerbosity > 1) {cout << oW << "Value number (" << i << "," << j << ") of label='" << label << "' does not seem to be a numeric number. value=" << val << endl;}
         ret[i].push_back(atof(val.c_str()));
      }
   }
   return ret;
}


vector<vector<int> > read_steer::getit(const string& label) {
   // get table values as integers
   vector<vector<int> > ret;
   vector<vector<string> > sf = getst(label);
   for (unsigned int i = 0 ; i<sf.size() ; i++) {
      ret.push_back(vector<int>());
      for (unsigned int j = 0 ; j<sf[i].size() ; j++) {
         string val = sf[i][j];
         if (!CheckInt(val.c_str()))
            if (fVerbosity > 1) {cout << oW << "Value number (" << i << "," << j << ") of label='" << label << "' does not seem to be an integer number. value=" << val << endl;}
         ret[i].push_back(atoi(val.c_str()));
      }
   }
   return ret;
}


string read_steer::gets(const string& label) {
   string ret = fstrings[label];
   if (ret=="")
      if (fVerbosity > 1) {cout  <<  oW << "Label '" <<   label  << "' was not found in list or has an empty value." <<  endl;}
   return ret;
}

double read_steer::getd(const string& label) {
   string val = gets(label);
   if (!CheckNumber(val.c_str()))
      if (fVerbosity > 1) {cout << oW << "Value of label='" << label << "' does not seem to be a numeric number. value=" << val << endl;}
   return atof(val.c_str());
}

int read_steer::geti(const string& label) {
   string val = gets(label);
   bool isnan = CheckInt(val.c_str());
   if (!isnan)
      if (fVerbosity > 1) {cout << oW << "Value of label='" << label << "' does not seem to be an integer number. value=" << val << endl;}
   return atoi(val.c_str());
}

bool read_steer::getb(const string& label) {
   return StringToBool(gets(label),label);
}

bool read_steer::StringToBool(const string& sval, const string& label) const {
   if (sval!="0" && sval!="1" && sval!="true" && sval!="false" && sval!="") {
      if (label=="") {
         if (fVerbosity > 1) {cout << oW << "Expecting value '0','1','true', 'false' or no value for boolean values.  value='" << sval << "'. Using 'true'." << endl;}
      } else {
         if (fVerbosity > 1) {cout << oW << "Expecting value '0','1','true', 'false' or no value for boolean values for label=" << label << " and its value='" << sval << "'. Using 'true'." << endl;}
      }
      return true;
   }
   if (sval=="true") return true;
   else if (sval=="false") return false;
   else if (sval=="") return false;
   else return atoi(sval.c_str());
}

void read_steer::AddLabel(const string& key, const string& value) {
   if (fstrings.count(key)>0)
      cout << " # read_steer. Replacing label '" << key << "' with value '" << value << "'." << endl;
   fstrings[key]      = value;
}


void read_steer::AddArray(const string& key, const vector<string>& values) {
   if (ffields.count(key)>0)
      cout << " # read_steer. Replacing label '" << key << "' with an array of size '" << values.size() << "'." << endl;
   ffields[key]      = values;
}


void read_steer::AddTable(const string& key, const vector<string>& header, const vector<vector<string> >& values) {
   if (ftables.count(key)>0 ||  ftableheaders.count(key)>0 )
      cout << " # read_steer. Replacing label '" << key << "' with a table of '" << values.size() << "' columns." << endl;
   ftableheaders[key]      = header;
   ftables[key] = values;
}


void read_steer::AppendToArray(const std::string& key, const std::string& entry) {
   // append one element to an array
   if (ffields.count(key)==0 )
      cout<<" # read_steer. Could not find array '"<<key<<"' yet. Ingoring call please use ::AddArray() instead."<<endl;
   else
      ffields[key].push_back(entry);
}


void read_steer::AppendToTable(const std::string& key, const std::vector<std::string>& entry){
   // append one element to a table
   if (ftableheaders.count(key)==0 )
      cout<<" # read_steer. Could not find table '"<<key<<"' yet. Ingoring call please use ::AddTable() instead."<<endl;
   else {
      //cout<<" # read_steer. Append row with "<<entry.size()<<" entries to table '"<<key<<"'."<<endl;
      ftables[key].push_back(entry);
   }
}


bool read_steer::CheckNumber(const string& str) {
   return str.find_first_of("-+1234567890")==0;
}


bool read_steer::CheckInt(const string& str) {
   return str.find_first_of(".eE")==string::npos && CheckNumber(str);
}



bool read_steer::ParseString(string line) {
   // target variables
   string label;
   string value;

   // keep the string for error messages
   const string orgl=line;

   // parsing statements enclosed in '"'
   // CKR: Put fParseTableMode>0 in brackets since negation otherwise acts on
   //      fParseTableMode variable (warning from clang compiler)!
   if (!(fParseTableMode>0) && !fParseFieldMode)
      value = ParseEnclosedString(line.c_str()); // old
   else
      while (EnclosedStringToOneEntity(line));   // new

   // count
   int i=0;

   // parsing line
   //   char* str = (char*)line.c_str();
   char str[20000];
   strcpy(str,line.c_str());
   for (char* pch=strtok(str,str_sep.c_str());
        pch!=NULL;
        pch=strtok(NULL,str_sep.c_str()),i++) {
      // look for include files
      if (ParseFindString(pch,str_inc)) {
         //cout << "Found Include File. pch=" << pch << endl;
         string incfile,ls,le;
         string spch(pch);
         separatetag(spch, incfile, str_inc);
         separatetag(incfile,ls,":");
         separatetag(ls,le,":");
         ReplaceVariables(incfile);
         unsigned int is = atoi(ls.c_str());
         unsigned int ie = atoi(le.c_str());
         ifstream incstrm;
         incstrm.open(incfile.c_str());
         if (!incstrm) {
            if ( incfile == "fastNLO-warmup.txt" ) {
               if (fVerbosity > 2) {cout << oI << " Could not include specified warmup file " << incfile << ", default warm-up file will be tried!" << endl;}
            } else {
               if (fVerbosity > 0) {cerr << oE << " Could not open include file ('" << incfile << "') from include  statement (" << str_inc << ")." << endl;}
            }
            return EXIT_FAILURE;
         }
         readstrm(incstrm,is,ie,true);
         incstrm.close();
         break;
      }
      // tables
      if (fParseTableMode>0) {
         if (ParseFindString(pch,str_tabend)) {
            fParseTableMode = 0;
            // CKR: The following check doesn't make sense. Usually, one has only ONE label for the innermost dimension
            //      like "----- Array of bin-grid ..." ==> message appears always.
            //            if (!ftablevalues.empty() && !ffieldvalues.empty() && ffieldvalues.size() != ftablevalues[0].size())
            //               cout <<  oI << "Expected a 'table' with " << ffieldvalues.size() << " columns for label '" << ffieldlabel << "', but found a differing number of entries in at least one row." << endl;
            ftableheaders[ffieldlabel] = ffieldvalues;
            ftables[ffieldlabel]     = ftablevalues;
            ffieldvalues.clear();
            ftablevalues.clear();
            ffieldlabel = "";
            return true;
         } else {
            if (fParseTableMode==2) {  // column names
               if (ParseFindString(pch,str_cmt)) break;
               //ffieldvalues.push_back(pch);
               string val = pch;
               ReplaceVariables(val);
               ffieldvalues.push_back(val);
            } else { // table values
               if (ParseFindString(pch,str_cmt))  {
                  if (i==0) return true;  //--fParseTableMode;
                  break;
               }
               if ((int)ftablevalues.size() < fParseTableMode-2)
                  ftablevalues.push_back(vector<string>());
               //ftablevalues[fParseTableMode-3].push_back(pch);
               string val = pch;
               ReplaceVariables(val);
               ftablevalues.back().push_back(val);
            }
         }
      }
      // arrays
      else if (fParseFieldMode) {
         if (ParseFindString(pch,str_arrend)) {   // store field
            fParseFieldMode = false;
            ffields[ffieldlabel] = ffieldvalues;
            ffieldvalues.clear();
            ffieldlabel="";
            return true;
         }
         if (value=="") {  // read single values
            if (ParseFindString(pch,str_cmt)) break;
            //ffieldvalues.push_back(pch);
            string val = pch;
            ReplaceVariables(val);
            ffieldvalues.push_back(val);
         } else { // read enclosed value
            //ffieldvalues.push_back(value);
            string val = value;
            ReplaceVariables(val);
            ffieldvalues.push_back(val);
            break;
         }
      } else {
         // look for a namespace
         if (ParseFindString(pch,str_nmspcend)) {
            return false;
         }
         if (ParseFindString(pch,str_nmspcbeg)) {
            if (fParseIncMode>1) {
               if (fVerbosity > 0) {cerr << oE << "It is not possible to define namespaces in #include(ed) files." << endl;}
            }
            read_steer::initnamespace(ffile,fcurrentfilename,label);
            label = "";
            continue;
         }
         // look for a table
         if (ParseFindString(pch,str_tabbeg)) {
            fParseTableMode = 1;
            if (label=="")
               if (fVerbosity > 1) {cout  <<  oW << "Table found, starting with ' " << str_tabbeg << "' but no label was found." <<  endl;}
            ffieldlabel = label;
            break;
         }
         // look for an array of values
         if (ParseFindString(pch,str_arrbeg)) {
            fParseFieldMode = true;
            if (label=="")
               if (fVerbosity > 1) {cout  <<  oW << "Array found, starting with ' " << str_arrbeg << "' but no label was found." <<  endl;}
            ffieldlabel = label;
            continue;
         }
         // look for comments
         if (ParseFindString(pch,str_cmt)) {
            if (i==1) {
               if (fVerbosity > 1) {cout <<  oW << "Found comment after label ('" << label << "'), but before a value." << endl;}
            }
            break;
         }

         // set label and value
         if (i==0)  label = pch;
         else if (i==1 && value != "") {   // value was already filled with enclosed string
            break;
         } else if (i==1 && value=="") { // set value
            value = pch;
         } else {
            cout  <<  " # read_steer. Error parsing string: "  <<  endl;
            cout  <<  "'"  <<  orgl  <<  "'" << endl;
            cout  <<  " #   Expect two values separated by 'empty spaces' or 'tabstop'." <<  endl;
            cout  <<  " #   Add comments starting with '!' character." << endl;
         }
      }
   } // for parse



   strcpy(str,line.c_str());
   char* pch=strtok(str,str_sep.c_str());
   if (fParseTableMode>0) {
      if (fParseTableMode>2 && pch==NULL) return true;
      fParseTableMode++;
   }
   // if (fParseTableMode>2) {   // check number of columns in table
   //    if (ftablevalues.size() > 1)
   //       // CKR: The following check doesn't make sense, since it´s perfectly fine to
   //       //      have a different number of bins per defining line in steering.
   //       if (ftablevalues.back().size() != ftablevalues[0].size())
   //          cout  << oI << "Table ('" << ffieldlabel << "'): row " << ftablevalues.size()
   //                << " has a different number of columns (n=" <<  ftablevalues.back().size()
   //                << ") than first row (n=" << ftablevalues[0].size() << ")." << endl;
   // }
   //

   // CKR: Change philosophy: Last set value (e.g. from steering) overwrites previously set one (e.g. default)
   if (!fParseFieldMode && fParseTableMode==0) {
      if (fstrings[label] == "") {
         // Set previously empty value
         ReplaceVariables(value);
         fstrings[label] = value;
      } else if (fstrings[label] == value) {
         // Inform that old and new values are identical
         if (fVerbosity > 2) {cout  <<  oI << "Label '" <<   label  << "' already set identically. Keeping value '" << fstrings[label] << "'." <<  endl;}
      } else {
         // Inform that new value overwrites previously set one
         if (fVerbosity > 2) {cout  <<  oI << "Label '" <<   label  << "' already found. Replacing previous value '" << fstrings[label] << "' by new one '" << value << "'." <<  endl;}
         ReplaceVariables(value);
         fstrings[label] = value;
      }
   }
   return true;
}

int read_steer::ReplaceVariables(string& str) {
   // replace all occurences of ${<sth>} by
   // the stringvalue of the label <sth>
   // return the number of replaced variables

   // remove all '$' which could come from enclosed statements in tablemode.

   int ret=0;
   size_t found = str.find(string("${"));
   while (found!=string::npos) {
      size_t end = str.find(string("}"),found);
      if (end==string::npos) {
         if (fVerbosity > 1) {cout  <<  oW << "Start of a variable found with '${', but termination with '}' is missing." << endl;}
         break;
      }
      string var = string(str,found+2,end-found-2);
      string val = fstrings[var];

      if (val=="")
         if (fVerbosity > 1) {cout << oW << "Value of variable ${" << var << "} is empty or was not defined." << endl;}
      str.replace(found,end-found+1,val);
      ret++;
      found = str.find(string("${"));
   }

   // erasing '$' if e.g. in tablemode enclosed strings were used.
   found = str.find(string("$&$"));
   while (found!=string::npos) {
      str.replace(found,3," ");
      found = str.find(string("$&$"));
   }

   // remove empty strings ""
   found = str.find(string("$$%$$"));
   while (found!=string::npos) {
      str.replace(found,5,"");
      found = str.find(string("$$%$$"));
   }

   return ret;
}


string read_steer::ParseEnclosedString(const string& str) const {
   vector<size_t> occ;
   for (size_t found = str.find_first_of('"'); found!=string::npos; found = str.find_first_of('"',found+1))
      occ.push_back(found+1);
   if (occ.size()>0 && occ[0] > str.find(string(str_cmt))) {
      //cout << "Return. '!' before enclosed string: \" in str=" << str << endl;
      return string();
   }
   if (occ.size() >= 2) {
      //if ( occ.size()!=2)
      //cout << oW << "Only lines with exactly two \" symbols are expected in substring '" << str << "'." << endl;
      if (occ[1]-occ[0]-1 > 0) return str.substr(occ[0],occ[1]-occ[0]-1);
      else return "$$%$$";
   } else {
      if (!occ.empty())
         if (fVerbosity > 1) {cout << oW << "Only lines more than two \" symbols are allowed in substring '" << str << "'." << endl;}
   }
   return string();
}

bool read_steer::EnclosedStringToOneEntity(string& str) const {
   // replace within first occurances of "" all empty-spaces and tab-stops by $
   // return false, of nothing was replaced
   // return true of, sth. was replaced

   size_t one = str.find('"');
   if (one==string::npos) return false;
   if (one!=string::npos) str.erase(one,1);
   size_t two = str.find('"');
   if (one!=string::npos) str.erase(two,1);
   string substr = str.substr(one,two-one);
   const string sub = str.substr(one,two-one);

   if (one==two) { // it is only "" -> we tag it.
      str.insert(two,"$$%$$");
   }

   while (substr.find_first_of(str_sep)!=string::npos) {
      size_t pos = substr.find_first_of(str_sep);
      substr.replace(pos,1,"$&$");
   }

   // replace ins str, sub with substr
   size_t pos = str.find(sub);
   size_t len = sub.size();
   str = str.replace(pos,len,substr);

   return true;
}


bool read_steer::ParseFindString(const string& str, const string& tag) const {
   //return strncmp(str,tag.c_str(),tag.size())==0;
   return (str.find(tag)==0);
}

int read_steer::separatetag(string& vallhs, string& valrhs, const string& sep) {
   // separate a string, according to separation string sep;
   // input:  vallhs, sep
   // output: vallhs, valrhs
   // if no separator was found return vallhs and valrhs untouched
   // otherwise vallhs is replaced by left-hand-side of separator, and
   //  valrhs is replace by right-hand-side of separator
   // return -1, if no separator was found,
   // return position of separato string in vallhs
   const size_t pos = vallhs.find(sep);
   const string temp=vallhs;
   if (pos!=string::npos) {
      vallhs=temp.substr(0,pos);
      valrhs=temp.substr(pos+sep.size(),temp.size());
      return (int)pos;
   }
   return -1;
}

int read_steer::cmdlinetag(const char* arg, string& label, string& value) {
   label = string(arg);
   int ret = separatetag(label,value,"=");
   return ret;
}


bool read_steer::parsecommandline(vector<string> argv) {
   vector<char*> argc;//(argv.size(), new char[400]);
   // C++11
   // for ( unsigned int i = 0 ; i<argv.size() ; i++ )
   //    argc.push_back((char*)&argv[0]); // C++11
   //  return ret = parsecommandline(argc.size(),&argc[0]);

   for ( unsigned int i = 0 ; i<argv.size() ; i++ ){
      cout<<"input: "<<argv[i]<<endl;
      argc.push_back(new char[argv[i].size()+1]);
      std::copy(argv[i].begin(), argv[i].end(), argc[i]);
      argc[i][argv[i].size()] = '\0';
   }
   bool ret = parsecommandline(argc.size(),&argc[0]);
   for ( unsigned int i=0 ; i<argc.size() ;i ++ ) delete[] argc[i];
   return ret;
}

bool read_steer::parsecommandline(int argc,char** argv) {
   bool gotfile = false;
   map<string,string> cmdvals;
   for (int i=0; i<argc; i++) {
      string val, lab;
      cout << "i=" << i << "\targv[i]=" << argv[i] << endl;
      int suc=cmdlinetag(argv[i],lab,val);
      if (suc>0) {
         cout << " found cmd line tag. label=" << lab << "\tvalue=" << val << endl;
         if (lab=="steerfile") {
            string fID = stdID;
            int pos = separatetag(val,fID,"->");
            cout << "strfile. "  << "val=" << val << "\tfID=" << fID << endl;
            if (pos>=0)  cout << " # read_steer. Info. Reading new steerfile '" << val << "'." << endl;;
            readfile(val,fID);
            gotfile=true;
         } else {
            cmdvals[lab] = val;
         }
      }
   }

   // if already StdId found -> replace values
   // else create new StdID namespace with cmd-line arguments
   if (instances==NULL) read_steer::Steering(stdID);
   for (map<string, string>::const_iterator ii=cmdvals.begin(); ii!=cmdvals.end(); ++ii) {
      //! use like key=value:namespace
      //string val = (*ii).second;
      //string fID = stdID;
      //separatetag(val, fID ,":");
      //! use like namespace:key=value
      string val = (*ii).second;
      string fID = stdID;
      string key = (*ii).first;
      separatetag(key, fID ,"::");
      if ( fID!=stdID )swap(key,fID);
      read_steer::Steering(fID)->AddLabel(key, val);
   }
   return gotfile;
}


set<string> read_steer::GetAvailableLabels() const {
   std::set<std::string> ret;
   std::map<std::string,std::string> map = fstrings;
   for (std::map<std::string,std::string>::iterator it=map.begin(); it!=map.end(); ++it)
      ret.insert(it->first);
   return ret;
}


set<string> read_steer::GetAvailableArrrays() const {
   std::set<std::string> ret;
   std::map<std::string,std::vector<std::string> > map = ffields;
   for (std::map<std::string,std::vector<std::string> >::iterator it=map.begin(); it!=map.end(); ++it)
      ret.insert(it->first);
   return ret;
}


set<string> read_steer::GetAvailableTables() const {
   std::set<std::string> ret;
   std::map<std::string,std::vector<std::vector<std::string> > > map = ftables;
   for (std::map<std::string,std::vector<std::vector<std::string> > > ::iterator it=map.begin(); it!=map.end(); ++it)
      ret.insert(it->first);
   return ret;
}
