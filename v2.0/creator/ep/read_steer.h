// Author: Daniel Britzger
// DESY, 17/07/2012
#ifndef __read_steer_h__
#define __read_steer_h__ 1

//**********************************************************************************
//
//     read_steer.h
//     Tiny reading tool to read steering values from one or more steering files.
//
//     This class reads in values, which are stored in a file. New variables
//     can be included without changes of the steering class. 
//
//     Features
//     ------------------------------
//       o  Following types are supported:
//            - Single values
//              bool, int, double, string (with empty spaces)
//            - Arrays
//              int, double, string (with empty spaces)
//            - Tables/matrices
//              int, double, string (no empty spaces)
//       o  Multiple files can be read in and handled individually or together. 
//       o  Namespaces can be defined.
//       o  Local variables within the steer-file can be defined.
//
//
//     Initalize the steering
//     ------------------------------
//     Set the filename and initilize the read_steer class by using:
//        read_steer::readfile(string filename)
//
//
//     Single values
//     ------------------------------
//     To access the values, simply use the adequate getter functions
//     for the desired variable type and the label of this variable.
//     To speed up the code and to avoid repeated string comparisions, 
//     use static variables, e.g if you want to access the value in 
//     your steering file with the label 'pi' or 'name', use:
//        static double pi   = read_steer::getdouble("pi");
//     or 
//        static string name = read_steer::getstring("name");
//        static int    age  = read_steer::getint("age");
//        static bool   sex  = read_steer::getbool("female");
//
//     Labels are case sensitive!
//
//
//     Syntax of steering file 
//     ------------------------------
//     The steering file can consist of an arbitrary number of lines, where
//     the syntax should follow:
//         <label>		<value>		[!comment]
//     where 'label' and 'value' are necessary tags and comments are
//     beginning with the '!' character and are ignored by the read_steer class.
//     As seperator between the <label> and the <value> empty spaces or tabstops
//     are recognized. Complete lines can beginn with '!' to mark comments.
//     If string values should contain empty spaces, enclose them in double quotes
//     like:
//         Name			"Peter Higgs"
//     Boolean values can be assigned by 0, 1, true or false, e.g.
//         WithHiggs		true
//
//
//     Arrays
//     ------------------------------
//     To read in an array of values, assign a label and enclose the following 
//     values in curly brackets { }, with a leading empty space [" {" and "}"].
//     Within curly brackets, each separated (by an empty space or tabstop)
//     value is stored in the array as an element, or each line when double quotes are
//     used (only one occurence of double quotes per line is recognized).
//     The steering file should look like:
//
//	  !Numbers 1-9 are read into an array (9 elements)
//        Array1 {			!array starts here
//          1 2 3 4 5			! integers from 1 to 5
//          6 7 8 9
//        }
//        !Eleven Names of famous musicians (11 elements)
//        FamousMusicians {
//          John Paul Ringo George	!The Beatles
//          Beethoven Bach Mozart	!Famous componists
//          Mick Keith Ron Charlie	!The Rolling Stones
//        }
//        !Full sentences or documentations (2 string-elements)
//        Array3 {
//          "Hello World!"		!Sentences are great
//          "Was the first scream."
//        }
//
//     To access the arrays use e.g.:
//         static vector<string> musicians = read_steer::getstringarray("FamousMusicians");
//         static vector<double> nums      = read_steer::getdoublearray("Array1");
//         static vector<int>    ints      = read_steer::getintarray("Array1");
//
//
//     Tables and matrices
//     ------------------------------
//     Tables and matrices are tagged by ' {{' and '}}' in the steering file.
//     No double quotes (e.g. "text") are allowed as table values. The first row
//     of a table is always expected to be the row-header. The row headers are
//     separated by empty spaces or tabstops. If matrices are necessary
//     keep the first line empty or add a comment there. The steering file should look like:
//    
//        Crossections {{
//		Q2min	Q2max	cs[pb]	stat[%]			!header tags should not contain emtpy spaces
//		100	200	22.12	1.2
//		200	300	12.72	2.7
//		300	500	23.22	5.3
//	  }}
//        Participants {{
//		Name		Surname		Country		! first line is always the header
//		Obama		Barack		U.S.A.
//		Merkel		Angela		Germany
//		Benedikt	XVI		Vatican
//	  }}
//        Matrix {{
//	        !the first line is ignored. Keep it empty.
//		11 12
//		21 22
//	   }}
//
//	To access the table use e.g.:
//	    static vector<vector<string> > guys = read_steer::getstringtable("Participants");
//	    static vector<vector<double> > cs   = read_steer::getdoubletable("Crossections");
//	    static vector<vector<int> >    mat  = read_steer::getinttable("Matrix");
//      To access the table header use:
//         static vector<string>	   head = read_steer::gettableheader("Crossections);
//	To access a single column of a table use:
//	    static vector<double> xs		= read_steer::getdoublecolumn("Crossections","cs[pb]");
//	    static vector<string> nick		= read_steer::getstringcolumn("Participants","Surname");
//
//
//     Multiple steering files.
//     ------------------------------
//     In case multiple steering files are necessary, each steering file must
//     be assigned a unique 'steerID' if variable names (labels)  are identically.
//          read_steer::readfile("file1.steer","file1")
//          read_steer::readfile("anotherfile.steer","constants")
//
//     To access values, pass the steerID to the getter methods, e.g.:
//          static double pi   = read_steer::getdouble("pi","constants");
//          static string name = read_steer::getstring("name","file1");
//	    static vector<vector<string> > ConfIchepNames = read_steer::getstringcolumn("Participants","Surname","file1")
//     You can access the values at any place within your code.
//
//     If different labels should be read in from multiple files, just call
//          read_steer::readfile("file1.steer");
//          read_steer::readfile("file2.steer");
//     and access the variables without the usage of the steerID.
//
//
//     Namespaces
//     ------------------------------
//     Instead of using multiple files for reading identical labels for
//     various occasions, one can use namespaces instead. Namespaces are
//     handled identically to multiple files, but can be defined within
//     one single steering file. Each namespace is assigned a steerID.
//     Namespaces are defined by a label, which is used as the steerID 
//     and start with the '{{{' tag and end with the '}}}' tag.
//     A steerfile could look like:
//	 HostInstitute		CERN		! standard variabel
//       ATLAS {{{				! namespace ATLAS starts here
//          length	45			! define variables as usual
//	    height	22
//          weight	7000
//	    Crossection {{			! also tables are possible
//	 	bin	cs[pb]	stat[%]
//	 	1	32.2	1.2
//		2	12.2	3.2
//          }}
//	 }}}					! namespace ATLAS ends here
//       CMS {{{
//          length	21
//	    height	16
//          weight	12500
//	    Crossection {{
//		bin	cs[pb]	stat[%]
//		1	33.1	0.8
//		2	13.6	3.4
//          }}
//	 }}}
//
//     To access the values, use the steerID which is the label of the namespace
//          static double ATLASheight	  = read_steer::getdouble("height","ATLAS");
//          static double CMSheight	  = read_steer::getdouble("height","CMS");
//          static vector<double> CMSxs	  = read_steer::getdoublecolumn("Crossection","cs[pb]","CMS");
//          static vector<double> ATLASxs = read_steer::getdoublecolumn("Crossection","cs[pb]","ATLAS");
//
//     Warning: Namespace steerID and file steerID might conflict if identically!
//     It is NOT possible to read in multiple files, wherein identical namespaces are define!
//
//
//     Script-like Variables
//     ------------------------------
//     It is often neessary to read in identical substrings, e.g. if
//     many different files are located in the same folder. To simplify the
//     structure of the steering file, 'script-like' variables can be used.
//
//     It is possible to access previously defined variables foo by ${foo}.
//     An example steering file can look like:
//	     !Home directories of famous physicists
//           HomeDir			/afs/cern.ch/user
//           UserEinstein		${HomeDir}/e/einstein
//           UserNewton			"${HomeDir}/i/isaac"
//
//
//     Printing
//     ------------------------------
//     Print all steering information in SingleFileMode by calling
//	    read_steer::print();
//     If multiple files are used, print all information using:
//          read_steer::printall();
//     or just the information of one steerID:
//          read_steer::print("constants");
//
//
//     Warning
//     ------------------------------
//     This tool is basically based on string comparisions to identify
//     the steering values. Do not call it too often within your code in order to
//     avoid speed problems. The best way to use the steering values is to assign the
//     values to static variables.
//
//
//     Not existing labels or steerID
//     ------------------------------
//     If you access an element which was not read in from the steering file,
//     this label is automatically added to the list of elements with value zero.
//     If a steerID is accesed, which is not identified, a new steerID is added
//     to the list without any values.
//
//
//     D. Britzger
//     daniel.britzger@desy.de
//
//
//**********************************************************************************
 
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
#include <vector>

using namespace std;

class read_steer {
private:
   static map<string,read_steer*> instances;
   static const string stdID;
   read_steer(){;};
   read_steer(const read_steer& ) {;};
public:
   ~read_steer() {;};
   static read_steer* Steering(string steerID=stdID);			// get an object!
   static void destroy();						// destroy all instances
   
   // static member function
   static void readfile(string filename,string steerID=stdID) {	// set the steer-filename
      read_steer::Steering(steerID)->inits(filename); }
   // getters
   // values
   static bool getbool(string label,string steerID=stdID) {
      return read_steer::Steering(steerID)->getb(label); }
   static int getint(string label,string steerID=stdID) {
      return read_steer::Steering(steerID)->geti(label); }
   static double getdouble(string label,string steerID=stdID) {
      return read_steer::Steering(steerID)->getd(label); }
   static string getstring(string label,string steerID=stdID) {
      return read_steer::Steering(steerID)->gets(label); }
   // arrays
   static vector<int> getintarray(string label,string steerID=stdID) {
      return read_steer::Steering(steerID)->getif(label); }
   static vector<double> getdoublearray(string label,string steerID=stdID) {
      return read_steer::Steering(steerID)->getdf(label); }
   static vector<string> getstringarray(string label,string steerID=stdID) {
      return read_steer::Steering(steerID)->getsf(label); }
   // tables header
   static vector<string> gettableheader(string label,string steerID=stdID) {
      return read_steer::Steering(steerID)->getsthead(label); }
   // tables/matrices
   static vector<vector<int> > getinttable(string label,string steerID=stdID) {
      return read_steer::Steering(steerID)->getit(label); }
   static vector<vector<double> > getdoubletable(string label,string steerID=stdID) {
      return read_steer::Steering(steerID)->getdt(label); }
   static vector<vector<string> > getstringtable(string label,string steerID=stdID) {
      return read_steer::Steering(steerID)->getst(label); }
   // table columns
   static vector<int> getintcolumn(string label,string column ,string steerID=stdID) {
      return read_steer::Steering(steerID)->getitcol(label,column); }
   static vector<double> getdoublecolumn(string label,string column ,string steerID=stdID) {
      return read_steer::Steering(steerID)->getdtcol(label,column); }
   static vector<string> getstringcolumn(string label,string column ,string steerID=stdID) {
      return read_steer::Steering(steerID)->getstcol(label,column); }
   
   const static void printall();						// print values of all files
   const static void print(string steerID=stdID);				// print values

public:
   // getters for single instance
   // values
   const bool getb(string label);
   const int geti(string label);
   const double getd(string label);
   const string gets(string label);
   // arrays
   const vector<int> getif(string label);
   const vector<double> getdf(string label);
   const vector<string> getsf(string label);
   // tables/matrices
   const vector<string> getsthead(string label);
   const vector<vector<int> > getit(string label);
   const vector<vector<double> > getdt(string label);
   const vector<vector<string> > getst(string label);
   const vector<int> getitcol(string label,string col);
   const vector<double> getdtcol(string label,string col);
   const vector<string> getstcol(string label,string col);
   
   // controls
   void inits(string filename);
   int initnmspc(ifstream& strm, string filename);
   const void prt();
   static void initnamespace(ifstream& strm,string filename, string steerID=stdID) {	// set the steer-filename
      read_steer::Steering(steerID)->initnmspc(strm,filename); }

private:
   int read_stdin(string filename);
   int readstrm(ifstream& strm);
   bool ParseString(string value);
   const bool ParseFindString(const string str, const string tag);
   const string ParseEnclosedString(const string);
   const int ReplaceVariables(string& value);
   const bool CheckNumber(const string str);
   const bool CheckInt(const string str);

   map<string,string> fstrings;
   map<string,vector<string> > ffields;
   map<string,vector<vector<string> > > ftables;
   map<string,vector<string> > ftableheaders;
   bool fParseFieldMode;
   int  fParseTableMode;
   string ffieldlabel;
   vector<string> ffieldvalues;
   vector<vector<string> > ftablevalues;
   string ffilename;
   string fcurrentfilename;
   ifstream ffile;

   static const string str_cmt;
   static const string str_sep;
   static const string str_arrbeg;
   static const string str_arrend;
   static const string str_tabbeg;
   static const string str_tabend;
   static const string str_nmspcbeg;
   static const string str_nmspcend;

   static const string oW, oI, oE;
};


map<string,read_steer*> read_steer::instances = map<string,read_steer*>();
const string read_steer::stdID = "SingleFileMode";
const string read_steer::str_sep=" \t";
const string read_steer::str_cmt="!";
const string read_steer::str_arrbeg="{";
const string read_steer::str_arrend="}";
const string read_steer::str_tabbeg="{{";
const string read_steer::str_tabend="}}";
const string read_steer::str_nmspcbeg="{{{";
const string read_steer::str_nmspcend="}}}";
const string read_steer::oW=" # read_steer. Warning. ";
const string read_steer::oI=" # read_steer. Info. ";
const string read_steer::oE=" # read_steer. ERROR. ";

read_steer* read_steer::Steering(string steerID)
{
   // get singleton class
   if ( !instances[steerID] ){
      if ( steerID.compare(stdID)!=0 ) 
	 cout<<oI<<"Initalizing new read_steer instance with steerID = '"<<steerID<<"'."<<endl;
      instances[steerID] = new read_steer();
   }
   return instances[steerID];
}


void read_steer::inits(string filename){
   //    if ( ffilename != "" )
   //       cout<<oW<<"Filename already set (old="<<ffilename<<", new="<<filename<<"). Is the used steerID unique?"<<endl;
   if ( filename == "" )
      cout<<oW<<"No filename specified."<<endl;
   if ( ffilename !="" ) ffilename+=", ";
   ffilename += filename;
   fcurrentfilename = filename;
   read_stdin(fcurrentfilename);
}


void read_steer::destroy()
{
   for( map<string, read_steer*>::iterator ii=instances.begin(); ii!=instances.end(); ++ii){
      if ( (*ii).second ) delete (*ii).second;
      instances.erase((*ii).first);
   }
}

const void read_steer::printall()
{
   const string linesep = " +----------------------------------------------------------------------------+\n";
   const string l = " | ";
   cout<<linesep;
   cout<<l<<"    read_steer. Printing all steering information.                         |"<<endl;
   cout<<linesep;
   for( map<string, read_steer*>::iterator ii=instances.begin(); ii!=instances.end(); ++ii)
      {
	 cout<<l<<endl;
	 cout<<l<<"steerID = '"<<(*ii).first<<"'"<<endl;
	 cout<<linesep;
	 (*ii).second->prt();
	 cout<<linesep;
      }
}


const void read_steer::print(string steerID)
{
   const string linesep = " +----------------------------------------------------------------------------+\n";
   const string l = " | ";
   cout<<linesep;
   cout<<l<<"    read_steer. Printing steering information of steerID = '"<<steerID<<"'"<<endl;
   cout<<linesep;
   read_steer::Steering(steerID)->prt();
   cout<<linesep;
}

const void read_steer::prt()
{
   const string l = " | ";
   // filename
   printf("%s%-30s\t%s\n",l.c_str(),"Filename(s)",ffilename.c_str());
   cout<<l<<endl;;
   //single values
   cout<<l<<"Single values"<<endl;
   for( map<string,string>::iterator ii=fstrings.begin(); ii!=fstrings.end(); ++ii)
      if ( (*ii).first!="") 
	 printf("%s   %-27s\t%s\n",l.c_str(),(*ii).first.c_str(),(*ii).second.c_str());
   // arrays
   if ( !ffields.empty() ) {
      cout<<l<<endl;
      cout<<l<<"Arrays"<<endl;
      for( map<string,vector<string> >::iterator ii=ffields.begin(); ii!=ffields.end(); ++ii){
	 cout<<l<<"   "<<(*ii).first<< " {"<<endl;
	 for ( unsigned int j = 0 ; j<(*ii).second.size() ; j++ )
	    cout <<l<<"     ["<<j<<"]\t"<<(*ii).second[j]<<endl;
	 cout <<l<<"   }"<<endl;
      }
   }
   // tables
   if ( !ftables.empty() ){
      cout<<l<<endl;
      cout<<l<<"Tables/Matrices"<<endl;
      for( map<string,vector<string> >::iterator ii=ftableheaders.begin(); ii!=ftableheaders.end(); ++ii){
	 cout<<l<<"   "<<(*ii).first<< " {{"<<endl;
	 cout<<l<<"     [H]\t";
	 for ( unsigned int j = 0 ; j<(*ii).second.size() ; j++ )
	    printf("%-10s",(*ii).second[j].c_str());
	 cout<< endl;
	 vector<vector<string> > tab = ftables[(*ii).first];
	 for ( unsigned int ll = 0 ; ll<tab.size() ; ll++ ){
	    cout<<l<<"     ["<<ll<<"]\t";
	    for ( unsigned int j = 0 ; j<tab[ll].size() ; j++ )
	       printf("%-10s",tab[ll][j].c_str());
	    cout<<endl;
	 }
	 cout <<l<<"   }}"<<endl;
      }
   }
}

int read_steer::initnmspc(ifstream& strm, string filename){
   if ( ffilename !="" ) ffilename+=", ";
   ffilename += filename;
   //ffilename = filename;
   return readstrm(strm);
}

int read_steer::readstrm(ifstream& strm){
   if (!strm){
      cerr<<oE<<"This is not a valid stream."<<endl;
      return EXIT_FAILURE;
   }
   string lineread;
   int nlines=0;
   fParseFieldMode = false;
   fParseTableMode = 0;
   ffieldlabel = "";
   while( std::getline(strm, lineread)) {
      bool goon = ParseString(lineread);
      nlines++;
      if ( !goon ) break;
   }
   return nlines; // return nlines including comments
}

int read_steer::read_stdin(string filename)
{
   //If the steering has alread been read -> do nothing
   ffile.open(filename.c_str());
   if (!ffile){
      cerr<<oE<<" Could not open file ('"<<filename<<"')."<<endl;
      return EXIT_FAILURE;
   }
   int n = readstrm(ffile); 
   ffile.close();
   return n;
};
  

const vector<int> read_steer::getif(string label){
   vector<int> ret;
   vector<string> sf = ffields[label];
   for ( unsigned int i = 0 ; i<sf.size() ; i++ ){
      string val = sf[i];
      bool isnan = CheckInt(val.c_str());
      if ( !isnan ) 
	 cout<<oW<<"Value number "<< i<<" of label='"<<label<<"' does not seem to be an integer number. value="<<val<<endl;
      ret.push_back( atoi(val.c_str()));
   }
   return ret;
}

const vector<double> read_steer::getdf(string label){
   vector<double> ret;
   vector<string> sf = ffields[label];
   for ( unsigned int i = 0 ; i<sf.size() ; i++ ){
      string val = sf[i];
      bool isnan = CheckNumber(val.c_str());
      if ( !isnan ) 
	 cout<<oW<<"Value number "<< i<<" of label='"<<label<<"' does not seem to be a numeric number. value="<<val<<endl;
      ret.push_back( atof(val.c_str()));
   }
   return ret;
}

const vector<string> read_steer::getsf(string label){
   vector<string> ret = ffields[label];
   if (ret.empty())
      cout << oW<<"Label '"<<  label <<"' was not found in list or has no values."<< endl;
   return ret;
}

const vector<string> read_steer::getstcol(string label,string col){
   // get column of a table with header 'col' as string values
   vector<string> ret;
   vector<string> head = getsthead(label);
   vector<vector<string> > tab = getst(label);
   for(vector<string>::size_type i = 0; i != head.size(); i++) {
      if ( col.compare(head[i])==0 ){
	 for(vector<string>::size_type j = 0; j != tab.size(); j++) {
	    ret.push_back(tab[j][i]);
	 }
	 return ret;
      }
   }
   cout<<oW<<"Column '"<<col<<"' was not found in table '"<<label<<"'."<<endl;
   return ret;
}

const vector<int> read_steer::getitcol(string label,string col){
   // get column of a table with header 'col' as string values
   vector<int> ret;
   vector<string> scol = getstcol(label,col);
   for(vector<string>::size_type i = 0; i != scol.size(); i++) {
      string val = scol[i];
      if ( !CheckInt(val.c_str()) ) 
	 cout<<oW<<"Value number "<<i<<" of table='"<<label
	     <<"' in column '"<<col<<"' does not seem to be an integer number. value="<<val<<endl;
      ret.push_back(atoi(val.c_str()));
   }
   return ret;
}

const vector<double> read_steer::getdtcol(string label,string col){
   // get column of a table with header 'col' as string values
   vector<double> ret;
   vector<string> scol = getstcol(label,col);
   for(vector<string>::size_type i = 0; i != scol.size(); i++) {
      string val = scol[i];
      if ( !CheckNumber(val.c_str()) ) 
	 cout<<oW<<"Value number "<<i<<" of table='"<<label
	     <<"' in column '"<<col<<"' does not seem to be a numeric number. value="<<val<<endl;
      ret.push_back(atof(val.c_str()));
   }
   return ret;
}

const vector<string> read_steer::getsthead(string label){
   // get table header
   vector<string> ret = ftableheaders[label];
   if (ret.empty())
      cout << oW<<"Label '"<<  label <<"' was not found in list or has no values."<< endl;
   return ret;
}

const vector<vector<string> > read_steer::getst(string label){
   // get table values as strings
   vector<vector<string> > ret = ftables[label];
   if (ret.empty())
      cout << oW<<"Label '"<<  label <<"' was not found in list or has no values."<< endl;
   return ret;
}

const vector<vector<double> > read_steer::getdt(string label){
   // get table values as doubles
   vector<vector<double> > ret;
   vector<vector<string> > sf = getst(label);
   for ( unsigned int i = 0 ; i<sf.size() ; i++ ){
      ret.push_back(vector<double>());
      for ( unsigned int j = 0 ; j<sf[i].size() ; j++ ){
	 string val = sf[i][j];
	 if ( !CheckNumber(val.c_str()) ) 
	    cout<<oW<<"Value number ("<<i<<","<<j<<") of label='"<<label<<"' does not seem to be a numeric number. value="<<val<<endl;
	 ret[i].push_back( atof(val.c_str()));
      }
   }
   return ret;
}


const vector<vector<int> > read_steer::getit(string label){
   // get table values as integers
   vector<vector<int> > ret;
   vector<vector<string> > sf = getst(label);
   for ( unsigned int i = 0 ; i<sf.size() ; i++ ){
      ret.push_back(vector<int>());
      for ( unsigned int j = 0 ; j<sf[i].size() ; j++ ){
	 string val = sf[i][j];
	 if ( !CheckInt(val.c_str()) )
	    cout<<oW<<"Value number ("<<i<<","<<j<<") of label='"<<label<<"' does not seem to be an integer number. value="<<val<<endl;
	 ret[i].push_back( atof(val.c_str()));
      }
   }
   return ret;
}


const string read_steer::gets(string label){
   string ret = fstrings[label];
   if (ret=="")
      cout << oW<<"Label '"<<  label <<"' was not found in list or has an empty value."<< endl;
   return ret;
}

const double read_steer::getd(string label){
   string val = gets(label);
   if ( !CheckNumber(val.c_str()) ) 
      cout<<oW<<"Value of label='"<<label<<"' does not seem to be a numeric number. value="<<val<<endl;
   return atof(val.c_str());
}

const int read_steer::geti(string label){
   string val = gets(label);
   bool isnan = CheckInt(val.c_str());
   if ( !isnan )
      cout<<oW<<"Value of label='"<<label<<"' does not seem to be an integer number. value="<<val<<endl;
   return atoi(val.c_str());
}

const bool read_steer::getb(string label){
   string sval = gets(label);
   int val = atoi(sval.c_str());
   if ( sval!="0" && sval!="1" && sval!="true" && sval!="false")
      cout<<oW<<"Expecting value '0','1','true' or 'false' for boolean values. label='"<<label<<"', value='"<<sval<<"'."<<endl;
   if ( sval=="true" ) return true;
   else if ( sval=="false") return false;
   else return val;
}


const bool read_steer::CheckNumber(const string str){
   return str.find_first_of("-+1234567890")==0;
}


const bool read_steer::CheckInt(const string str){
   return str.find_first_of(".eE")==string::npos && CheckNumber(str);
}



bool read_steer::ParseString(string line)
{
   // target variables
   string label;
   string value;

   // keep the string for error messages
   const string orgl=line;

   // parsing statements enclosed in '"' 
   value = ParseEnclosedString(line.c_str());

   // count
   int i=0;

   // parsing line
   char* str = (char*)line.c_str();
   for(char* pch=strtok (str,str_sep.c_str());
       pch!=NULL;
       pch=strtok(NULL,str_sep.c_str()),i++)
      {
	 if ( fParseTableMode>0 ) {
	    if ( ParseFindString(pch,str_tabend) ) { // store table
	       fParseTableMode = 0;
	       if ( !ftablevalues.empty() && !ffieldvalues.empty() && ffieldvalues.size() != ftablevalues[0].size() )
		  cout<< oW<<"Table ('"<<ffieldlabel<<"'): header has a different number of columns (n="
		      <<ffieldvalues.size()<<") than table (n="<<ftablevalues[0].size()<<")."<<endl;
	       ftableheaders[ffieldlabel] = ffieldvalues;
	       ftables[ffieldlabel]	= ftablevalues;
	       ffieldvalues.clear();
	       ftablevalues.clear();
	       ffieldlabel = "";
	       return true;
	    } else {
	       if ( fParseTableMode==2 ){ // column names
		  if ( ParseFindString(pch,str_cmt) ) break;
		  //ffieldvalues.push_back(pch);
		  string val = pch;
		  ReplaceVariables(val);
		  ffieldvalues.push_back(val);
	       }
	       else { // table values
		  if ( (int)ftablevalues.size() < fParseTableMode-2 )
		     ftablevalues.push_back(vector<string>());
		  //ftablevalues[fParseTableMode-3].push_back(pch);
		  string val = pch;
		  ReplaceVariables(val);
		  ftablevalues[fParseTableMode-3].push_back(val);
	       }
	    }
	 }
	 else if ( fParseFieldMode ) {
	    if ( ParseFindString(pch,str_arrend) ) { // store field
	       fParseFieldMode = false;
	       ffields[ffieldlabel] = ffieldvalues;
	       ffieldvalues.clear();
	       ffieldlabel="";
	       return true;
	    }
	    if ( value=="" ){ // read single values
	       if ( ParseFindString(pch,str_cmt) ) break;
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
	 } 
	 else {
	    // look for a namespace
	    if ( ParseFindString(pch,str_nmspcend) ){
	       return false;
	    }
	    if ( ParseFindString(pch,str_nmspcbeg) ){
	       read_steer::initnamespace(ffile,fcurrentfilename,label);
	       label = "";
	       continue;
	    }
	    // look for a table
	    if ( ParseFindString(pch,str_tabbeg) ){
	       fParseTableMode = 1;
	       if ( label=="" ) 
		  cout << oW<<"Table found, starting with ' "<<str_tabbeg<<"' but no label was found."<< endl;
	       ffieldlabel = label;
	       break;
	    }
	    // look for an array of values
	    if (ParseFindString(pch,str_arrbeg) ){
	       fParseFieldMode = true;
	       if ( label=="" ) 
		  cout << oW<<"Array found, starting with ' "<<str_arrbeg<<"' but no label was found."<< endl;
	       ffieldlabel = label;
	       break;
	    }
	    // look for comments
	    if ( ParseFindString(pch,str_cmt) ) {
	       if ( i==1 ) { cout<< oW<<"Found comment after label ('"<<label<<"'), but before a value."<<endl;}
	       break;
	    }

	    // set label and value
	    if ( i==0 )  label = pch;
	    else if ( i==1 && value != "" ) { // value was already filled with enclosed string
	       break;
	    } 
	    else if ( i==1 && value=="" ){ // set value
	       value = pch;
	    }
	    else {
	       cout << " # read_steer. Error parsing string: " << endl;
	       cout << "'" << orgl << "'"<<endl;
	       cout << " #   Expect two values separated by 'empty spaces' or 'tabstop'."<< endl;
	       cout << " #   Add comments starting with '!' character."<<endl;
	    }
	 }
      }

   if ( fParseTableMode>2 ) { // check number of columns in table
      if ( ftablevalues.size() > 1 ) 
	 if ( ftablevalues.back().size() != ftablevalues[0].size() )
	    cout <<oW<<"Table ('"<<ffieldlabel<<"'): row "<<ftablevalues.size()+1
		 <<" has a different number of columns (n="<< ftablevalues.back().size()
		 <<") than first row (n="<<ftablevalues[0].size()<<")."<<endl;
   }
   if ( fParseTableMode>0 ) fParseTableMode++;
   if ( !fParseFieldMode && fParseTableMode==0) {
      if ( fstrings[label] == "" ) {
	 ReplaceVariables(value);
	 fstrings[label] = value;
      }  else 
	 cout << oW<<"Label '"<<  label <<"' already found. Ignoring value ('"<<value<<"'.)"<< endl;
   }
   return true;
}

const int read_steer::ReplaceVariables(string& str){
   // replace all occurences of ${<sth>} by 
   // the stringvalue of the label <sth> 
   // return the number of replaced variables
   int ret=0;
   size_t found = str.find(string("${"));
   while ( found!=string::npos ) {
      size_t end = str.find(string("}"));
      if ( end==string::npos ) {
	 cout << oW<<"Start of a variable found with '${', but termination with '}' is missing."<<endl;
	 break;
      }
      string var = string(str,found+2,end-found-2);
      string val = fstrings[var];
      if ( val=="" )
	 cout<<oW<<"Value of variable ${"<<var<<"} is empty or was not defined."<<endl;
      str.replace(found,end-found+1,val);
      ret++;
      found = str.find(string("${"));
   }
   return ret;
}


const string read_steer::ParseEnclosedString(const string str){
   vector<size_t> occ;
   for ( char* pch=strchr(str.c_str(),'"');pch!=NULL;pch=strchr(pch+1,'"'))
      occ.push_back(pch-str.c_str()+1);
   if ( occ.size() == 2 )
      return str.substr(occ[0],occ[1]-occ[0]-1);
   else {
      if ( !occ.empty() ) 
	 cout<<oW<<"Only lines with exactly two \" symbols are allowed in substring '"<<str<<"'."<<endl;
      return string();
   }
}


const bool read_steer::ParseFindString(const string str, const string tag)
{
   //return strncmp(str,tag.c_str(),tag.size())==0;
   return ( str.find(tag)==0 );
}

#endif
