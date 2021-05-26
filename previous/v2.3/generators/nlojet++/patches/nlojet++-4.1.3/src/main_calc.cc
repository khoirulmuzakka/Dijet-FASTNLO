//  Copyright (C) 2002 Zoltan Nagy
//
//  This program is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#include <cstdlib>
#include <ctime>
#include <getopt.h>

#include <string>
#include <cstring>
#include <iostream>
#include <sstream>


//----- the used namespaces -----
using namespace std;

#include "nlojet++.h"
#include "ltdl.h"


//----- print out the help -----
void main_calc_help()
{
  cout
    <<"Usage: nlojet++ --calcualte -n name -c born|nlo|full -u module.la \n" 
    <<"       [-d dir] [--saving-mode txt|bin] [-m | --mchel]\n"
    <<"       [-A num | --alpha num] [--max-event nevent]\n"
    <<"       [--save-after nsave] [-T num, --time-rate num]\n"
    <<"\n       nlojet++ -- calculate -h|--help \n\n"
    
    <<"Options:\n" 
    <<"  -h, --help                 print this message\n"
    <<"  -n name                    name of the run\n"
    <<"  -c born|nlo|full           contributions\n"
    <<"  -u module.la               module file contains the user defined functions\n"
    <<"  -d dir                     output directory\n"
    <<"  --txt                      formated output (default is binary output) \n"
    <<"  -m, --mchel                turn off the Monte Carlo helicity sum\n"
    <<"  -A num, --alpha num        non-physical cut in the dipole terms (0,1]\n"
    <<"  --max-event number         max number of the events\n"
    <<"  --save-after number        save after every 'number' events\n"
    <<"  -T num, --time-rate num    time rate between the n+1 and n parton\n"
    <<"                             integral (num = 1,...10 integer, default: 1)"
    <<endl;
}


int main_calc(int argc, char **argv)
{
  //----- long options -----
  option longopt[] =
    {{"txt",        required_argument, 0,  0 },
     {"mchel",      no_argument,       0, 'm'},
     {"alpha",      required_argument, 0, 'A'},
     {"max-event",  required_argument, 0,  0 },
     {"save-after", required_argument, 0,  0 },
     {"time-rate",  required_argument, 0, 'T'},
     {"help",       no_argument,       0, 'h'},
     {0, 0, 0, 0}
    };

  //------ default options -----
  main_input in;
  const char *name = "run", *dir= "./output", *contr = contbl[0], *usr = 0;
  unsigned int seed = ((unsigned int) time(NULL)) & 5313465U; 
 
  //------ argument parsing -----
  while(1) {
    int optidx = 0;      
    int c = getopt_long(argc, argv, "mhc:n:d:u:s:A:T:", longopt, &optidx);
    if(c == -1) break;
    
    switch(c) {
    case 'h': main_calc_help(); exit(0); break;
    case 'n': name = optarg; break;
    case 'd': dir  = optarg; break;
    case 'u': usr  = optarg; break;
    case 's': seed = atoi(optarg); break;
    case 'm': in.mchel = false; break;
    case 'A': in.alpha = atof(optarg); break;
    case 'T': in.time  = atoi(optarg); break;
    case 'c': contr = optarg; break;
    case 0:
      switch(optidx) {
      case 0: in.txtout = true; break;
      case 3: in.nevent = atoi(optarg); break;
      case 4: in.nsave  = atoi(optarg); break;
      }
      break;
    case '?': break;
    default: throw; break;
    }
  }
  
  //----- user module -----
  if(!usr) throw "no user module specified";

  //----- contribution -----
  for(unsigned ic = 0; contbl[ic]; ic++)
    if(strcmp(contbl[ic], contr) == 0) {
      in.contr = ic;
      break;
    }
  
  if(in.contr == -1) {
    cerr<<"Unknown contribution type `"<<contr<<"`. "
	<<"The deafault contribution (`"
	<<contbl[0]<<"`) will be calculated."<<endl;
    in.contr = 0;
  }
  
  //----- make the output directory -----
  make_dir(dir);
  
  //----- output file name template -----
  stringstream strs;
  strs<<dir; if(dir[strlen(dir)-1] != '/') strs<<"/";
  strs<<name<<"-@PROC@-"<<contbl[in.contr]<<"-@NJET@jet"<<'\0';
  in.outfile = strs.str();  
  
  //  the time rate parameter must be greater than zero 
  //  and should be less than 10
  if(in.time == 0U) in.time = 1U;
  if(in.time > 10U) in.time = 10U;
  
  //  initialize the random seed
  srand(seed);
  
  // KR: Add random seed printout
  cout<<"User module        : "<<usr<<"\n"
      <<"Contribution       : "<<contbl[in.contr]<<"\n"
      <<"MC helicity sum    : "<<(in.mchel ? "true" : "false")<<"\n"
      <<"Dipole cut         : "<<in.alpha<<"\n"
      <<"Number of events   : "<<in.nevent<<"\n"
      <<"Output directory   : "<<dir<<"\n"
      <<"Saving mode        : "<<(in.txtout ? "text" : "binary")<<"\n"
      <<"Save after         : "<<in.nsave<<" events\n"
      <<"Time rate In+1:In  : "<<in.time<<":1\n"
      <<"Random seed        : "<<seed<<"\n"
      <<endl;
  
  //----- open the user defined file -----
  if(lt_dlinit() != 0) throw lt_dlerror();
  
  lt_dlhandle handle = lt_dlopen(usr);
  if (!handle) {
    cerr<<"can't open the module "<<usr<<"!\n";
    throw lt_dlerror();
  }
  
  const lt_dlinfo *info = lt_dlgetinfo(handle);
  if (!info) {
    cerr<<"can't get module info: "<<lt_dlerror()<<endl;
    return 1;
  }
  
  if(info->name) cout<<"module name: "<<info->name<<"\n";
  else cout<<"module is not a libtool module\n";
  cout<<"module filename: "<<info->filename<<"\n";
  cout<<"module reference count: "<<info->ref_count<<endl;
  
  //----- try to find the symbol table -----
  cout<<"\nTrying to find the symbol table :\n"
      <<"user_defined_functions[]        : ";
  
  in.symtable = (nlo_symbol *) lt_dlsym(handle, "user_defined_functions");
  if(in.symtable) cout<<"     OK\n"<<endl;
  else {
    cout<<"     FAIL\n"<<endl;
    throw "There is no well defined user module!";
  }
  
  //----- find the process type -----
  const char *proc = (const char *) in.find_symbol("procindex");
  in.outfile.replace(in.outfile.find("@PROC@"), 6, proc);
  
  int procidx = -1;
  while(proctbl[++procidx].name)
    if(std::strcmp(proc, proctbl[procidx].name) == 0) break;
  
  if(proctbl[procidx].name == 0) 
    throw "Undefined process. Check your user module!";
  
  //----- starting the calculation -----
  (*(proctbl[procidx].module_calc))(in);
  
  //   close the dlopend files
  lt_dlclose(handle);
  lt_dlexit();
  
  return 0;
}

