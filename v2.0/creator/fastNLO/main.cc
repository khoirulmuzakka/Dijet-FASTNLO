// D. Britzger 28.06.2013

// #include <cfloat>
// #include <cmath>
// #include <cstdlib>
// #include <iomanip>
#include <vector>
#include <utility>
#include <iostream>
#include "fastNLOInterpolCatmulRom.h"
#include "fastNLOTable.h" 
#include "fastNLOCreate.h"

#include "fastNLOLHAPDF.h"

//__________________________________________________________________________________________________________________________________
int main()
{
   // namespaces
   using namespace std;

   
   fastNLOCreate tt("fnlttbar000.str");
   tt.Print();
   return 0;
   /*
   fastNLOTable table();
   
   fastNLOTable tab("../ep/finaltables/fnh5001.tab");
   //fastNLOTable tab("../hadron/tt.tab");
   tab.Print();

   fastNLOTable tab2("/afs/desy.de/user/b/britzger/www/FastNLOForKumar/fastNLOv2/fnl1014_v2_all.tab");
   tab2.Print();

   fastNLOTable tab3("../hadron/output/files_fnl5002ak06/fnl5002ak06_final_20G_2jet.tab");
   tab3.Print();

   */
   //fastNLOReader r1("../ep/finaltables/fnh5001.tab");
   //fastNLOLHAPDF  r1("../ep/finaltables/fnh5001.tab","CT10.LHgrid",0);
   //fastNLOLHAPDF  r1("/afs/desy.de/user/b/britzger/www/FastNLOForKumar/fastNLOv2/fnl1014_v2_all.tab","CT10.LHgrid",0);

   //fastNLOLHAPDF  r1("../hadron/output/files_fnl5002ak06/fnl5002ak06_final_20G_2jet.tab","CT10.LHgrid",0);
   

   fastNLOLHAPDF  r1("/afs/desy.de/user/b/britzger/www/FastNLOForKumar/fastNLOv2/fnl1014_v2_all.tab","CT10.LHgrid",0);
   r1.CalcCrossSection();
   r1.PrintCrossSections();
   
//    r1.Print();
//    r1.SetUnits(fastNLO::kAbsoluteUnits);
//    //   r1.SetContributionON(fastNLO::kFixedOrder, 1, false);
//    //r1.RunFastNLODemo();

   return 0;

   fastNLOCreate t("FastNLOSteering.str");
   

   fastNLOInterpolCatmulRom inter(10,100);
   inter.MakeGrids(fastNLOGrid::kLinear, 20);
   vector<pair<int,double> > p = inter.GetNodeValues(30);
   for ( int i = 0 ; i<4 ; i++ )
      cout<<p[i].first<<"\t"<<p[i].second<<endl;

   cout<< ">>> " << 10 <<endl;
   p = inter.GetNodeValues(10);
   for ( int i = 0 ; i<4 ; i++ )
      cout<<p[i].first<<"\t"<<p[i].second<<endl;

   cout<< ">>> " << 10.01 <<endl;
   p = inter.GetNodeValues(10.01);
   for ( int i = 0 ; i<4 ; i++ )
      cout<<p[i].first<<"\t"<<p[i].second<<endl;

   cout<< ">>> " << 11 <<endl;
   p = inter.GetNodeValues(11);
   for ( int i = 0 ; i<4 ; i++ )
      cout<<p[i].first<<"\t"<<p[i].second<<endl;

   cout<< ">>> " << 99 <<endl;
   p = inter.GetNodeValues(99);
   for ( int i = 0 ; i<4 ; i++ )
      cout<<p[i].first<<"\t"<<p[i].second<<endl;
   cout<< ">>> " << 100 <<endl;
   p = inter.GetNodeValues(100);
   for ( int i = 0 ; i<4 ; i++ )
      cout<<p[i].first<<"\t"<<p[i].second<<endl;

   cout<<"sanity!"<<endl;
   cout<< ">>> " << 0 <<endl;
   p = inter.GetNodeValues(0);
   for ( int i = 0 ; i<4 ; i++ )
      cout<<p[i].first<<"\t"<<p[i].second<<endl;

   cout<< ">>> " << 10000 <<endl;
   p = inter.GetNodeValues(10000);
   for ( int i = 0 ; i<4 ; i++ )
      cout<<p[i].first<<"\t"<<p[i].second<<endl;
   
   return 0;
}

