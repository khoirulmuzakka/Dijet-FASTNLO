{
#include "colors.h";

//
// Load needed libraries
//
// cout << LIGHTBLUE << "\nReading gnu scientific library."<<NOCOLOR<<endl;
// gSystem.Load("~/lib/libgsl.so");
cout << LIGHTBLUE << "Reading headers and classes from 'libfastnlo' into memory."<<NOCOLOR<<endl;
gSystem.Load("~/lib/libfastnlo.so");
puts("Done.");


//
// Define and set my graphics style
//
cout << LIGHTBLUE <<"Setting graphical defaults."<<NOCOLOR<<endl;
TStyle *Style = gROOT->GetStyle("Plain");
Style->SetHistLineWidth(1);
Style->SetTextFont(43);
Style->SetTextSizePixels(12.0);
Style->SetLabelFont(42,"XYZ");
Style->SetLabelSize(0.05,"YZ");
Style->SetLabelSize(0.05,"X");
Style->SetNdivisions(505,"XYZ");
Style->SetStatFont(42);
Style->SetStatFontSize(.03);
Style->SetTitleFont(43,"XYZ");
Style->SetTitleSize(15,"YZ");
Style->SetTitleSize(15,"X");
Style->SetTitleOffset(2.3,"Y");
Style->SetTitleOffset(2.2,"X");
Style->SetTitleOffset(2.3,"Z");
Style->SetPadBorderSize(0);
Style->SetCanvasBorderSize(0);
Style->SetOptStat(0);
Style->SetOptFit(1);
Style->SetOptTitle(0);
Style->SetTitleX(0.05);
Style->SetTitleW(0.9);
Style->SetTitleH(0.08);
Style->SetTitleFont(43);
Style->SetTitleFontSize(12.0);
Style->SetPadTopMargin(0.1);
Style->SetPadBottomMargin(.14);
Style->SetPadLeftMargin(.14);
Style->SetPadRightMargin(.04);

Style->SetStatW(0.28);
Style->SetStripDecimals(kFALSE);
Style->SetStatFormat("5.3g");
Style->SetPaintTextFormat("4.1g");
Style->SetFillColor(10);
Style->SetPaperSize(21.,29.5);

Style->SetPalette(1,0);

TStyle *NewStyle = new TStyle("Big","Big");
Style->Copy(*NewStyle);
NewStyle->SetLabelSize(0.2,"XYZ");
NewStyle->SetTitleX(1.0);

gROOT->SetStyle("Plain");
gROOT->ForceStyle(kTRUE);
puts("Done.");

cout << LIGHTGREEN << "What may I do for you?\n" << NOCOLOR << endl;
}

