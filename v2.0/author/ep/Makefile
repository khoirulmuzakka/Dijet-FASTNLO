
# H1 HERA I inclusive jets, DIS

fnh1001.la: pdf-cteq6.h fnh1001.cc 
	create-nlojet-user fnh1001.cc alew.f -o fnh1001 -Wx,-O0,-Wall,-I${HOME}/include/lhpdf,-I../../common -Wl,-L${HOME}/lib/fastnlo,-llhpdf_cteq6,-lfastnloNlojet

# H1 HERA I inclusive jets, low Q2 DIS

fnh2003.la: pdf-cteq6.h fnh2003.cc 
	create-nlojet-user fnh2003.cc alew.f -o fnh2003 -Wx,-O0,-Wall,-I${HOME}/include/lhpdf,-I../../common -Wl,-L${HOME}/lib/fastnlo,-llhpdf_cteq6,-lfastnloNlojet


# H1 DESY 02-225, photo production

fnh2101d.la: pdf-cteq6.h pdf-photo-dummy.h fnh2101d.cc 
	create-nlojet-user fnh2101d.cc  -o fnh2101d -Wx,-O0,-Wall,-I${HOME}/include/lhpdf,-I../../common -Wl,-L${HOME}/lib/fastnlo,-llhpdf_cteq6,-lfastnloNlojet
fnh2101r.la: pdf-grv-cteq6.h pdf-hhc-dummy.h fnh2101r.cc
	create-nlojet-user fnh2101r.cc pdf-grv-cteq6.cc grvphoton.f -o fnh2101r -Wx,-O0,-Wall,-I${HOME}/include/lhpdf,-I../../common -Wl,-L${HOME}/lib/fastnlo,-llhpdf_cteq6,-lfastnloNlojet 

# ZEUS DESY 02-228, photo production

fnh2102d.la: pdf-cteq6.h pdf-photo-dummy.h fnh2102d.cc 
	create-nlojet-user fnh2102d.cc -o fnh2102d -Wx,-O0,-Wall,-I${HOME}/include/lhpdf,-I../../common -Wl,-L${HOME}/lib/fastnlo,-llhpdf_cteq6,-lfastnloNlojet
fnh2102r.la: pdf-grv-cteq6.h pdf-hhc-dummy.h fnh2102r.cc
	create-nlojet-user fnh2102r.cc pdf-grv-cteq6.cc grvphoton.f -o fnh2102r -Wx,-O0,-Wall,-I${HOME}/include/lhpdf,-I../../common -Wl,-L${HOME}/lib/fastnlo,-llhpdf_cteq6,-lfastnloNlojet 
