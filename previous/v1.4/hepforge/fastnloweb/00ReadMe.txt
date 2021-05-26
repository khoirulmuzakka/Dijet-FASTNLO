
****************************************************************
* special feature in code at CEDAR:  alpha_s in common-block
****************************************************************

To have a flexible alpha_s(Mz) value, the CEDAR web-code has 
been modified. The alphas value is passed from the main file
"main-cedar.f" through the commonblock /FNCEDAR/ to
the interface routine "fn-interface-cedar.f".


main-cedar.f:      COMMON /FNCEDAR/ ASMZ
fn-interface-cedar.f:      COMMON /FNCEDAR/ ASMZ


