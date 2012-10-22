if ( $?TEXINPUTS ) then
    setenv TEXINPUTS "${TEXINPUTS}:./eps//:./jpg//:./pdf//:./png//:./texinputs//:"
else
    setenv TEXINPUTS "./eps//:./jpg//:./pdf//:./png//:./texinputs//:"
endif
