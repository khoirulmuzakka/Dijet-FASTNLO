if [ $TEXINPUTS ]; then
    export TEXINPUTS=${TEXINPUTS}:./eps//:./jpg//:./pdf//:./png//:./texinputs//:
  else
    export TEXINPUTS=./eps//:./jpg//:./pdf//:./png//:./texinputs//:
fi
