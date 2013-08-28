#include "fastNLOCoeffAddEvalBase.h"


fastNLOCoeffAddEvalBase::fastNLOCoeffAddEvalBase(int NObsBin, fastNLOReader* reader) : fastNLOCoeffAddBase(NObsBin) , fReader(reader) {

}

fastNLOCoeffAddEvalBase::fastNLOCoeffAddEvalBase(const fastNLOCoeffBase& base , int iLOord, fastNLOReader* reader) : fastNLOCoeffAddBase(base) , fReader(reader) {

}

fastNLOCoeffAddEvalBase::fastNLOCoeffAddEvalBase() : fastNLOCoeffAddBase() {
}
