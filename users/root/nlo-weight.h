#ifndef __NLO_NLO_WEIGHT_H__
#define __NLO_NLO_WEIGHT_H__ 1

//   standard headers
#include <cmath>
#include <cstring>



namespace nlo {
  

  template<unsigned int _Size>
  class weight
  {
  public:
    //   constructors
    weight() {
      std::memset(_M_sub, 0, _Size*sizeof(double));
    }
    
    explicit weight(const double *arr) {
      std::memcpy(_M_sub, arr, _Size*sizeof(double));
    }

    //   copy constructor
    weight(const weight<_Size>& w) {
      std::memcpy(_M_sub, w._M_sub, _Size*sizeof(double));
    }
    
    //   assignment
    weight<_Size>& operator=(const weight<_Size>& w) {
      if(this != &w) std::memcpy(_M_sub, w._M_sub, _Size*sizeof(double));
      return *this;
    }
    
    //   computed assignments
    weight& operator*=(const double& x) { 
      for(double *i = _M_sub; i < _M_sub+_Size; *(i++) *= x);
       return *this;
    }
    
    weight& operator/=(const double& x) { 
      for(double *i = _M_sub; i < _M_sub+_Size; *(i++) /= x);
       return *this;
    }
    
    weight& operator*=(const weight& x) { 
      for(unsigned i = 0; i < _Size; i++) 
	_M_sub[i] *= x[i];
      return *this;
    }
 
   weight& operator/=(const weight& x) { 
      for(unsigned i = 0; i < _Size; i++) 
	_M_sub[i] /= x[i];
      return *this;
    }
 
    weight& operator+=(const weight& x) { 
      for(unsigned i = 0; i < _Size; i++) 
	_M_sub[i] += x[i];
      return *this;
    }

    weight& operator-=(const weight& x) { 
      for(unsigned i = 0; i < _Size; i++) 
	_M_sub[i] -= x[i];
      return *this;
    }

    //   iterators
    double * begin() { return _M_sub;} 
    double * end() { return _M_sub;} 

    const double * begin() const { return _M_sub;} 
    const double * end() const { return _M_sub;} 

    //   element access
    const double& operator[](unsigned int n) const { return _M_sub[n];}
    double& operator[](unsigned int n) { return _M_sub[n];}
    
    //  size of the weight
    unsigned size() const { return _Size;}
    
  private:
    double _M_sub[_Size];
  };

  
  template<unsigned int _Size>
  inline weight<_Size> operator-(const weight<_Size>& a) {
    return weight<_Size>(a) *= -1.0;
  }
  
  template<unsigned int _Size> inline 
  weight<_Size> operator+(const weight<_Size>& a, const weight<_Size>& b) {
    return weight<_Size>(a) += b;
  }

  template<unsigned int _Size> inline 
  weight<_Size> operator-(const weight<_Size>& a, const weight<_Size>& b) {
    return weight<_Size>(a) -= b;
  }

  template<unsigned int _Size> inline 
  weight<_Size> operator/(const weight<_Size>& a, const weight<_Size>& b) {
    return weight<_Size>(a) /= b;
  }

  template<unsigned int _Size> inline 
  weight<_Size> operator*(const weight<_Size>& a, const weight<_Size>& b) {
    return weight<_Size>(a) *= b;
  }
  
  template<unsigned int _Size> inline 
  weight<_Size> operator*(const weight<_Size>& a, double b) {
    return weight<_Size>(a) *= b;
  }
  
  template<unsigned int _Size> inline 
  weight<_Size> operator*(double b, const weight<_Size>& a) {
    return weight<_Size>(a) *= b;
  }

  template<unsigned int _Size> inline 
  weight<_Size> operator/(const weight<_Size>& a, double b) {
    return weight<_Size>(a) /= b;
  }
  
  template<unsigned int _Size> 
  weight<_Size> operator/(double a, const weight<_Size>& b) 
  { 
    weight<_Size> res;
    for(unsigned i = 0; i < _Size; i++) 
      res[i] = a/b[i];
    return res;
  }

  template<unsigned int _Size> 
  weight<_Size> apply(const weight<_Size>& x, double f(double)) 
  { 
    weight<_Size> res;
    for(unsigned i = 0; i < _Size; i++) 
      res[i] = f(x[i]);
    return res;
  }
  
  template<unsigned int _Size> 
  weight<_Size> apply(const weight<_Size>& x, double f(const double&)) 
  { 
    weight<_Size> res;
    for(unsigned i = 0; i < _Size; i++) 
      res[i] = f(x[i]);
    return res;
  }
  
  
  template<unsigned int _Size> 
  inline weight<_Size> abs(const weight<_Size>& x) {
    return apply(x, std::abs);
  }
  
  template<unsigned int _Size> 
  inline weight<_Size> acos(const weight<_Size>& x) {
    return apply(x, std::acos);
  }
  
  template<unsigned int _Size> 
  inline weight<_Size> asin(const weight<_Size>& x) {
    return apply(x, std::asin);
  }
  
  template<unsigned int _Size> 
  inline weight<_Size> atan(const weight<_Size>& x) {
    return apply(x, std::atan);
  }
  
  template<unsigned int _Size> 
  inline weight<_Size> cos(const weight<_Size>& x) {
    return apply(x, std::cos);
  }
  
  template<unsigned int _Size> 
  inline weight<_Size> cosh(const weight<_Size>& x) {
    return apply(x, std::cosh);
  }
  
  template<unsigned int _Size> 
  inline weight<_Size> exp(const weight<_Size>& x) {
    return apply(x, std::exp);
  }
  
  template<unsigned int _Size> 
  inline weight<_Size> log(const weight<_Size>& x) {
    return apply(x, std::log);
  }
  
  template<unsigned int _Size> 
  inline weight<_Size> log10(const weight<_Size>& x) {
    return apply(x, std::log10);
  }
  
  template<unsigned int _Size> 
  inline weight<_Size> sin(const weight<_Size>& x) {
    return apply(x, std::sin);
  }
  
  template<unsigned int _Size> 
  inline weight<_Size> sinh(const weight<_Size>& x) {
    return apply(x, std::sinh);
  }
  
  template<unsigned int _Size> 
  inline weight<_Size> sqrt(const weight<_Size>& x) {
    return apply(x, std::sqrt);
  }
  
  template<unsigned int _Size> 
  inline weight<_Size> tan(const weight<_Size>& x) {
    return apply(x, std::tan);
  }
  
  template<unsigned int _Size> 
  inline weight<_Size> tanh(const weight<_Size>& x) {
    return apply(x, std::tanh);
  }
  
  template<unsigned int _Size> 
  weight<_Size> atan2(const weight<_Size>& x, const weight<_Size>& y) 
  { 
    weight<_Size> res;
    for(unsigned i = 0; i < _Size; i++) 
      res[i] = std::atan2(x[i], y[i]);
    return res;
  }

  template<unsigned int _Size> 
  weight<_Size> atan2(const weight<_Size>& x, const double& y)
  { 
    weight<_Size> res;
    for(unsigned i = 0; i < _Size; i++) 
      res[i] = std::atan2(x[i], y);
    return res;
  }
  
  template<unsigned int _Size> 
  weight<_Size> atan2(const double& x, const weight<_Size>& y) 
  { 
    weight<_Size> res;
    for(unsigned i = 0; i < _Size; i++) 
      res[i] = std::atan2(x, y[i]);
    return res;
  }
  
  template<unsigned int _Size> 
  weight<_Size> pow(const weight<_Size>& x, const weight<_Size>& y)
  { 
    weight<_Size> res;
    for(unsigned i = 0; i < _Size; i++) 
      res[i] = std::pow(x[i], y[i]);
    return res;
  }
   
  template<unsigned int _Size> 
  weight<_Size> pow(const weight<_Size>& x, const double& y) 
  { 
    weight<_Size> res;
    for(unsigned i = 0; i < _Size; i++) 
      res[i] = std::pow(x[i], y);
    return res;
  }
  
  template<unsigned int _Size> 
  weight<_Size> pow(const double& x, const weight<_Size>& y)
  { 
    weight<_Size> res;
    for(unsigned i = 0; i < _Size; i++) 
      res[i] = std::pow(x, y[i]);
    return res;
  }
  
  template<unsigned int_Size> 
  weight<_Size> pow(const weight<_Size>& x, int y) 
  { 
    weight<_Size> res;
    for(unsigned i = 0; i < _Size; i++) 
      res[i] = std::pow(x[i], y);
    return res;
  }

  template<unsigned int _Size>
  inline bool IsFinite(const weight<_Size> & x) 
  {
    for(unsigned int i = 0; i < _Size; i++)
      if(!::finite(x[i])) return false;
    return true;
  }
  

}   //  namespace nlo


#endif

