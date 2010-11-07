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
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
#ifndef __NLO_NLO_BASIC_USER_H__
#define __NLO_NLO_BASIC_USER_H__ 1

//   Standard includes
#include <typeinfo>
#include <iostream>
#include <fstream>
#include <ctime>

//   nlojet++ includes
#include <distribution.h>
#include <bits/hep-func.h>
#include <bits/nlo-jetfunc.h>


namespace nlo {

  //
  //   Abstract class for defining basic_user types 
  //
  template<class _Jet>
  class basic_user_base : public _Jet
  {
  public:
    //   destructor
    virtual ~basic_user_base() {}
  
    //   set the output (file name saving mode)
    virtual void phys_output(const char *fname, unsigned long nsave = 10000UL, bool txt = false) {
      _M_ofile = fname; _M_txt = txt; _M_nsave = nsave;
    }
    
    virtual void phys_output(const std::basic_string<char>& fname, unsigned long nsave = 10000UL, bool txt = false) {
      _M_ofile = fname; _M_txt = txt; _M_nsave = nsave;
    }
    
    virtual void phys_output(unsigned long nsave, bool txt) {
      _M_txt = txt; _M_nsave = nsave;
    }
     
    // operations at the end of the event and saving the results
    virtual void end_of_event();
    
  protected:
    //   do any operetaions at the end of the event
    virtual void operations_at_the_end_of_event() = 0;
    
    //   write some information about the type of the sample
    virtual void write_typeinfo(std::ostream&) = 0;
    
    //   save the results in txt and binary mode
    virtual void write_result_txt(std::ostream&) = 0;
    virtual void write_result_bin(std::ostream&) = 0;
    
    //  constructors and destructor
    basic_user_base() 
      : _M_txt(false), _M_nsave(10000UL), _M_ofile("nlojet++.out"),
	_M_iter(0UL), _M_start_time(std::time(0)){}
    
  private:
    //    data members
    bool                      _M_txt;    
    unsigned long             _M_nsave;
    std::basic_string<char>   _M_ofile;
    
    unsigned long _M_iter;
    time_t _M_start_time;
  };
  
  template<class _Jet>
  void basic_user_base<_Jet>::end_of_event()
  {
    this -> operations_at_the_end_of_event();
    
    if((++_M_iter) % _M_nsave == 0) {
      std::fstream ofile(_M_ofile.c_str(), std::ios::out);
      
      ofile<<_M_txt<<"\n";
      ofile<<_M_iter<<"\n";
      this -> write_typeinfo(ofile);
      
      if(_M_txt) {
		ofile.setf(std::ios::scientific, std::ios::floatfield);
		ofile.precision(14);
		this -> write_result_txt(ofile);
      } else this -> write_result_bin(ofile);
      
      ofile.flush();
      ofile.close();
      
      time_t hh, mm, ss = std::time(0) - _M_start_time;
      
      hh  = ss/3600L;
      ss -= hh*3600L;
      mm  = ss/60L;
      ss -= mm*60L;
      
      std::cout<<"--->     "<<_M_iter<<"   "
	       <<(hh < 10 ? "0" : "")<<hh
	       <<(mm < 10 ? ":0" : ":")<<mm
	       <<(ss < 10 ? ":0" : ":")<<ss<<std::endl;
    }
  }
  

  template<class _From, class _To>
  struct default_conversion 
    : public std::unary_function<_From, _To> 
  {
    _To operator()(const _From& x) { return _To(x);}
  };
  
  template<class _Tp>
  struct default_conversion<_Tp, _Tp>
    : public std::unary_function<_Tp, _Tp> 
  {
    const _Tp& operator()(const _Tp& x) { return x;}
  };
  
  
  //
  //    The general basic_user type
  //
  template<class _Jet, typename _Point, typename _Weight = typename _Jet::amp_type::weight_type,
           typename _Sample = _Weight, template<class _Xp, class _Yp> class _Conv = default_conversion,
           class _Traits = sample_traits<_Sample>, class _PointTraits = distpoint_traits<_Point> >
  class basic_user 
    : public virtual basic_user_base<_Jet>
  {
  public:
    //   public types
    typedef _Point  point_type;
    typedef _Weight weight_type;
    typedef _Sample sample_type;
    typedef _Traits sample_traits;
    typedef _PointTraits point_traits;

    typedef distbook<_Sample, _Point, _Traits, _PointTraits> distbook_type;
    typedef _Conv<_Weight, _Sample> conversion_type;
    
    //   destructor
    virtual ~basic_user() {}
    
  protected:
    //   do any operetaions at the end of the event
    void operations_at_the_end_of_event() {
      ++_M_dist;
    }
    
    //   save the results in txt and binary mode
    void write_result_txt(std::ostream& os) {
      os<<_M_dist<<"\n";
    }
    
    void write_result_bin(std::ostream& os) {
      write(os, _M_dist);
    }
   
    void write_typeinfo(std::ostream& os) {
      os<<typeid(distbook_type).name()<<"\n";
    }

  private:
    //    data members
    distbook<_Sample, _Point, _Traits, _PointTraits> _M_dist;
 
  protected:
    //   It is visible for the inhereted class. One can have some advanced 
    //   conversion that requires special initialization.  
    static conversion_type _S_conv;

  protected:
    //   create a distribution
    void phys(int id, const char *name, unsigned int nbin, 
	      const point_type *base) 
    {
      _M_dist.create(id, name, nbin, base);
    } 
    
  private:
    template<class _Ptr>
    struct unary_func_conv 
    {
      typedef typename _Ptr::argument_type argument_type;
      typedef typename _Ptr::result_type result_type;
      
      unary_func_conv(const _Ptr& f) : _M_func(f) {}
      
      _Sample operator()(argument_type x) const {
		return _S_conv(_M_func(x));
      }
 
      const _Ptr& _M_func;
    };

  protected:
    //   fill the distribution
    template<class _Arg, class _Res> 
    void physfill(int id, const std::pointer_to_unary_function<_Arg, _Res>& f) {
      typedef unary_func_conv<std::pointer_to_unary_function<_Arg, _Res> > func_conv_type;
      _M_dist.accumulate(id, ptr_fun(func_conv_type(f), &func_conv_type::operator()));
    }
    
    template<class _Sp, class _Arg, class _Res> 
    void physfill(int id, const pointer_to_unary_member_function<_Sp, _Arg, _Res>& f) {
      typedef unary_func_conv<pointer_to_unary_member_function<_Sp, _Arg, _Res> > func_conv_type;
      _M_dist.accumulate(id, ptr_fun(func_conv_type(f), &func_conv_type::operator()));
    }
    
    template<class _Sp, class _Arg, class _Res> 
    void physfill(int id, const pointer_to_unary_const_member_function<_Sp, _Arg, _Res>& f) {
      typedef unary_func_conv<pointer_to_unary_const_member_function<_Sp, _Arg, _Res> > func_conv_type;
      _M_dist.accumulate(id, ptr_fun(func_conv_type(f), &func_conv_type::operator()));
    }
   
    template<class _Arg> 
    void physfill(int id, const std::pointer_to_unary_function<_Arg, _Sample>& f) {
      _M_dist.accumulate(id, f);
    }
    
    template<class _Sp, class _Arg> 
    void physfill(int id, const pointer_to_unary_member_function<_Sp, _Arg, _Sample>& f) {
      _M_dist.accumulate(id, f);
    }
    
    template<class _Sp, class _Arg> 
    void physfill(int id, const pointer_to_unary_const_member_function<_Sp, _Arg, _Sample>& f) {
      _M_dist.accumulate(id, f);
    }
    
    template<class _Arg1, class _Arg2, class _Res> 
    void physfill(int id, const std::pointer_to_binary_function<_Arg1, _Arg2, _Res>& f, _Arg1 x, const _Weight& w) {
      _M_dist.accumulate(id, f, x, _S_conv(w));
    }
    
    template<class _Sp, class _Arg1, class _Arg2, class _Res> 
    void physfill(int id, const pointer_to_binary_member_function<_Sp, _Arg1, _Arg2, _Res>& f, _Arg1 x, const _Weight& w) {
      _M_dist.accumulate(id, f, x, _S_conv(w));
    }
    
    template<class _Sp, class _Arg1, class _Arg2, class _Res> 
    void physfill(int id, const pointer_to_binary_const_member_function<_Sp, _Arg1, _Arg2, _Res>& f, _Arg1 x, const _Weight& w) {
      _M_dist.accumulate(id, f, x, _S_conv(w));
    }
  };
  
  
  template<class _Jet, typename _Point, typename _Weight, typename _Sample, template<class _Xp, class _Yp> class _Conv, class _Traits, class _PointTraits>
  _Conv<_Weight, _Sample> basic_user<_Jet, _Point, _Weight, _Sample, _Conv, _Traits, _PointTraits>::_S_conv;
  
  
  //  Specialization:
  //         Scalar distribution (In this case the _Point = void)
  //
  template<class _Jet, typename _Weight, typename _Sample, template<class _Xp, class _Yp> class _Conv, class _Traits>
  class basic_user<_Jet, void, _Weight, _Sample, _Conv, _Traits, distpoint_traits<void> >
    : public virtual basic_user_base<_Jet>
  {
  public:
    //   public types
    typedef void    point_type;
    typedef _Weight weight_type;
    typedef _Sample sample_type;
    typedef _Traits sample_traits;
    typedef distbook<_Sample, void, _Traits, distpoint_traits<void> > distbook_type;
    typedef _Conv<_Weight, _Sample> conversion_type;
		 
    //   destructor
    virtual ~basic_user() {}
		 
  private:
    //    data members
    distbook_type _M_dist;   
		 
  protected:
    //   It is visible for the inhereted class. One can have some advanced 
    //   conversion that requires special initialization.  
    static conversion_type _S_conv;
		 
  protected:
    //   do any operetaions at the end of the event
    void operations_at_the_end_of_event() {
      ++_M_dist;
    }
		 
    //   save the results in txt and binary mode
    void write_result_txt(std::ostream& os) {
      os<<_M_dist<<"\n";
    }
		 
    void write_result_bin(std::ostream& os) {
      write(os, _M_dist);
    }
		 
    void write_typeinfo(std::ostream& os) {
      os<<typeid(distbook_type).name()<<"\n";
    }
		 
    //   create a scalar distribution
    void phys(int id, const char *name) {
      _M_dist.create(id, name);
    }
		 
    //   fill the scalar distributions
    void physfill(int id, const _Weight& w) {
      _M_dist.accumulate(id, _S_conv(w));
    }
  };
  
  template<class _Jet, typename _Weight, typename _Sample, template<class _Xp, class _Yp> class _Conv, class _Traits>
  _Conv<_Weight, _Sample> basic_user<_Jet, void, _Weight, _Sample, _Conv, _Traits, distpoint_traits<void> >::_S_conv;
  
  
  //  Specialization:
  //         One dimensional distribution (In this case the _Point = double)
  //
  
  //   Helper functions to compute Dirac delta distribution and step function distribution
  namespace __helper_basic_user {
    extern std::pointer_to_binary_function<double, double, double> _G_step_dbl;
  }  //  namespace __helper_basic_user
 
  
  template<class _Jet, typename _Weight, typename _Sample, template<class _Xp, class _Yp> class _Conv, class _Traits, class _PointTraits>
  class basic_user<_Jet, double, _Weight, _Sample, _Conv, _Traits, _PointTraits> 
    : public virtual basic_user_base<_Jet>
  {
  public:
    //   public types
    typedef double  point_type;
    typedef _Weight weight_type;
    typedef _Sample sample_type;
    typedef _Traits sample_traits;
    typedef _PointTraits point_traits;
    typedef distbook<_Sample, double, _Traits, _PointTraits> distbook_type;
    typedef _Conv<_Weight, _Sample> conversion_type;
    
    //   destructor
    virtual ~basic_user() {}
    
  protected:
    //   do any operetaions at the end of the event
    void operations_at_the_end_of_event() {
      ++_M_dist; 
    }
    
    //   save the results in txt and binary mode
    void write_result_txt(std::ostream& os) {
      os<<_M_dist<<"\n";
    }
    
    void write_result_bin(std::ostream& os) {
      write(os, _M_dist);
    }
   
    void write_typeinfo(std::ostream& os) {
      os<<typeid(distbook_type).name()<<"\n";
    }
  
  private:
    //    data members
    distbook<_Sample, double, _Traits, _PointTraits> _M_dist;   

  protected:
    //   It is visible for the inhereted class. One can have some advanced 
    //   conversion that requires special initialization.  
    static conversion_type _S_conv;
    
  protected:
    //   create an one dimensional distribution
    void phys(int, const char *, unsigned int, double, double);
    void phys(int id, const char *name, unsigned int nbin, const double *base) {   
      _M_dist.create(id, name, nbin, base);
    }
      
    void physfillh(int id, double x, const _Weight& w) {
      _M_dist.accumulate(id, __helper_basic_user::_G_step_dbl, x, _S_conv(w));
    }
    
  private:
    template<class _Ptr>
    struct unary_func_conv 
    {
      typedef typename _Ptr::argument_type argument_type;
      typedef typename _Ptr::result_type result_type;
      
      unary_func_conv(const _Ptr& f)
		: _M_func(f) {}
      
      _Sample operator()(argument_type x) const {
		return _S_conv(_M_func(x));
      }
      
      const _Ptr& _M_func;
    };

  protected:
    template<class _Arg, class _Res> 
    void physfill(int id, const std::pointer_to_unary_function<_Arg, _Res>& f) {
      typedef unary_func_conv<std::pointer_to_unary_function<_Arg, _Res> > func_conv_type;
      _M_dist.accumulate(id, ptr_fun(func_conv_type(f), &func_conv_type::operator()));
    }
    
    template<class _Sp, class _Arg, class _Res> 
    void physfill(int id, const pointer_to_unary_member_function<_Sp, _Arg, _Res>& f) {
      typedef unary_func_conv<pointer_to_unary_member_function<_Sp, _Arg, _Res> > func_conv_type;
      _M_dist.accumulate(id, ptr_fun(func_conv_type(f), &func_conv_type::operator()));
    }
    
    template<class _Sp, class _Arg, class _Res> 
    void physfill(int id, const pointer_to_unary_const_member_function<_Sp, _Arg, _Res>& f) {
      typedef unary_func_conv<pointer_to_unary_const_member_function<_Sp, _Arg, _Res> > func_conv_type;
      _M_dist.accumulate(id, ptr_fun(func_conv_type(f), &func_conv_type::operator()));
    }
   
 
    template<class _Arg> 
    void physfill(int id, const std::pointer_to_unary_function<_Arg, _Sample>& f) {
      _M_dist.accumulate(id, f);
    }
    
    template<class _Sp, class _Arg> 
    void physfill(int id, const pointer_to_unary_member_function<_Sp, _Arg, _Sample>& f) {
      _M_dist.accumulate(id, f);
    }
    
    template<class _Sp, class _Arg> 
    void physfill(int id, const pointer_to_unary_const_member_function<_Sp, _Arg, _Sample>& f) {
      _M_dist.accumulate(id, f);
    }
    
    template<class _Arg1, class _Arg2, class _Res> 
    void physfill(int id, const std::pointer_to_binary_function<_Arg1, _Arg2, _Res>& f, _Arg1 x, const _Weight& w) {
      _M_dist.accumulate(id, f, x, _S_conv(w));
    }
    
    template<class _Sp, class _Arg1, class _Arg2, class _Res> 
    void physfill(int id, const pointer_to_binary_member_function<_Sp, _Arg1, _Arg2, _Res>& f, _Arg1 x, const _Weight& w) {
      _M_dist.accumulate(id, f, x, _S_conv(w));
    }
    
    template<class _Sp, class _Arg1, class _Arg2, class _Res> 
    void physfill(int id, const pointer_to_binary_const_member_function<_Sp, _Arg1, _Arg2, _Res>& f, _Arg1 x, const _Weight& w) {
      _M_dist.accumulate(id, f, x, _S_conv(w));
    }
  };

  template<class _Jet, typename _Weight, typename _Sample, template<class _Xp, class _Yp> class _Conv, class _Traits, class _PointTraits>
  _Conv<_Weight, _Sample> basic_user<_Jet, double, _Weight, _Sample, _Conv, _Traits, _PointTraits>::_S_conv;
  
  template<class _Jet, typename _Weight, typename _Sample, template<class _Xp, class _Yp>  class _Conv, class _Traits, class _PointTraits> 
  void basic_user<_Jet, double, _Weight, _Sample, _Conv, _Traits, _PointTraits>::
  phys(int id, const char *name, unsigned int npoints, double min, double max)
  {
    double *p = new double[npoints+1];
    double dx = (max-min)/npoints;
    
    for(unsigned int i = 0; i <= npoints; i++) 
      p[i] = min + i*dx; 
    
    _M_dist.create(id, name, npoints+1, p);
    
    delete [] p;
  }
 


  //  Specialization:
  //         One dimensional distribution with histogram support
  //                                 (_Point = histpoint1d)

  //   Helper functions to compute Dirac delta distribution and step function distribution
  namespace __helper_basic_user {
    extern std::pointer_to_binary_function<double, const histpoint1d&, double> _G_dirac;
    extern std::pointer_to_binary_function<double, const histpoint1d&, double> _G_step;
  }  //  namespace __helper_basic_user
 
  
  template<class _Jet, typename _Weight, typename _Sample, template<class _Xp, class _Yp> class _Conv, class _Traits, class _PointTraits>
  class basic_user<_Jet, histpoint1d, _Weight, _Sample, _Conv, _Traits, _PointTraits> 
    : public virtual basic_user_base<_Jet>
  {
  public:
    //   public types
    typedef histpoint1d point_type;
    typedef _Weight     weight_type;
    typedef _Sample     sample_type;
    typedef _Traits     sample_traits;
    typedef _PointTraits point_traits;
    typedef distbook<_Sample, histpoint1d, _Traits, _PointTraits> distbook_type;
    typedef _Conv<_Weight, _Sample> conversion_type;
    
    //   destructor
    virtual ~basic_user() {}
    
  protected:
    //   do any operetaions at the end of the event
    void operations_at_the_end_of_event() {
      _M_sdirac.accumulate();
      ++_M_dist; 
    }
    
    //   save the results in txt and binary mode
    void write_result_txt(std::ostream& os) {
      os<<_M_dist<<"\n";
    }
    
    void write_result_bin(std::ostream& os) {
      write(os, _M_dist);
    }
    
    void write_typeinfo(std::ostream& os) {
      os<<typeid(distbook_type).name()<<"\n";
    }
    
  private:
    //    data members
    distbook<_Sample, histpoint1d, _Traits> _M_dist;   
    
  protected:
    //   It is visible for the inhereted class. One can have some advanced 
    //   conversion that requires special initialization.  
    static conversion_type _S_conv;
    
    //   Help to compute a better dirac delta distribution
    //    (I hope!)    
    struct sphysreal {
      explicit sphysreal(const _Weight& wt) : weight(wt) {} 
      const _Weight& weight;
    };
    
    struct sphyssub {
      explicit sphyssub(const _Weight& wt) : weight(wt) {} 
      const _Weight& weight;
    };
    
  private:
    //   helper class to implement the smoother Dirac-delta
    class _PhysDirac
    {
    private:
      //  collect the real and dipole contributions
      struct smooth {
	explicit smooth(double delta) : dy(delta), Ny(0U) {}
	
	double dy, yr, ys;
	unsigned int Ny;
	_Sample weight;
      };
      
      //   constructor
      explicit _PhysDirac(distbook_type *dbook) 
	: _M_dbook(dbook) {}
      
      void create(int id, double delta) {
	_M_smooth.insert(std::pair<int, smooth>(id, smooth(delta)));
      }
      
      //  fill the weight at the and of the event
      void accumulate();
      
      //   these operators fill the weights
      void operator() (int id, double x, const _Weight& wt) {
	_M_dbook -> accumulate(id, __helper_basic_user::_G_dirac,
			       x, _S_conv(wt));
      }
      
      void operator() (int id, double x, const sphysreal& wt) 
      {
	smooth& smth = _M_smooth.find(id) -> second;
	smth.Ny = 1U; 
	smth.yr = x; smth.ys = 0.0;
	_Traits::assign(smth.weight, _S_conv(wt.weight));
      }
      
      void operator() (int id, double x, const sphyssub& wt) 
      {
	smooth& smth = _M_smooth.find(id) -> second;
	
	if(smth.Ny != 0U)
	  if(std::abs(smth.yr - x) < smth.dy) {
	    smth.Ny++; smth.ys += x; 
	    _Traits::assadd(smth.weight, _S_conv(wt.weight));
	    return;
	  }
        
	_M_dbook -> accumulate(id, __helper_basic_user::_G_dirac, 
			       x, _S_conv(wt.weight));
      }
      
      //   data members
      distbook_type *_M_dbook;
      std::map<int, smooth> _M_smooth;
      
      //   friend declarations
      friend class basic_user;
    };
    
  protected:
    //  constructors and destructor
    basic_user() : _M_sdirac(&_M_dist) {}
    
    //   create an one dimensional distribution
    void phys(int, const char *, unsigned int, double, double, double = 0.1);
    void phys(int, const char *, unsigned int, const double *, double, double = 0.1);
    void phys(int, const char *, unsigned int, const point_type *, double = 0.1); 
    
    //   fill the one dimensional distributions (dirac delta and step function) 
    void physfilld(int id, double x, const _Weight& w) {
      _M_dist.accumulate(id, __helper_basic_user::_G_dirac, x, _S_conv(w));
    }
    
    void physfillh(int id, double x, const _Weight& w) {
      _M_dist.accumulate(id, __helper_basic_user::_G_step, x, _S_conv(w));
    }
    
    void sphysfilld(int id, double x, const _Weight& wt) {
      _M_sdirac(id, x, _S_conv(wt));
    }
    
    void sphysfilld(int id, double x, const sphysreal& wt) {
      _M_sdirac(id, x, wt);
    }
    
    void sphysfilld(int id, double x, const sphyssub& wt) {
      _M_sdirac(id, x, wt);
    }
    
  private:
    template<class _Ptr>
    struct unary_func_conv 
    {
      typedef typename _Ptr::argument_type argument_type;
      typedef typename _Ptr::result_type result_type;
      
      unary_func_conv(const _Ptr& f) : _M_func(f) {}
      
      _Sample operator()(argument_type x) const {
	return _S_conv(_M_func(x));
      }
      
      const _Ptr& _M_func;
    };
    
  protected:
    template<class _Arg, class _Res> 
    void physfill(int id, const std::pointer_to_unary_function<_Arg, _Res>& f) {
      typedef unary_func_conv<std::pointer_to_unary_function<_Arg, _Res> > func_conv_type;
      _M_dist.accumulate(id, ptr_fun(func_conv_type(f), &func_conv_type::operator()));
    }
    
    template<class _Sp, class _Arg, class _Res> 
    void physfill(int id, const pointer_to_unary_member_function<_Sp, _Arg, _Res>& f) {
      typedef unary_func_conv<pointer_to_unary_member_function<_Sp, _Arg, _Res> > func_conv_type;
      _M_dist.accumulate(id, ptr_fun(func_conv_type(f), &func_conv_type::operator()));
    }
    
    template<class _Sp, class _Arg, class _Res> 
    void physfill(int id, const pointer_to_unary_const_member_function<_Sp, _Arg, _Res>& f) {
      typedef unary_func_conv<pointer_to_unary_const_member_function<_Sp, _Arg, _Res> > func_conv_type;
      _M_dist.accumulate(id, ptr_fun(func_conv_type(f), &func_conv_type::operator()));
    }
    
    
    template<class _Arg> 
    void physfill(int id, const std::pointer_to_unary_function<_Arg, _Sample>& f) {
      _M_dist.accumulate(id, f);
    }
    
    template<class _Sp, class _Arg> 
    void physfill(int id, const pointer_to_unary_member_function<_Sp, _Arg, _Sample>& f) {
      _M_dist.accumulate(id, f);
    }
    
    template<class _Sp, class _Arg> 
    void physfill(int id, const pointer_to_unary_const_member_function<_Sp, _Arg, _Sample>& f) {
      _M_dist.accumulate(id, f);
    }
    
    template<class _Arg1, class _Arg2, class _Res> 
    void physfill(int id, const std::pointer_to_binary_function<_Arg1, _Arg2, _Res>& f, _Arg1 x, const _Weight& w) {
      _M_dist.accumulate(id, f, x, _S_conv(w));
    }
    
    template<class _Sp, class _Arg1, class _Arg2, class _Res> 
    void physfill(int id, const pointer_to_binary_member_function<_Sp, _Arg1, _Arg2, _Res>& f, _Arg1 x, const _Weight& w) {
      _M_dist.accumulate(id, f, x, _S_conv(w));
    }
    
    template<class _Sp, class _Arg1, class _Arg2, class _Res> 
    void physfill(int id, const pointer_to_binary_const_member_function<_Sp, _Arg1, _Arg2, _Res>& f, _Arg1 x, const _Weight& w) {
      _M_dist.accumulate(id, f, x, _S_conv(w));
    }
    
  private:
    _PhysDirac _M_sdirac;
  };
  
  template<class _Jet, typename _Weight, typename _Sample,template<class _Xp, class _Yp> class _Conv, class _Traits, class _PointTraits>
  _Conv<_Weight, _Sample> basic_user<_Jet, histpoint1d, _Weight, _Sample, _Conv, _Traits, _PointTraits>::_S_conv;
  	
  template<class _Jet, typename _Weight, typename _Sample, template<class _Xp, class _Yp>  class _Conv, class _Traits, class _PointTraits> 
  void basic_user<_Jet, histpoint1d, _Weight, _Sample, _Conv, _Traits, _PointTraits>::_PhysDirac::accumulate() 
  {
    typename std::map<int, smooth>::iterator iter = _M_smooth.begin();
    unsigned int ny;
    double ya;
    
    while(iter != _M_smooth.end()) {
      ny = iter -> second.Ny;
      smooth& smth = iter -> second;
      
      if(ny != 0U) { 
	ya = (ny == 1U ? smth.yr : 0.5*(smth.yr + smth.ys/(ny-1U)));
	_M_dbook -> accumulate(iter -> first, __helper_basic_user::_G_dirac, ya, smth.weight);
	iter -> second.Ny = 0U;
      }
      ++iter;
    }
  }
  
  template<class _Jet, typename _Weight, typename _Sample, template<class _Xp, class _Yp>  class _Conv, class _Traits, class _PointTraits> 
  void basic_user<_Jet, histpoint1d, _Weight, _Sample, _Conv, _Traits, _PointTraits>::
  phys(int id, const char *name, unsigned int nbin, double min, double max, double eps)
  {
    point_type *p = new point_type[nbin];
    double dx = 0.5*(max-min)/nbin;
    
    p[0].xmin = min;      
    for(unsigned int i = 0; i < nbin-1; i++) {
      p[i].xmid = p[i].xmin + dx;
      p[i].xmax = p[i+1].xmin = p[i].xmid + dx;
    }
    
    p[nbin-1].xmid = p[nbin-1].xmin + dx;
    p[nbin-1].xmax = p[nbin-1].xmid + dx;
    
    _M_dist.create(id, name, nbin, p);
    _M_sdirac.create(id, 2.0*dx*eps);
    
    delete [] p;
  }
  
  template<class _Jet, typename _Weight, typename _Sample, template<class _Xp, class _Yp>  class _Conv, class _Traits, class _PointTraits> 
  void basic_user<_Jet, histpoint1d, _Weight, _Sample, _Conv, _Traits, _PointTraits>::
  phys(int id, const char *name, unsigned int nbin, const double *base, double round, double eps) 
  {
    point_type *p = new point_type[nbin];
    
    for(unsigned int i = 0; i < nbin; i++) {
      p[i].xmin = base[i] - round;
      p[i].xmid = base[i];
      p[i].xmax = base[i] + round;
    }
    
    _M_dist.create(id, name, nbin, p);
    _M_sdirac.create(id, round*eps);
    
    delete [] p;
  }
  
  template<class _Jet, typename _Weight, typename _Sample, template<class _Xp, class _Yp> class _Conv, class _Traits, class _PointTraits> 
  void basic_user<_Jet, histpoint1d, _Weight, _Sample, _Conv, _Traits, _PointTraits>::
  phys(int id, const char *name, unsigned int nbin, const histpoint1d *base, double eps) 
  {
    double tmp, round = base[0].xmax - base[0].xmin;
    for(unsigned i = 1; i < nbin; i++) {
      tmp = base[i].xmax - base[i].xmin;
      if(tmp < round) round = tmp;
    }
    
    _M_dist.create(id, name, nbin, base);
    _M_sdirac.create(id, round*eps);
  }


  //  Specialization:
  //         Two dimensional distribution with histogram support
  //                                 (_Point = histpoint2d)
  
  //   Helper functions to compute Dirac delta distribution and step function distribution
  namespace __helper_basic_user {
   
    struct __helper_point2d {
      __helper_point2d(const double& __x, const double& __y) : x(__x), y(__y) {}
      double x,y;
    };

    extern std::pointer_to_binary_function<__helper_point2d, const histpoint2d&, double> _G_hist2d_dd;
    extern std::pointer_to_binary_function<__helper_point2d, const histpoint2d&, double> _G_hist2d_dh;
    extern std::pointer_to_binary_function<__helper_point2d, const histpoint2d&, double> _G_hist2d_hd;
    extern std::pointer_to_binary_function<__helper_point2d, const histpoint2d&, double> _G_hist2d_hh;
  }  //  namespace __helper_basic_user
 

  
  template<class _Jet, typename _Weight, typename _Sample, template<class _Xp, class _Yp> class _Conv, class _Traits, class _PointTraits>
  class basic_user<_Jet, histpoint2d, _Weight, _Sample, _Conv, _Traits, _PointTraits> 
    : public virtual basic_user_base<_Jet>
  {
  public:
    //   public types
    typedef histpoint2d point_type;
    typedef _Weight     weight_type;
    typedef _Sample     sample_type;
    typedef _Traits     sample_traits;
    typedef _PointTraits point_traits;
    typedef distbook<_Sample, histpoint2d, _Traits, _PointTraits> distbook_type;
    typedef _Conv<_Weight, _Sample> conversion_type;
    
    //   destructor
    virtual ~basic_user() {}
    
  protected:
    //   do any operetaions at the end of the event
    void operations_at_the_end_of_event() {
      ++_M_dist; 
    }
    
    //   save the results in txt and binary mode
    void write_result_txt(std::ostream& os) {
      os<<_M_dist<<"\n";
    }
    
    void write_result_bin(std::ostream& os) {
      write(os, _M_dist);
    }
   
    void write_typeinfo(std::ostream& os) {
      os<<typeid(distbook_type).name()<<"\n";
    }
    
  private:
    //    data members
    distbook<_Sample, histpoint2d, _Traits, _PointTraits> _M_dist;   

  protected:
    //   It is visible for the inhereted class. One can have some advanced 
    //   conversion that requires special initialization.  
    static conversion_type _S_conv;
        
    //   create an two dimensional distribution
    void phys(int, const char *, unsigned int, double, double, unsigned int, double, double);
    void phys(int, const char *, unsigned int, const double *, double, unsigned int, const double *, double);
    void phys(int id, const char *name, unsigned int nbin, const histpoint2d *base) {
      _M_dist.create(id, name, nbin, base);
    } 
    
    
    //   fill the two dimensional distributions (dirac delta and step function) 
    void physfilldd(int id, double x, double y, const _Weight& w) {
      _M_dist.accumulate(id, __helper_basic_user::_G_hist2d_dd, __helper_basic_user::__helper_point2d(x,y), _S_conv(w));
    }
    
    void physfilldh(int id, double x, double y, const _Weight& w) {
      _M_dist.accumulate(id, __helper_basic_user::_G_hist2d_dh, __helper_basic_user::__helper_point2d(x,y), _S_conv(w));
    }
    
    void physfillhd(int id, double x, double y, const _Weight& w) {
      _M_dist.accumulate(id, __helper_basic_user::_G_hist2d_hd, __helper_basic_user::__helper_point2d(x,y), _S_conv(w));
    }
    
    void physfillhh(int id, double x, double y, const _Weight& w) {
      _M_dist.accumulate(id, __helper_basic_user::_G_hist2d_hh, __helper_basic_user::__helper_point2d(x,y), _S_conv(w));
    }

  private:
    template<class _Ptr>
    struct unary_func_conv 
    {
      typedef typename _Ptr::argument_type argument_type;
      typedef typename _Ptr::result_type result_type;
      
      unary_func_conv(const _Ptr& f)
	: _M_func(f) {}
      
      _Sample operator()(argument_type x) const {
	return _S_conv(_M_func(x));
      }
 
      const _Ptr& _M_func;
    };

  protected:
    //   fill the distribution
    template<class _Arg, class _Res> 
    void physfill(int id, const std::pointer_to_unary_function<_Arg, _Res>& f) {
      typedef unary_func_conv<std::pointer_to_unary_function<_Arg, _Res> > func_conv_type;
      _M_dist.accumulate(id, ptr_fun(func_conv_type(f), &func_conv_type::operator()));
    }
    
    template<class _Sp, class _Arg, class _Res> 
    void physfill(int id, const pointer_to_unary_member_function<_Sp, _Arg, _Res>& f) {
      typedef unary_func_conv<pointer_to_unary_member_function<_Sp, _Arg, _Res> > func_conv_type;
      _M_dist.accumulate(id, ptr_fun(func_conv_type(f), &func_conv_type::operator()));
    }
    
    template<class _Sp, class _Arg, class _Res> 
    void physfill(int id, const pointer_to_unary_const_member_function<_Sp, _Arg, _Res>& f) {
      typedef unary_func_conv<pointer_to_unary_const_member_function<_Sp, _Arg, _Res> > func_conv_type;
      _M_dist.accumulate(id, ptr_fun(func_conv_type(f), &func_conv_type::operator()));
    }
   
    template<class _Arg> 
    void physfill(int id, const std::pointer_to_unary_function<_Arg, _Sample>& f) {
      _M_dist.accumulate(id, f);
    }
    
    template<class _Sp, class _Arg> 
    void physfill(int id, const pointer_to_unary_member_function<_Sp, _Arg, _Sample>& f) {
      _M_dist.accumulate(id, f);
    }
    
    template<class _Sp, class _Arg> 
    void physfill(int id, const pointer_to_unary_const_member_function<_Sp, _Arg, _Sample>& f) {
      _M_dist.accumulate(id, f);
    }
    
    template<class _Arg1, class _Arg2, class _Res> 
    void physfill(int id, const std::pointer_to_binary_function<_Arg1, _Arg2, _Res>& f, _Arg1 x, const _Weight& w) {
      _M_dist.accumulate(id, f, x, _S_conv(w));
    }
    
    template<class _Sp, class _Arg1, class _Arg2, class _Res> 
    void physfill(int id, const pointer_to_binary_member_function<_Sp, _Arg1, _Arg2, _Res>& f, _Arg1 x, const _Weight& w) {
      _M_dist.accumulate(id, f, x, _S_conv(w));
    }
    
    template<class _Sp, class _Arg1, class _Arg2, class _Res> 
    void physfill(int id, const pointer_to_binary_const_member_function<_Sp, _Arg1, _Arg2, _Res>& f, _Arg1 x, const _Weight& w) {
      _M_dist.accumulate(id, f, x, _S_conv(w));
    }
  };

	
  template<class _Jet, typename _Weight, typename _Sample, template<class _Xp, class _Yp> class _Conv, class _Traits, class _PointTraits>
  _Conv<_Weight, _Sample> basic_user<_Jet, histpoint2d, _Weight, _Sample, _Conv, _Traits, _PointTraits>::_S_conv;
  
  
  template<class _Jet, typename _Weight, typename _Sample, template<class _Xp, class _Yp>  class _Conv, class _Traits, class _PointTraits> 
  void basic_user<_Jet, histpoint2d, _Weight, _Sample, _Conv, _Traits, _PointTraits>::
  phys(int id, const char *name, unsigned int nx, double xmin, double xmax, unsigned int ny, double ymin, double ymax)
  {
    point_type *p = new point_type[nx*ny];
    double xn,xd, xx, yn, yd, yx, dx = 0.5*(xmax-xmin)/nx, dy = 0.5*(ymax-ymin)/ny;
    unsigned int idx;
    
    for(unsigned int i = 0; i < nx; i++) {
      xn = xmin + i*2.0*dx; xd = xn + dx; xx = xd + dx;
      for(unsigned int j = 0; j < ny; j++) {
	idx = ny*i+j;
	yn = ymin + j*2.0*dy; yd = yn + dy; yx = yd + dy;
	p[idx].xmin = xn; p[idx].xmid = xd; p[idx].xmax = xx;
	p[idx].ymin = yn; p[idx].ymid = yd; p[idx].ymax = yx;
      }
    }
    
    _M_dist.create(id, name, nx*ny, p);
    
    delete [] p;
  }
  
  template<class _Jet, typename _Weight, typename _Sample, template<class _Xp, class _Yp>  class _Conv, class _Traits, class _PointTraits> 
  void basic_user<_Jet, histpoint2d, _Weight, _Sample, _Conv, _Traits, _PointTraits>::
  phys(int id, const char *name, unsigned int nx, const double *xbase, double rx, unsigned int ny, const double *ybase, double ry) 
  {
    point_type *p = new point_type[nx*ny];
    unsigned int idx;
    
    for(unsigned int i = 0; i < nx; i++)
      for(unsigned int j = 0; j < ny; j++) {
	idx = ny*i+j;
	p[idx].xmin = xbase[i]-rx; 
	p[idx].xmid = xbase[i]; 
	p[idx].xmax = xbase[i]+rx;
	
	p[idx].ymin = ybase[j]-ry; 
	p[idx].ymid = ybase[j]; 
	p[idx].ymax = ybase[j]+ry;
      }
    
    _M_dist.create(id, name, nx*ny, p);
    
    delete [] p;
  }


}  //   namespace nlo

#endif
