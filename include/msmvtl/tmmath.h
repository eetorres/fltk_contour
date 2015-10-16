//========================================================================
// Mini Scientific Matrix and Vector Template Library - MSMVTL
// MSMVTL is a simple  matrix and vector template library
//
// FILE: tmmath.h
//       math operation for the matrix template
//
// Copyrigth 2002 by Edmanuel Torres, eetorres@gmail.com
//
// Don't hesitate to contact me for any question, suggestion,
// modification or bug. Get the latest version at:
// http://msmvtl.sourceforge.net
//
// This library is free  software;  you  can  redistribute  it and/or
// modify it  under  the  terms  of  the   GNU Library General Public
// License  as  published  by  the  Free  Software Foundation; either
// version 2 of the License,  or  (at your option)  any later version.
//
// This  library  is  distributed  in the hope that it will be useful,
// but  WITHOUT ANY WARRANTY;  without  even  the  implied warranty of
// MERCHANTABILITY  or FITNESS FOR A PARTICULAR PURPOSE.   See the GNU
// Library General Public License for more details.
//
// You should have  received a copy  of the GNU Library General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA.
//
//========================================================================


#ifndef _TMMATH_H_
#define _TMMATH_H_ 1

#include <msmvtl/tvmath.h>

template<class _T> inline TVector<real> mtrace(_T _m){
    TVector<real> _tv(_m.rows());
    if(_m.is_square()){
    	for(unsigned int i=0; i<_tv.size(); i++)
	_tv[i] = _m[i][i];
    }
    return _tv;
}

template<class _T> inline _T mnormalize( _T _m, _T _p){
    real min_scale, max_scale, _k;
    TVector<real> _t;
    _T _tm(_m.rows(),0);
    for(unsigned int i=0; i<_m.cols(); i++){
        min_scale = _p[i][0];
        max_scale = _p[i][1];
	_t = _m.get_col(i);
	if(-min_scale != max_scale){
	    _k = 0.5*(max_scale+min_scale);
	    _t = (_t - _k);
	    _t = (1.0/_k)*_t;
	}else
	    _t = (1.0/max_scale)*_t;
	_tm.add_col(_t);
    }
    return _tm;
}

template<class _T> inline _T minterval(_T _m){
    _T _intv(_m.cols(),2);
    TVector<real> _t;
    for(unsigned int i=0; i<_m.cols(); i++){
	_t = _m.get_col(i);
	_intv[i][0] =_t.get_min();
	_intv[i][1] =_t.get_max();
    }
    return _intv;
}

template<class _T> _T mnormrange(_T _m1, _T _m2){
    TVector<real> _t;
    _T _m(_m1.rows(),0);
    for(unsigned int i=0; i<_m1.cols(); i++){
	_t = _m1.get_col(i);
	_t = vnormrange(_t,_m2[i]);
	_m.add_col(_t);
    }
    return _m;
}

template<class _T> _T mrange(_T _m1, _T _m2){
    TVector<real> _t;
    _T _m(_m1.rows(),0);
    for(unsigned int i=0; i<_m1.cols(); i++){
	_t = _m1.get_col(i);
	_t = vrange(_t,_m2[i]);
	_m.add_col(_t);
    }
    return _m;
}

template<class _T> _T mscalexy(_T _m, _T _ml, real _ss){
    unsigned int i;
    TVector<double> _t(_m.rows());
    _T _tm(_m.rows(),0);
    for(i=0; i<_m.cols()-1; i++){
	_t = _m.get_col(i);
	_t = vscale(_t,0.5*(_ml[i][0]+_ml[i][1]),_ss);
	_tm.add_col(_t);
    }
    _t = _m.get_col(i);
    _t = vscale(_t-0.5,0.5*(_ml[i][0]+_ml[i][1]),_ml[i][1]-_ml[i][0]);
    _tm.add_col(_t);
    return _tm;
}

template<class _T> _T mnormxy(_T _m, _T _ml, real _ss){
    unsigned int i;
    TVector<double> _t(_m.rows());
    _T _tm(_m.rows(),0);
    for(i=0; i<_m.cols()-1; i++){
	_t = _m.get_col(i);
	_t = vnormalize(_t,0.5*(_ml[i][0]+_ml[i][1]),_ss);
	_tm.add_col(_t);
    }
    _t = _m.get_col(i);
    _t = 0.5+vnormalize(_t,0.5*(_ml[i][0]+_ml[i][1]),(_ml[i][1]-_ml[i][0]));
    _tm.add_col(_t);
    return _tm;
}

template<class _T> _T mNormXYmat(_T _m){
    unsigned int i;
    TVector<double> _t(_m.rows());
    _T _tm(_m.rows(),0);
    _T _mlim;
    real _ss1, _ss2;
    
    _mlim = minterval(_m); // getting scaling parameters
    _mlim.redim(2,2);
    _ss1 = (_mlim[0][1] - _mlim[0][0]);
    _ss2 = (_mlim[1][1] - _mlim[1][0]);
    real _ss = maxi(_ss1,_ss2);
    for(i=0; i<_m.cols()-1; i++){
	_t = _m.get_col(i);
	_t = vnormalize(_t,0.5*(_mlim[i][0]+_mlim[i][1]),_ss);
	_tm.add_col(_t);
    }
    _t = _m.get_col(i);
    _t = vnorm(_t);
    _tm.add_col(_t);
    return _tm;
}

#endif

//

