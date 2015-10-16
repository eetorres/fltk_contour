//========================================================================
// Mini Scientific Matrix and Vector Template Library - MSMVTL
// MSMVTL is a simple  matrix and vector template library
//
// FILE: tvmath.h
//       math operation for the vector template
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


#ifndef _TVMATH_H_
#define _TVMATH_H_ 1

#include <msmvtl/tvector.h>

//////////////////////////////////////////////////////////////////////////////

template<class _T, class _Tu> _T vscale(_T _v, _Tu  _shift, _Tu _scale){
    _T _t(_v.size());
    _t = ((_v * _scale)+_shift);
    return _t;
}

template<class _T, class _Tu> _T vnormalize(_T _v, _Tu  _shift, _Tu _scale){
    _T _t(_v.size());
    _t = (_v - _shift)/_scale;
    return _t;
}

template<class _T, class _Tu> _T vrange(_T _v, _Tu _min_v, _Tu _max_v){
    _T _t;
    double min_scale=_v.get_min(), max_scale=_v.get_max();
    _t = (_v - min_scale)/(max_scale - min_scale);
    _t = _t*(_max_v-_min_v);
    _t = _t+_min_v;
    return _t;
}

template<class _T> _T vnorm(_T _v){
    _T _t(_v.size());
    double min_scale=_v.get_min(), max_scale=_v.get_max();
    _t = (_v - min_scale)/(max_scale - min_scale);
    return _t;
}

template<class _T> _T vnormrange(_T _v, _T _v2){
    _T _t=_v;
    real _p;
    _p = (0.5*(_v2[1]+_v2[0]));
    _t = (_t - _p);
    _p = (_v2[1]-_v2[0]);
    _t = (_t/_p);
    return _t;
}

template<class _T> _T vrange(_T _v, _T _v2){
    _T _t=_v;
    real _p;
    _p = (_v2[1]-_v2[0]);
    _t = (_t*_p);
    _p = (0.5*(_v2[1]+_v2[0]));
    _t = (_t + _p);
    return _t;
}

template<class _T> _T vcentre( _T _v, _T _intv){
    real min_scale = _intv[0];
    real max_scale = _intv[1];
    _T _t;
    _t = (_v - 0.5*(min_scale+max_scale))/(0.5*(max_scale + min_scale));
    return _t;
}

template<class _T> _T acos(_T _v){
    _T _t(_v.size());
    for(unsigned int i=0; i<_v.size(); i++)
    	_t[i] = acos(_v[i]);
    return _t;
}

#endif

//
