// for numerical programming.
// 
// tmmanip.h
//
// Copyright 2002 by Edmanuel Torres A. (eetorres@gmail.com)
//
// This widget is  free  software;  you  can  redistribute  it and/or
// modify it  under  the  terms  of  the   GNU Library General Public
// License  as  published  by  the  Free  Software Foundation; either
// version 2 of the License,  or  (at your option)  any later version
// or much better FLTK license, which allow you static linking.
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
// Sent me suggestion, modification or bugs. Don't hesitate to contact
// me for any question, I will be very grateful with your feedbacks.
// Get the last version at http://eetorres.googlepages.com/fltk_en


#ifndef _TMMANIP_H_
#define _TMMANIP_H_ 1

#include <msmvtl/tvmath.h>

template<class _T> inline _T m_get_intervals(_T _m){
    _T _intv(_m.cols(),2);
    TVector<real> _t;
    for(unsigned int i=0; i<_m.cols(); i++){
	_t = _m.get_col(i);
	_intv[i][0] =_t.get_min();
	_intv[i][1] =_t.get_max();
    }
    return _intv;
}

template<class _T> _T m_normalize_ranges(_T _m1, _T _m2){
    TVector<real> _t;
    _T _m(_m1.rows(),0);
    for(unsigned int i=0; i<_m1.cols(); i++){
	_t = _m1.get_col(i);
	_t = vnormrange(_t,_m2[i]);
	_m.add_col(_t);
    }
    return _m;
}

template<class _T> _T m_rescale(_T _m1, _T _m2){
    TVector<real> _t;
    _T _m(_m1.rows(),0);
    for(unsigned int i=0; i<_m1.cols(); i++){
	_t = _m1.get_col(i);
	_t = vrange(_t,_m2[i]);
	_m.add_col(_t);
    }
    return _m;
}

#endif

//
    
