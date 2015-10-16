//========================================================================
// Mini Scientific Matrix and Vector Template Library - MSMVTL
// MSMVTL is a simple  matrix and vector template library
//
// FILE: libdef.h
//       MSMVTL definitions
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


#ifndef _DEF_H_
#define _DEF_H_ 1

#include <float.h>
#include <math.h>
#include <fstream>
#include <iostream>
#include <string.h>
#include <string>

#define __CHECK_BOUNDS
#define __GET_MESSAGE


#ifndef NUM
#define NUM	0x023
#endif

#ifndef TIL
#define TIL	0x07E
#endif

#ifndef NLN
#define NLN	0x00A
#endif

#ifndef SPC
#define SPC	0x020
#endif

#ifndef TAB
#define TAB	0x009
#endif

#if defined(WIN32) || defined(WIN64)
#define POW(x,y)	powf(x,y);
#define ABS(x)		abs(x);
#else
#define POW(x,y)	pow(x,y);
#define ABS(x)		fabs(x);
#endif

#ifndef mini
    #define mini(a,b) (((a)<(b))?(a):(b))
#endif

#ifndef maxi
    #define maxi(a,b) (((a)>(b))?(a):(b))
#endif

#ifndef real
typedef double real;
#endif

#ifndef lreal
typedef long double lreal;
#endif

#ifndef uint
typedef unsigned int uint;
#endif

enum runtime_error { END=0, EXIT=-1 };

#ifndef real
typedef struct {
  real r,g,b;
} rgb;
#endif

//using namespace std;

#endif

