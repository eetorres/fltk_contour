// Latest modification 20/08/2006
// ========================================================================
// version: 0.7
// name:    geometry.h
//
// Copyrigth 2002 by Edmanuel Torres A. (eetorres@gmail.com)
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
// Sent me suggestions, modifications or bugs. Don't hesitate to contact
// me for any question, I will be very grateful with your feedbacks.
// Get the last version at http://eetorres.googlepages.com/fltk_en
//
//========================================================================


#ifndef _GEOMETRY_H_
#define _GEOMETRY_H_


#include <float.h>
#include <vector>
#include <iostream>

using namespace std;

#define GM_EPSILON DBL_EPSILON

#ifndef GM_TRUE
#define GM_TRUE true
#endif

#ifndef GM_FALSE
#define GM_FALSE false
#endif

#ifndef GM_RIGHT
#define GM_RIGHT true
#endif

#ifndef GM_LEFT
#define GM_LEFT false
#endif

#ifndef GM_3D
#define GM_3D true
#endif

#ifndef GM_2D
#define GM_2D false
#endif


typedef double 	gm_real;
typedef short  	gm_short;
typedef int    	gm_int;
typedef unsigned int gm_uint;
typedef bool 	gm_bool;

template <class T, gm_uint pdim> class Geom_Vertex{

private:
    gm_uint _psize;
    vector<T>  pdata;
	
public:
    Geom_Vertex<T,pdim> (){
        _psize = pdim;
        pdata.resize(pdim);
    }
	
    ~Geom_Vertex<T,pdim> (){};

    void point2d(T x0, T y0){x() = x0; y() = y0;}
    void point3d(T x0, T y0, T z0){x()=x0; y()=y0; z()=z0;}
    void pointnd(Geom_Vertex<T,pdim>& xyz){ x() = xyz.x(); y() = xyz.y(); z() = xyz.z();}
    /*
    Geom_Vertex<T,pdim>& operator = (const Geom_Vertex<T,pdim>& p1){
        if (this!=&p1) {
	    _psize = p1._psize;
	    for(gm_uint i=0; i<_psize; i++)
		pdata[i] = p1.pdata[i];}
	return (*this);
    }*/
    T& operator[] (gm_uint i){
#ifdef __CHECK_BOUNDS
	if (i<0 || i>=pdim) cerr << "Data out of range " << i << EXIT;
#endif
	return (pdata[i]);
    };
    //T& operator[] (gm_uint i){
    Geom_Vertex& operator= (const Geom_Vertex& _v){
	pdata = _v.pdata;
	return *this;
    }
    
    vector<T> get_vec(void){
	return pdata;
    }
    
    inline T& _xyz(gm_uint d){ return pdata[d];};
    inline T& x(void){ return _xyz(0);};
    inline T& y(void){ return _xyz(1);};
    inline T& z(void){ return _xyz(2);};
    inline void x(T d){ _xyz(0)=d;};
    inline void y(T d){ _xyz(1)=d;};
    inline void z(T d){ _xyz(2)=d;};
    //friend class Geom_Triangle;
};

//============================================================//

template <class T, gm_uint vsize, gm_uint pdim> class Geom_Edges{

private:
    gm_uint _vsize, nv;
    vector< Geom_Vertex<T, pdim> > vdata;

public:

    Geom_Edges<T, vsize, pdim>(){
    	nv = 0; _vsize = vsize;
	vdata.resize(_vsize);
    }
    ~Geom_Edges<T, vsize, pdim>(){};
	
    void resize(gm_uint _vdim){
	_vsize = _vdim;
	vdata.resize(_vsize);}
    inline void null2d(gm_uint _vpos){
	vdata[_vpos].x() = 0; vdata[_vpos].y() = 0;}
    inline void null(gm_uint _vpos){
        for(gm_uint i=0; i<pdim;i++)
    	    vdata[_vpos][i]=0;
    }    
    inline void push_back(Geom_Vertex<T,pdim> p0){
    	_vsize++; nv = _vsize; vdata.push_back(p0);}
    inline void addlnpnt(Geom_Vertex<T,pdim> p0){
	vdata[nv] = p0; if(nv<(_vsize-1)) nv++;}
    inline void addlnpnt(Geom_Vertex<T,pdim> p0, gm_uint pos){
	vdata[pos] = p0; if(pos>nv) nv=pos;}
    inline void add3(gm_uint _p, Geom_Vertex<T,pdim> p0, Geom_Vertex<T,pdim> p1, Geom_Vertex<T,pdim> p2){
	vdata[_p] = p0; vdata[_p+1] = p1; vdata[_p+2] = p2; nv = _p + 3;}
    inline void reset(void){nv=0;}
    inline bool comp_xx(gm_uint e0, gm_uint e1){
      T data = (int)(vdata[e0].x() - vdata[e1].x());
      return abs(data)<= GM_EPSILON;
    }
    inline bool comp_yy(gm_uint e0,unsigned  int e1){
      T data = (int)(vdata[e0].y() - vdata[e1].y());
      return abs(data)<= GM_EPSILON;
    }
    inline bool comp_xy(int e0, int e1){
      T data = (int)(vdata[e0].x() - vdata[e1].y());
      return abs(data)<= GM_EPSILON;
    }
    inline bool nonull(int e){
	return(vdata[e].x() != 0 && vdata[e].y() != 0);}
    inline void size(int _dim){nv = _dim;}
    gm_uint size(void){return _vsize;}
    inline void addvertex(gm_uint n){nv = nv+n;}
    void delete_last(void){ vdata.pop_back(); nv--; _vsize--;};
    inline T& x(gm_uint p){return vdata[p].x();}
    inline T& y(gm_uint p){return vdata[p].y();}
    inline T& z(gm_uint p){return vdata[p].z();}
    
    Geom_Edges<T, vsize, pdim>& operator = (const Geom_Edges<T, vsize, pdim>& p1){
	if (this!=&p1){
	    nv = p1.nv; _vsize = p1._vsize;
    	    for(gm_uint i=0; i<_vsize; i++)
	        vdata[i] = p1.vdata[i];
	}
	return (*this);
    }
    
    Geom_Vertex<T,pdim>& operator[] (gm_uint r){
#ifdef __CHECK_BOUNDS
	if (r<0 || r>=_vsize)
	    cerr << "Vector fuera de rango " << r << EXIT;
#endif
        return (vdata[r]);
    }
};

//////////////////////////////////////////////////////////////////////

template <class T, gm_uint vdim> class Geom_Contour{
private:
	gm_uint csize;
	vector< Geom_Edges<T, 2, vdim> > cdata;

public:
	Geom_Contour(){
		csize = 0;
		cdata.resize(csize);
	}

	Geom_Contour<T, vdim>(gm_uint cs){
		csize = cs;
		cdata.resize(csize);
	}
	~Geom_Contour(){};

	void add(Geom_Edges<T, 2, vdim>& _v){
		cdata.push_back(_v);
		csize++;
	}
	gm_uint size(void){ 
	    return cdata.size();
	}
	void resize(gm_uint cs){
		csize = cs;
		cdata.resize(cs);
	}
	void clear(void){
		cdata.clear();
		csize = 0;
	}
	Geom_Edges<T, 2, vdim>& operator[] (gm_uint r) {
#ifdef __CHECK_BOUNDS
		if (r<0 || r>=csize)
			cerr << "Fuera de rango " << r << EXIT;
#endif
		return (cdata[r]);
	}
};

//============================================================//

template <class T> class Geom_Triangle {

private:
	Geom_Edges<T, 3, 3> tdata;

public:
    Geom_Triangle(){};
	
    ~Geom_Triangle(){};
	
    void data(T x0, T y0, T x1, T y1, T x2, T y2){
	x(0)=x0; x(1)=x1; x(2)=x2; y(0)=y0; y(1)=y1; y(2)=y2;
    }	
    inline void data(Geom_Vertex<T,3> p0, Geom_Vertex<T,3> p1, Geom_Vertex<T,3> p2){
	tdata[0] = p0; tdata[1] = p1; tdata[2] = p2;
    }
    inline bool collinear(void){
	if((fabs(y(0)-y(1)) <= GM_EPSILON) && (fabs(y(1)-y(2)) <= GM_EPSILON)){
	    return GM_TRUE;}
	return GM_FALSE;
    }
    inline T& x(int p){return tdata[p].x();}
    inline T& y(int p){return tdata[p].y();}
    inline T& z(int p){return tdata[p].z();}
    Geom_Vertex<T,3>& operator[] (gm_uint r) {
#ifdef __CHECK_BOUNDS
	if (r<0 || r>=3)
	cerr << "Fila fuera de rango " << r << EXIT;
#endif
	return (tdata[r]);
    }
};

template <class T> class Geom_2D_Triangle {

private:
    Geom_Edges<T, 3, 2> tdata;

public:
    Geom_2D_Triangle(){};

    ~Geom_2D_Triangle(){};
	
    void data(T x0, T y0, T x1, T y1, T x2, T y2){
	x(0)=x0; x(1)=x1; x(2)=x2; y(0)=y0; y(1)=y1; y(2)=y2;
    }
	
    inline void data(Geom_Vertex<T,2> p0, Geom_Vertex<T,2> p1, Geom_Vertex<T,2> p2){
	tdata[0] = p0; tdata[1] = p1; tdata[2] = p2;
    }

    inline bool collinear(void){
	if((fabs(y(0)-y(1)) <= GM_EPSILON) && (fabs(y(1)-y(2)) <= GM_EPSILON)){
	    return GM_TRUE;}
	return GM_FALSE;
    }

    inline T& x(int p){return tdata[p].x();}
    inline T& y(int p){return tdata[p].y();}
    
    Geom_Vertex<T,2>& operator[] (gm_uint r) {
#ifdef __CHECK_BOUNDS
	if (r<0 || r>=3)
	cerr << "Fila fuera de rango " << r << EXIT;
#endif
	return (tdata[r]);
    }
    
    Geom_2D_Triangle<T>& operator = (const Geom_2D_Triangle<T>& p1){
        //if (this!=&p1) tdata=p1.data;
	return *this;
    }

};

template <class T> class Geom_3D_Triangle {

private:
    Geom_Edges<T, 3, 3> tdata;

public:
    Geom_3D_Triangle(){};
    ~Geom_3D_Triangle(){};
	
    void data(T x0, T y0, T x1, T y1, T x2, T y2){
    	x(0)=x0; x(1)=x1; x(2)=x2; y(0)=y0; y(1)=y1; y(2)=y2;}
	
    inline void data(Geom_Vertex<T,3> p0, Geom_Vertex<T,3> p1, Geom_Vertex<T,3> p2){
	tdata[0] = p0; tdata[1] = p1; tdata[2] = p2;}

    inline bool collinear(void){
	if((fabs(y(0)-y(1)) <= GM_EPSILON) && (fabs(y(1)-y(2)) <= GM_EPSILON)){
	    return GM_TRUE;}
	return GM_FALSE;}

    inline T& x(int p){return tdata[p].x();}
    inline T& y(int p){return tdata[p].y();}
    inline T& z(int p){return tdata[p].z();}
    
    Geom_3D_Triangle<T>& operator = (const Geom_3D_Triangle<T>& p1){
        if (this!=&p1) tdata=p1.tdata;
	return *this;
    }
    
    Geom_Vertex<T,3>& operator[] (gm_uint r) {
#ifdef __CHECK_BOUNDS
	if (r<0 || r>=3)
	cerr << "Fila fuera de rango " << r << EXIT;
#endif
	return (tdata[r]);
    }
};

//////////////////////////////////////////////////////////////////////

template <class T, gm_uint vsize, gm_uint vdim> class Geom_Group{

private:
	gm_uint _csize;
	//Geom_Vertex<T,vdim> tdata[3];
	vector< Geom_Edges<T, vsize, vdim> > mdata;
	//Geom_Vertex<T,vdim> cc;

public:
    Geom_Group(){
	_csize = 0;
	mdata.resize(_csize);
    }

    Geom_Group<T, vsize, vdim>(gm_uint cs){
	_csize = cs;
	mdata.resize(_csize);
    }
	
    ~Geom_Group(){};

    void add(Geom_Edges<T, vsize, vdim>& _v){
	mdata.push_back(_v);
	_csize++;
    }
    gm_uint size(void){
	return mdata.size();}
    void resize(gm_uint cs){
    	_csize = cs;
	mdata.resize(cs);
    }
    void clear(void){
	mdata.clear();
	_csize = 0;
    }

    Geom_Edges<T, vsize, vdim>& operator[] (gm_uint r) {
#ifdef __CHECK_BOUNDS
	if (r<0 || r>=_csize)
	    cerr << "Fuera de rango " << r << EXIT;
#endif
	return (mdata[r]);
    }
};

/////////////////////////////////////////////////////////////////////////////////

template <class T> class Geom_3D_Mesh {
	gm_uint msize;
	vector< Geom_3D_Triangle<T> > mdata;
	
public:
	Geom_3D_Mesh<T>(){
	};
	
	~Geom_3D_Mesh(){};
	
	void resize(gm_int _ms){
		mdata.resize(msize);
	}

	void push_back(Geom_3D_Triangle<T>& c){
		mdata.push_back(c);
	}

	void clear(){
		mdata.clear();
	}

	gm_uint size(void){
		return mdata.size();
	}
	
	Geom_3D_Triangle<T>& operator[] (gm_uint r) {
#ifdef __CHECK_BOUNDS
		if (r<0 || r>=mdata.size())
			cerr << "Fila fuera de rango " << r << EXIT;
#endif
		return mdata[r];
	}
};

#endif

////
