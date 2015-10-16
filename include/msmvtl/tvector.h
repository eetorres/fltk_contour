//========================================================================
// Mini Scientific Matrix and Vector Template Library - MSMVTL
// MSMVTL is a simple  matrix and vector template library
//
// FILE: tvector.h
//       the vector template head file
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


#ifndef _TVECTOR_H_
#define _TVECTOR_H_ 1

#include <msmvtl/libdef.h>
#include <valarray>
#include <vector>

const short min_cap=1;

template <class T> class TVector{

protected:
  uint u_rr;
  std::valarray<T> v_data;

private:
  char ch;
  uint i_prs;

  void resize_(uint u){
    unsigned int u_d;
    std::valarray<T> v_t = v_data;
    v_data.resize(u);
    u_d = mini(v_t.size(),u);
    for(unsigned int i=0; i<u_d; i++)
      v_data[i]=v_t[i];
  }

  bool isndata(uint u){
    bool b;
    b = (u<0 || u>=v_data.size());
    if(b){
#ifdef __GET_MESSAGE
      std::cerr << "Index out of range - "<<u<<" - EXIT Vector "<<EXIT<<std::endl;
#endif
      exit(EXIT);
    }
    return b;
  }

  bool isndata(uint u) const{
    return isndata(u);
  }

  int read_file_rows(const char* f){
    int i=0;
    std::ifstream infile(f);
    if (infile.bad()){
      std::cerr<<" The file can't be opened "<<f<< EXIT;
      return false;
    }
    if (!infile.is_open()){
      std::cerr<<" The file can't be opened "<<f<< EXIT;
      return false;
    }else{
      while(!infile.eof()){
        ch = infile.peek();
        while((ch == TIL) || (ch == NUM)){
          infile.ignore(1024,NLN);
          ch = infile.peek();
        }
        while(ch == SPC){
          infile.get();
          ch = infile.peek();
          if(ch == NLN){
            infile.get();
            ch = infile.peek();
          }
        }
        if(((ch >= '0') && (ch <= '9')) || ((ch == '-') || (ch == '+')))
          i++;
        else
          break;
        infile.ignore(1024,NLN);
      }
    }
    infile.close();
    return i;
  }

  void read_(const char* f, uint r){
    T d;
    resize_(r);
    std::ifstream infile(f);
    for(unsigned int i=0; i<r; i++){
      ch = infile.peek();
      while((ch == TIL) || (ch == NUM)){ 
        infile.ignore(1024,NLN);
        ch = infile.peek();
      }
      while(ch == SPC){
        infile.get();
        ch = infile.peek();
        if(ch == NLN){
          infile.get();
          ch = infile.peek();
        }
      }
      infile >> d;
      v_data[i]=d;
      infile.ignore(1024,NLN);
    }
    infile.close();
  }

public:
  //
  TVector(){
    resize(0);
    init();
  }

  TVector(uint u){
    resize(u);
    init();
  }

  TVector(const TVector<T>& v){
    resize(v.v_data.size());
    for(uint i=0;i<v.v_data.size();i++)
      v_data[i] = v.v_data[i];
    init();
  }

  TVector(char* f){
    read_file(f,0);
    init();
  }

  TVector(char* f, int i){
    read_file(f,i);
    init();
  }

  void init(void){
    i_prs=6;
  }

  unsigned int size(void){
    return v_data.size();
  }

  void zero(void){
    constant(0);
  }

  void clear(void){
    v_data.resize(0);
  }

  void resize(uint u){
    resize_(u);
  }

  void push_back(T d){
    resize_(size()+1);
    v_data[size()-1]=d;
  }

  void push_back(TVector<T> v){
    for(uint i=0; i<v.size(); i++){
      resize_(size()+1);
      v_data[size()-1]=v.v_data[i];
    }
  }

  void vrand(void){
    for(uint i=0; i<size(); i++){
      v_data[i] = ((real)rand()/(RAND_MAX/2.0))-1.0;
    }
  }

  void svrand(uint u){
    srand(u);
    vrand();
  }

  void constant(real d){
    for (uint i=0; i<size(); i++)
      v_data[i] = d;
  }

  void sort(void){
    T d;
    for(int i=0;i<size()-1;i++){
      for(int j=i+1;j<size();j++){
        if(v_data[i] > v_data[j]){
          d = v_data[i];
          v_data[i] = v_data[j];
          v_data[j] = d;
        }
      }
    }
  }

  T sum(void){
    T d=0;
    for (unsigned int _i=0; _i<size(); _i++)
      d+=v_data[_i];
    return d;
  }

  T abs(void){
    T d=0;
    for (unsigned int _i=0; _i<size(); _i++)
      d+=ABS(v_data[_i]);
    return d;
  }

  T norm(void){
    T d=0;
    for (unsigned int i=0; i<size(); i++)
      d+=(v_data[i]*v_data[i]);
    return d;
  }

  T magnitude(void){
    T _magnitude = sqrt(norm());
    return _magnitude;
  }

  void remove(int i){
    if(!isndata(i)){
      v_data[i] = v_data[size()-1];
      resize_(size()-1);
    }
  }

  void remove(void){
    remove(size()-1);
  }

  T get_min(void){
    T _min=size()>0?v_data[0]:0;
    for (unsigned int i=1; i<size(); i++)
      _min=mini(v_data[i],_min);
    return _min;
  }

  T get_max(void){
    T _max=size()>0?v_data[0]:0;
    for (unsigned int i=1; i<size(); i++)
      _max=maxi(_max,v_data[i]);
    return _max;
  }

  T first(void){
    return v_data[0];
  }

  T last(void){
    return v_data[size()-1];
  }

  void precision(int i){
    i_prs=i;
  }

  bool read_file(char* f){
    return read_file(f, 0);
  }

  bool read_file(const char* f, uint u){
    resize_(0);
    std::ifstream infile(f);
    if (infile.bad()){
      std::cerr << " The file can't be opened "<<f<< EXIT;
      return false;
    }
    if (!infile.is_open()){
      std::cerr << " The file can't be opened "<<f<< EXIT;
      return false;
    }
    if(u>0){
      read_(f,u);
    }else{
      u_rr = read_file_rows(f);
      read_(f,u_rr);
    }
    infile.close();
    return true;
  }

  void write_file(char* f){
    write_file(f, 9);
  }

  void write_file(char* f, int u){
    //uint i;
    std::ofstream outfile(f);
    outfile.precision(u);
    //outfile.setf(iosstd::right, iosstd::adjustfield);
    outfile.setf(std::ios::fixed, std::ios::floatfield);
    if (outfile.bad()){
      std::cerr<<"The file can't be created "<<f<<std::endl<<" EXIT "<<EXIT<<std::endl;
      exit(EXIT);
    }
    for (uint i=0; i<size(); i++){
      //outfile.width(_sci);
      outfile << v_data[i] << "\n";
    }
    outfile.close();
#ifdef _DEBUG
    printf("The file was successful saved\n");
#endif
  }

  T& operator[] (uint u){
#ifdef __CHECK_BOUNDS
    isndata(u);
#endif
    return v_data[u];
  }

  const T& operator[] (uint u) const{
    return v_data[u];
  }

  TVector<T>& operator= (const TVector<T>& v){
    if (this != &v){
      v_data.resize(v.v_data.size());
      v_data = v.v_data;
    }
    return *this;
  }

  TVector<T>& operator+= (const TVector<T>& v){
    if (this != &v)
      for (unsigned int i=0; i<size(); i++)
        v_data[i] = (v_data[i] + v.v_data[i]);
    return *this;
  }

  TVector<T>& operator-= (const TVector<T>& v){
    if (this != &v)
     for (unsigned int i=0; i<size(); i++)
       v_data[i] = (v_data[i] - v.v_data[i]);
    return *this;
  }

  friend TVector<T> operator+ (const TVector<T>& v1, const TVector<T>& v2){
    TVector<T> v_t;
    v_t.resize(v1.v_data.size());
    for (unsigned int i=0; i<v1.v_data.size(); i++)
      v_t.v_data[i] = (v1.v_data[i] + v2.v_data[i]);
    return v_t;
  }

  friend TVector<T> operator+ (const TVector<T>& v, const T d){
    TVector<T> v_t;
    v_t.v_data.resize(v.v_data.size());
    for (unsigned int i=0; i<v.v_data.size(); i++)
      v_t[i] = (v[i] + d);
    return v_t;
  }

  friend TVector<T> operator+ (T d, const TVector<T>& v){
    TVector<T> v_t;
    v_t = (v + d);
    return v_t;
  }

  friend TVector<T> operator- (const TVector<T>& v1, const TVector<T>& v2){
    TVector<T> v_t;
    v_t.v_data.resize(v1.v_data.size());
    for (unsigned int i=0; i<v1.v_data.size(); i++)
      v_t.v_data[i] = (v1.v_data[i]-v2.v_data[i]);
    return v_t;
  }

  friend TVector<T> operator- (const TVector<T>& v, const T d){
    TVector<T> v_t;
    v_t.v_data.resize(v.v_data.size());
    for (unsigned int i=0; i<v.v_data.size(); i++)
      v_t.v_data[i] = (v.v_data[i]-d);
    return v_t;
  }

  friend TVector<T> operator- (const T d, const TVector<T>& v){
    TVector<T> v_t;
    v_t = (v - d);
    return v_t;
  }

  friend TVector<T> operator- (const TVector<T>& v){
    TVector<T> v_t;
    v_t.v_data.resize(v.v_data.size());
    for (unsigned int i=0; i<v.v_data.size(); i++)
      v_t.v_data[i] = -v.v_data[i];
    return v_t;
  }

  friend T operator* (const TVector<T>& v1, const TVector<T>& v2){
    double s=0;
    for(unsigned int i=0; i<v1.v_data.size(); i++)
      s += (v1.v_data[i]*v2.v_data[i]);
    return s;
  }

  friend TVector<T> operator* (const T d, const TVector<T>& v){
    TVector<T> v_t;
    v_t.v_data.resize(v.v_data.size());
    for (uint i=0; i<v.v_data.size(); i++)
      v_t.v_data[i] = (d*v.v_data[i]);
    return v_t;
  }

  friend TVector<T> operator* (const TVector<T>& v, const T d){
    TVector<T> v_t;
    v_t.v_data.resize(v.v_data.size());
    for (uint i=0; i<v.v_data.size(); i++)
      v_t.v_data[i] = (d*v.v_data[i]);
    return v_t;
  }

  friend TVector<T> operator/ (const TVector<T>& v1, const TVector<T>& v2){
    TVector<T> v_t;
    v_t.v_data.resize(v1.v_data.size());
    for (unsigned int i=0; i<v1.v_data.size(); i++)
      v_t.v_data[i] = (v1.v_data[i]/v2.v_data[i]);
    return v_t;
  }

  friend TVector<T> operator/ (const TVector<T>& v, const T d){
    TVector<T> v_t;
    v_t.v_data.resize(v.v_data.size());
    for (unsigned int i=0; i<v.v_data.size(); i++)
      v_t.v_data[i] = (v.v_data[i]/d);
    return v_t;
  }

  friend std::ostream& operator<< (std::ostream& o, TVector<T>& v){
    o <<" {"<<v.size()<<"} ="<<std::endl;
    o.precision(v.i_prs);
    if(v.size()){
      o<<"  "; // Offset
      for (unsigned int i=0; i < v.size(); i++){
        o <<std::fixed<<std::right<< v.v_data[i];
        if(i<v.size()-1) o<<"  ";
      }
      o <<std::endl;
    }else{
      o<<"empty"<<std::endl;
    }
    return o;
  }

};

#endif
