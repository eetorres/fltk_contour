//========================================================================
// Mini Scientific Matrix and Vector Template Library - MSMVTL
// MSMVTL is a simple  matrix and vector template library
//
// FILE: tmatrix.h
//       the matrix template head file
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


#ifndef _TMATRIX_H_
#define _TMATRIX_H_ 1

//#define _DEBUG  1

#include <msmvtl/tvector.h>

const unsigned int min_rows=1;
const unsigned int min_cols=1;

template <class T> class TMatrix {

private:
    //
  unsigned int u_rr, u_cc;
  int i_rows, i_cols, i_rows_capacity, i_precision, i_width;
  char _buff[1024], c_ch;
  std::valarray< TVector<T> > m_data;

  void set(unsigned int r, unsigned int c) {
    i_cols = c;
    i_rows = r;
    m_data.resize(i_rows);
    for(unsigned int i=0; i<rows(); i++)
      m_data[i].resize(i_cols);
  }

  void reset(unsigned int r, unsigned int c) {
    unsigned int i, j;
    unsigned int d, w;
    std::valarray< TVector<T> > v_t;
    v_t.resize(rows());
    for(i=0; i<rows(); i++){
      v_t[i].resize(cols());
      v_t[i] = m_data[i];
    }
    d = mini(r,rows());
    w = mini(c,cols());
    set(r,c);
    for(i=0; i<d; i++)
      for(j=0; j<w; j++)
         m_data[i][j] = v_t[i][j];
  }

  bool isnrow(unsigned int idx){
    bool inm;
    inm = (idx<0 || idx>=rows());
    if(inm){
#ifdef __GET_MESSAGE
      std::cerr << "Matrix index out of columns range "<<idx<<std::endl<<" EXIT "<<EXIT<<std::endl;
#endif
      exit(EXIT);
    }
    return inm;
  }

  bool isncol(unsigned int idx){
    bool inm;
    inm = (idx<0 || idx>=cols());
    if(inm){
#ifdef __GET_MESSAGE
      std::cerr << "Matrix index out of rows range - "<<idx<<" - EXIT "<<EXIT<<std::endl;
#endif
      exit(EXIT);
    }
    return inm;
  }

  int read_file_rows(const char* f){
    int i=0;
    std::ifstream infile(f);
    if (infile.bad()){
#ifdef _DEBUG
      std::cerr << "The file can't be opened "<<f<<EXIT<<std::endl;
#endif
      return false;
    }
    if (!infile.is_open()){
#ifdef _DEBUG
      std::cerr << "The file can't be opened" <<f<< EXIT<<std::endl;
#endif
      return false;
    }else{
      while(!infile.eof()){
        c_ch = infile.peek();
        while((c_ch == TIL) || (c_ch == NUM)){
          infile.ignore(2048,NLN);
          c_ch = infile.peek();
        }
        while(c_ch == SPC){
          infile.get();
          c_ch = infile.peek();
          if(c_ch == NLN){
            infile.get();
            c_ch = infile.peek();
          }
        }
        if(((c_ch >= '0') && (c_ch <= '9')) || ((c_ch == '-') || (c_ch == '+') || (c_ch == '.') ))
          i++;
        else
          break;
        infile.ignore(2048,NLN);
      }
    }
    infile.close();
    return i;
  }

  int read_file_cols(const char* f){
    int i=0;
    double d;
    std::ifstream infile(f);
    if (!infile.is_open()){
#ifdef _DEBUG
      std::cerr << "The file "<<f<<" can't be opened"<< EXIT<<std::endl;
#endif
      return false;
    }else{
      c_ch = infile.peek();
      while((c_ch == TIL) || (c_ch == NUM)){
        infile.ignore(1024,NLN);
        c_ch = infile.peek();
      }
      while(c_ch == SPC){
        infile.get();
        c_ch = infile.peek();
        if(c_ch == NLN){
          infile.get();
          c_ch = infile.peek();
        }
      }
      while(c_ch != NLN){
        i++;
        infile >> d;
        c_ch = infile.peek();
        while(c_ch == SPC){
          infile.get();
          c_ch = infile.peek();
          if(c_ch == NLN)
            break;
        }
      }
    }
    infile.close();
    return i;
  }

  void read_(const char* f, unsigned int r, unsigned int c){
    T d;
    resize(r,c);
    std::ifstream infile(f);
    for(unsigned int i=0; i<r; i++){
      c_ch = infile.peek();
      while((c_ch == TIL) || (c_ch == NUM)){
        infile.ignore(1024,NLN);
        c_ch = infile.peek();
      }
      while(c_ch == SPC){
        infile.get();
        c_ch = infile.peek();
        if(c_ch == NLN){
          infile.get();
          c_ch = infile.peek();
        }
      }
      for(unsigned int j=0; j<c; j++){
        infile >> d;
        m_data[i][j]=d;
      }
      infile.ignore(1024,NLN);
    }
    infile.close();
  }

public:
  //
  TMatrix(){
    init();
    set(0,0);
  }

  TMatrix(unsigned int r, unsigned int c){
    init();
    set(r,c);
  }

  TMatrix(const char* f){
    init();
    read_file(f);
  }

  TMatrix(const TMatrix<T>& m){
    init();
    set(m.i_rows, m.i_cols);
    for(unsigned int i=0; i<rows(); i++)
      m_data[i] = m.m_data[i];
  }

  void init(void){
    i_precision=16;
    i_width=12;
  }

  bool read_file(const char* f){
    return read_file(f,0,0);
  }

  bool read_file(const char* f, unsigned int r, unsigned int c){
    std::ifstream infile(f);
    if(!infile.is_open()){
#ifdef _DEBUG
      std::cerr << "The file "<<f<<" can not be opened "<< EXIT<<std::endl;
#endif
      return false;
    }
    infile.close();
#ifdef _DEBUG
    u_rr = read_file_rows(f);
    printf("rows in file %i\n",u_rr);
#endif
#ifdef _DEBUG
    u_cc = read_file_cols(f);
    printf("cols in file = %i\n",u_cc);
#endif
    if(r <= u_rr && c <= u_cc){
      if((r>0) && (c>0)){
        read_(f,r,c);
      }else if((r>0) && !(c>0)){
        u_cc = read_file_cols(f);
        read_(f,r,u_cc);
      }else if(!(r>0) && (c>0)){
        u_rr = read_file_rows(f);
        read_(f,u_rr,c);
      }else{
        u_rr = read_file_rows(f);
        u_cc = read_file_cols(f);
        read_(f,u_rr,u_cc);
      }
    }else{
      return false;
    }
    return true;
  }

  void write_file(char* f){
    write_file(f, i_precision);
  }


  void write_file(char* f, int _p){
    std::ofstream outfile (f);
    precision(_p);
    outfile.precision(i_precision);
    //outfile.setf(std::ios::right, std::ios::adjustfield);
    outfile.setf(std::ios::fixed, std::ios::floatfield);
    if (outfile.bad()){
      std::cerr<<" It's not enable to open the file "<<f<<std::endl;
      std::cerr<<" EXIT "<<EXIT<<std::endl;
      exit(EXIT);
    }
    for(unsigned int i=0; i<rows(); i++){
      for(unsigned int j=0; j<cols(); j++){
         outfile<< m_data[i][j];
         if(j<cols()-1)
           outfile<<"\t";
      }
      outfile<< "\n";
    }
    outfile.close();
#ifdef _DEBUG
    std::cerr<<" The file was successful saved"<<std::endl;
#endif
  }

  // 15/09/2011: obsolete function use the next one
  void zero(void){
    clear();
  }

  void clear(void){
    set(0,0);
  }

  void copy(const real _m[3][3], unsigned int r, unsigned int c){
    resize(r,c);
    for(unsigned int i=0; i<r; i++)
      for(unsigned int j=0; j<c; j++)
        m_data[i][j]=_m[i][j];
  }

  void mrand(void){
    for(unsigned int i=0; i<rows(); i++){
      srand((unsigned)time(0)+i);
      m_data[i].vrand();
    }
  }

  // Wed Nov 23 17:11:13 MST 2011: obsolete function, use the next one
  bool if_square(void){
    return is_square();
  }

  bool is_square(void){
    return (rows() == cols());
  }

  // Wed Nov 23 17:11:13 MST 2011: obsolete function, use the next one
  bool if_empty(void){
    return is_empty();
  }

  bool is_empty(void){
    return bool(rows()*cols());
  }

  TMatrix<T> inverse(void){
    T q;
    TMatrix<T> m_t;
    m_t = (*this);
    for(unsigned int i= 0; i < rows(); i++){
      //q = ABS(m_t[i][i]);
      //if(q <= _EPSILON) return 0;
      q = (1.0/m_t[i][i]);
      m_t[i][i] = 1.0;
      for(unsigned int k = 0; k<rows(); k++){
        m_t[i][k] *= q;
      }
      for(unsigned int j = 0; j<rows(); j++){
        if (j!=i){
          q = m_t[j][i];
          m_t[j][i] = 0.0;
          for(unsigned int k = 0; k < rows(); k++){
          m_t[j][k] -= (q * m_t[i][k]);
          }
        }
      }
    }
  return m_t;
  }

  // BETA.... last test 08/03/2011
  // working good.
  bool transpose(void){
    real r_t;
    TMatrix<T> m_t;
    if( if_square() ){
      for(unsigned int i=0; i<rows(); i++){
        for(unsigned int j=0; j<i; j++){
          if(i!=j){
            r_t = m_data[i][j];
            m_data[i][j] = m_data[j][i];
            m_data[j][i] = r_t;
          }
        }
      }
    }else{
      m_t = (*this);
      resize(m_t.cols(),m_t.rows());
      for(unsigned int i=0; i<m_t.rows(); i++){
        for(unsigned int j=0; j<m_t.cols(); j++){
          //if(i!=j)
            //_tmp = m_data[i][j];
          m_data[j][i] = m_t[i][j];
          //m_data[j][i] = _tmp;
        }
      }
    }
    return (true);
  }

  void mconst(real r){
    for(unsigned int i=0; i<rows(); i++)
      m_data[i].vconst(r);
  }

  void main_diagonal(T r){
    if(if_square())
      for(unsigned int i=0; i<rows(); i++)
        m_data[i][i]=r;
  }

  // alpha: added 17/03/2011
  unsigned int rows(void) const {
    return rows();
  }

  unsigned int rows(void){
    return i_rows;
  }

  // alpha: added 17/03/2011
  unsigned int cols(void) const {
    return cols();
  }

  unsigned int cols(void){
    return i_cols;
  }

  TVector<T> abs_cols(void){
    TVector<T> v_t(rows());
    for(unsigned int i=0; i<rows(); i++)
      v_t[i] = m_data[i].abs();
    return v_t;
  }

  TVector<T> norm_cols(void){
    TVector<T> v_t(rows());
    for(unsigned int i=0; i<rows(); i++)
      v_t[i] = m_data[i].norm();
    return v_t;
  }

  TVector<T> sum_rows(void){
    TVector<T> v_t1(cols());
    TVector<T> v_t2(rows());
    for(unsigned int i=0; i<cols(); i++){
      v_t2 = get_col(i);
      v_t1[i] = v_t2.sum();
    }
    return v_t1;
  }

  TVector<T> sum_cols(void){
    TVector<T> v_t(rows());
    for(unsigned int i=0; i<rows(); i++)
      v_t[i] = m_data[i].sum();
    return v_t;
  }

  void sort_col(unsigned int c){
    TVector<T> v_t;
    for(unsigned int i=0; i<rows(); i++){
      for(unsigned int j=0; j<rows(); j++){
        if(m_data[i][c] > m_data[j][c]){
          v_t = m_data[i];
          m_data[i] = m_data[j];
          m_data[j] = v_t;
        }
      }
    }
  }

  void add_row(TVector<T>& v){
    if( cols() == v.size() ){
      reset(rows()+1, cols());
      m_data[rows()-1]=v;
    }else{
#ifdef __GET_MESSAGE
      printf("number of columns must be equals\n");
#endif
    }
  }

  void add_rows(TMatrix<T>& m){
    unsigned int u = rows();
    if( cols() == m.cols() ){
      reset(rows()+m.rows(), cols());
      for(unsigned int i=0; i<m.rows(); i++)
        m_data[u+i]=m[i];
    }else{
#ifdef __GET_MESSAGE
      printf("number of columns must be equals\n");
#endif
    }
  }

  void add_col(TVector<T>& v){
    if( rows() == v.size() ){
    for(unsigned int i=0; i<rows(); i++)
      m_data[i].push_back(v[i]);
      i_cols++;
    }else{
#ifdef __GET_MESSAGE
      printf("number of rows must be equals\n");
#endif
    }
  }

  void add_cols(TMatrix<T>& m){
    TVector<T> v_t;
    if( rows() == m.rows() ){
      for(unsigned int i=0; i<m.cols(); i++){
        v_t = m.get_col(i);
        add_col(v_t);
      }
    }else{
#ifdef __GET_MESSAGE
      printf("number of rows must be equals\n");
#endif
    }
  }

  TVector<T> get_col(int i) {
    TVector<T> v_t(rows());
    if(!isncol(i)){
      for(unsigned int r=0; r<rows(); r++)
        v_t[r] = m_data[r][i];
    }
    return v_t;
  }

  TVector<T> get_col(int i) const{
    return get_col(i);
  }

  TVector<T> get_diag(void){
    TVector<real> v_t(rows());
    if(if_square())
      for(int r=0; r<rows(); r++)
        v_t[r] = m_data[r][r];
    return v_t;
  }

  void del_row(unsigned int u){
    for(unsigned int i=u; i<rows()-1; i++)
      m_data[i] =m_data[i+1];
    i_rows--;
    reset(i_rows,cols());
  }

  void del_row(void){
    del_row(rows()-1);
  }

  void del_col(unsigned int c){
    for(unsigned int r=0; r<rows(); r++){
      m_data[r][c]=m_data[r][cols()-1];
      m_data[r].remove();
    }
    i_cols--;
  }

  void precision(int i){
    i_precision=i;
  }

  void width(int i){
    i_width=i;
  }

  void del_col(void){
    del_col(cols()-1);
  }

  void resize(int new_rows, int new_cols){
    set(new_rows, new_cols);
  }

  void redim(int new_rows, int new_cols){
    reset(new_rows, new_cols);
  }

  // alpha: added 20/03/2011
  T col_max(unsigned int i){
    TVector<T> _v1=get_col(i);
    T _max = _v1.get_max();
    return _max;
  }

  // alpha: added 20/03/2011
  T col_min(unsigned int i){
    TVector<T> _v1=get_col(i);
    T _min = _v1.get_min();
    return _min;
  }

  T get_max(void){
    double _max = rows()>0?m_data[0].get_max():0;
    for(unsigned int i=1; i<rows(); i++)
      _max = maxi(_max,m_data[i].get_max());
    return (_max);
  }

  T get_min(void){
    double _mm, _min;
    _min = rows()>0?m_data[0].get_min():0;
    for(unsigned int i=1; i<rows(); i++){
      _mm = m_data[i].get_min();
      //printf("%i)_mm = %f\n",i,_mm);
      _min = mini(_min,_mm);
    }
    return (_min);
  }

  std::string dimension(void){
    std::string o;
    if(if_empty()){
      sprintf(_buff,"{ %u x %u }",rows(),cols());
      o = _buff;
    }else{ o="empty";}
    return o;
  }

  TVector<T>& operator[] (const int i){
#ifdef __CHECK_BOUNDS
    isnrow(i);
#endif
    return m_data[i];
  }

  // beta: added 08/03/2011
  const TVector<T>& operator[] (const int i) const{
    return m_data[i];
  }

  TMatrix<T>& operator= (const TMatrix<T>& m){
    if(this != &m){
      if(this->i_rows != m.i_rows || this->i_cols != m.i_cols)
        set(m.i_rows, m.i_cols);
      m_data = m.m_data;
    }
    return *this;
  }

  // modified 08/03/2011
  friend TMatrix<T> operator* (const TMatrix<T>& m1, const TMatrix<T>& m2){
    TMatrix<T> v_t(m1.i_rows, m2.i_cols), m_t;
    //TVector<T> _v1;
    m_t = m2;
    m_t.transpose();
    for (int i=0; i<m1.i_rows; i++)
      for (int j=0; j<m2.i_cols; j++){
        //_v1 = _m2.get_col(j);
        v_t[i][j] = (m1[i]*m_t[j]);
      }
    return v_t;
  }

  // 16/03/2011: inconsistent operation, vectors are single column matrices
  friend TVector<T> operator* (const TMatrix<T>& m, const TVector<T>& v){
    TVector<T> v_t(m.i_rows);
    //if(_m1.cols == v_t.size()){
    for (int i=0; i<m.i_rows; i++)
      v_t[i] = (m[i] * v);
    return v_t;
    //}
    //return 0;
  }

  // added 16/03/2011
  friend TVector<T> operator* (const TVector<T>& v, const TMatrix<T>& m){
    TVector<T> v_t(m.i_rows);
    TMatrix<T> m_t;
    m_t = m;
    m_t.transpose();
    for (int i=0; i<m.i_rows; i++)
      v_t[i] = (v*m_t[i]);
    return v_t;
  }

  friend TMatrix<T> operator* (const T d, const TMatrix<T>& m){
    TMatrix<T> m_t;
    m_t.set(m.i_rows, m.i_cols);
    for(int i=0; i<m.i_rows; i++)
      m_t.m_data[i] = (d*m.m_data[i]);
    return m_t;
   }

  friend TMatrix<T> operator/ (const TMatrix<T>& m, const T d){
    TMatrix<T> v_t;
    v_t.set(m.i_rows,m.i_cols);
    T inv_d = 1.0/d;
    for(int i=0; i<m.i_rows; i++)
      v_t.m_data[i] = m.m_data[i]*inv_d;
    return v_t;
  }

  friend TMatrix<T> operator+ (const TMatrix<T>& m1, const TMatrix<T>& m2){
    TMatrix<T> v_t;
    v_t.set(m1.i_rows, m1.i_cols);
    for(int i=0; i<m1.i_rows; i++)
      v_t.m_data[i] = m1.m_data[i] + m2.m_data[i];
    return v_t;
  }

  friend TMatrix<T> operator+ (const TMatrix<T>& m, const TVector<T>& v){
    TMatrix<T> v_t;
    v_t.set(m.i_rows, m.i_cols);
    for(int i=0; i<m.i_rows; i++)
      v_t.m_data[i] = m.m_data[i] + v;
    return v_t;
  }

  friend TMatrix<T> operator- (const TMatrix<T>& m, T r){
    TMatrix<T> v_t;
    v_t.set(m.i_rows, m.i_cols);
    for(int i=0; i<m.i_rows; i++)
      v_t.m_data[i] = m.m_data[i]-r;
    return v_t;
  }

  friend TMatrix<T> operator- (const TMatrix<T>& m, const TVector<T>& v){
    TMatrix<T> v_t;
    v_t.set(m.i_rows, m.i_cols);
    for(int i=0; i<m.i_rows; i++)
      v_t.m_data[i] = m.m_data[i] - v;
    return v_t;
  }

  friend TMatrix<T> operator- (const TMatrix<T>& m1, const TMatrix<T>& m2){
    TMatrix<T> v_t;
    v_t.set(m1.i_rows, m1.i_cols);
    for(int i=0; i<m1.i_rows; i++)
      v_t.m_data[i] = m1.m_data[i] - m2.m_data[i];
    return v_t;
  }

  friend std::ostream& operator << (std::ostream& o, TMatrix<T>& m){
    o <<" {"<<m.i_rows<<"x"<<m.i_cols<<"} ="<<std::endl;
    o.precision(m.i_precision);
    for (int i=0; i < m.i_rows; i++){
      o << " ";
      for (int j=0; j < m.i_cols; j++){
        //if(_m.m_data[i][j]>=0) _o << " ";
        o.width(m.i_width);
        o <<std::fixed<<std::right<< m.m_data[i][j]<<" ";
      }
      o << std::endl;
    }
    return o;
  }

};

#endif
