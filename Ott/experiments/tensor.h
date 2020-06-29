#ifndef TENSOR_H
#define TENSOR_H

#undef INDEX_CHECK
#include <iostream>
#include <cstdlib>

class vector
{
	int n;
	double *v;
public:
	vector();
	vector(int n);
	vector(int n, double initial);
	vector(const vector &a);
	~vector();
	int size() const;
	void resize(int n);
	double *get_v();
	vector &operator = (const vector &rhs);
	vector &operator = (const double &rhs);
	double &operator[] (const int &index);
	double &operator() (const int &index);
	const double &operator[] (const int &index) const;
	const double &operator() (const int &index) const;
private:
	void index_error() const;
};

std::ostream& operator <<(std::ostream& os, const vector &a);
std::istream& operator >>(std::istream& is, vector &a);
vector operator + (const vector &a, const vector &b);
vector operator - (const vector &a, const vector &b);
vector operator * (double scalar, const vector &a);
vector operator * (const vector &a, double scalar);
vector operator / (const vector &a, double scalar);
vector &operator += (vector &a, const vector &b);
vector &operator -= (vector &a, const vector &b);
vector &operator *= (vector &a, double scalar);
vector &operator /= (vector &a, double scalar);
double operator * (const vector &a, const vector &b);
void add(const vector &a, const vector &b, vector &c);
void sub(const vector &a, const vector &b, vector &c);
void mul(double scalar, const vector &a, vector &b);
void mul(double scalar, vector &a);
void div(const vector &a, double scalar, vector &b);
void div(vector &a, double scalar);
double dot(const vector &a, const vector &b);
void vector_error();


class matrix
{
	int nrows;
	vector *m;
public:
	matrix();
	matrix(int nrows);
	matrix(int nrows, int ncols);
	matrix(int nrows, int nclos, double initial);
	matrix(const matrix &a);
	~matrix();
	int size() const;
	int size_row() const;
	int size_col() const;
	void resize(int nrows);
	void resize(int nrows, int ncols);
	vector *get_m();
	matrix &operator = (const matrix &rhs);
	matrix &operator = (const double &rhs);
	vector &operator[] (const int &index);
	vector &operator() (const int &index);
	const vector &operator[] (const int &index) const;
	const vector &operator() (const int &index) const;
	double &operator() ( const int &index_row, const int &index_col );
	const double &operator() ( const int &index_row, const int &index_col ) const;
private:
	void index_error() const;
};

std::ostream& operator <<(std::ostream& os, const matrix &a);
std::istream& operator >>(std::istream& is, matrix &a);
matrix operator + (const matrix &a, const matrix &b);
matrix operator - (const matrix &a, const matrix &b);
matrix operator * (double scalar, const matrix &a);
matrix operator * (const matrix &a, double scalar);
vector operator * (const matrix &m, const vector &v);
vector operator * (const vector &v, const matrix &m);
matrix operator / (const matrix &a, double scalar);
matrix &operator += (matrix &a, const matrix &b);
matrix &operator -= (matrix &a, const matrix &b);
matrix &operator *= (matrix &a, double scalar);
matrix &operator /= (matrix &a, double scalar);
void add(const matrix &a, const matrix &b, matrix &c);
void sub(const matrix &a, const matrix &b, matrix &c);
void mul(double scalar, const matrix &a, matrix &b);
void mul(double scalar, matrix &a);
void mul(const matrix &m, const vector &va, vector &vb);
void mul(const vector &va, const matrix &m, vector &vb);
void div(const matrix &a, double scalar, matrix &b);
void div(matrix &a, double scalar);
void transpose( const matrix &a, matrix &aTrans );
void transpose( matrix &a );
double trace( const matrix &m );
double norm( const matrix &m );
void matrix_error();

/////////////////////////////////////// inlne function definitions ////////////////////////////////////////////

// vector functions

inline vector::vector() : n(0), v(0) {}

inline vector::vector(int n) : n(n) { v = new double[n]; }

inline vector::vector(int n, double initial) : n(n)
{
	v = new double[n];
	for(int i=0; i<n; ++i) v[i] = initial;
}

inline vector::vector(const vector &a) : n(a.n)
{
	v = new double[n];
	for(int i=0; i<n; ++i) v[i] = a.v[i];
}

inline vector::~vector() { delete[] v; }

inline int vector::size() const { return n; }

inline void vector::resize(int n)
{
	if(vector::n != n) { delete[] v; vector::n = n; v = new double[n]; }
}

inline double *vector::get_v()
{
	return v;
}

inline double &vector::operator[] (const int &index)
{
	#ifdef INDEX_CHECK
	if (index<0 || index>=n) index_error();
	#endif
	return v[index];
}

inline double &vector::operator() (const int &index)
{
	#ifdef INDEX_CHECK
	if (index<0 || index>=n) index_error();
	#endif
	return v[index];
}

inline const double &vector::operator[] (const int &index) const
{
	#ifdef INDEX_CHECK
	if (index<0 || index>=n) index_error();
	#endif
	return v[index];
}

inline const double &vector::operator() (const int &index) const
{
	#ifdef INDEX_CHECK
	if (index<0 || index>=n) index_error();
	#endif
	return v[index];
}

// matrix functions

inline matrix::matrix() : nrows(0), m(0) {};

inline matrix::matrix(int nrows) : nrows(nrows)
{
	m = new vector[nrows];
}

inline matrix::matrix(int nrows, int ncols) : nrows(nrows)
{
	m = new vector[nrows];
	for(int i=0; i<nrows; ++i) m[i].resize(ncols);
}

inline matrix::matrix(int nrows, int ncols, double initial) : nrows(nrows)
{
	m = new vector[nrows];
	for(int i=0; i<nrows; ++i) { m[i].resize(ncols); m[i] = initial; }
}

inline matrix::matrix(const matrix &a) : nrows(a.nrows)
{
	m = new vector[nrows];

	for(int i=0; i<nrows; ++i) m[i] = a.m[i];
}

inline matrix::~matrix() { delete[] m; }

inline int matrix::size() const { return nrows; }

inline int matrix::size_row() const { return nrows; }

inline int matrix::size_col() const { return ( nrows == 0 ? 0 : m[0].size() ); }

inline void matrix::resize(int nrows)
{
	if(matrix::nrows != nrows) {
		delete[] m;
		matrix::nrows = nrows;
		m = new vector[nrows];
	}
}

inline void matrix::resize(int nrows, int ncols)
{
	if(matrix::nrows != nrows) {
		delete[] m;
		matrix::nrows = nrows;
		m = new vector[nrows];
	}
	for(int i=0; i<nrows; ++i) m[i].resize(ncols);
}

inline vector *matrix::get_m()
{
	return m;
}

inline vector &matrix::operator[] (const int &index)
{
	#ifdef INDEX_CHECK
	if (index<0 || index>=nrows) index_error();
	#endif
	return m[index];
}

inline vector &matrix::operator() (const int &index)
{
	#ifdef INDEX_CHECK
	if (index<0 || index>=nrows) index_error();
	#endif
	return m[index];
}

inline const vector &matrix::operator[] (const int &index) const
{
	#ifdef INDEX_CHECK
	if (index<0 || index>=nrows) index_error();
	#endif
	return m[index];
}

inline const vector &matrix::operator() (const int &index) const
{
	#ifdef INDEX_CHECK
	if (index<0 || index>=nrows) index_error();
	#endif
	return m[index];
}

inline double &matrix::operator() ( const int &index_row, const int &index_col )
{
	#ifdef INDEX_CHECK
	if( index_row<0 || index_row>=nrows ) index_error();
	#endif    
	return m[ index_row ][ index_col ];
}

inline const double &matrix::operator() ( const int &index_row, const int &index_col ) const
{
	#ifdef INDEX_CHECK
	if( index_row<0 || index_row>=nrows ) index_error();
	#endif    
	return m[ index_row ][ index_col ];
}

#endif
