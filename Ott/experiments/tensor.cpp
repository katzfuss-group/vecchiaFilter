#include <cmath>
#include "tensor.h"

// vector functions

vector &vector::operator = (const vector &rhs)
{
	resize(rhs.n);

	for(int i=0; i<n; ++i) v[i] = rhs.v[i];

	return *this;
}

vector &vector::operator = (const double &rhs)
{
	for(int i=0; i<n; ++i) v[i] = rhs;

	return *this;
}

void vector::index_error() const
{
	std::cerr << "vector index out of range" << '\n'; exit(1);
}

std::ostream& operator <<(std::ostream& os, const vector &a)
{
	int n = a.size();
	for( int i=0; i<n; ++i ) os << a[i] << '\t';

	return os;
}

std::istream& operator >>(std::istream& is, vector &a)
{
	int n = a.size();
	for( int i=0; i<n; ++i ) is >> a[i];
	
	return is;
}

vector operator + (const vector &a, const vector &b)
{
	int na = a.size();
	#ifdef INDEX_CHECK
	int nb = b.size();
	if (na != nb) vector_error();
	#endif
	vector c(na);
	for(int i=0; i<na; ++i) c[i] = a[i] + b[i];

	return c;
}

vector operator - (const vector &a, const vector &b)
{
	int na = a.size();
	#ifdef INDEX_CHECK
	int nb = b.size();
	if (na != nb) vector_error();
	#endif
	vector c(na);
	for(int i=0; i<na; ++i) c[i] = a[i] - b[i];

	return c;
}

vector operator * (double scalar, const vector &a)
{
	int n = a.size();
	vector b(n);
	for(int i=0; i<n; ++i) b[i] = scalar * a[i];

	return b;
}

vector operator * (const vector &a, double scalar)
{
	int n = a.size();
	vector b(n);
	for(int i=0; i<n; ++i) b[i] = scalar * a[i];

	return b;
}

vector operator / (const vector &a, double scalar)
{
	int n = a.size();
	vector b(n);
	for(int i=0; i<n; ++i) b[i] = a[i] / scalar;

	return b;
}

vector &operator += (vector &a, const vector &b)
{
	int na = a.size();
	#ifdef INDEX_CHECK
	int nb = b.size();
	if (na != nb) vector_error();
	#endif
	for(int i=0; i<na; ++i) a[i] += b[i];

	return a;
}

vector &operator -= (vector &a, const vector &b)
{
	int na = a.size();
	#ifdef INDEX_CHECK
	int nb = b.size();
	if (na != nb) vector_error();
	#endif
	for(int i=0; i<na; ++i) a[i] -= b[i];

	return a;
}

vector &operator *= (vector &a, double scalar)
{
	int n = a.size();
	for(int i=0; i<n; ++i) a[i] *= scalar;    

	return a;
}

vector &operator /= (vector &a, double scalar)
{
	int n = a.size();
	for(int i=0; i<n; ++i) a[i] /= scalar;    

	return a;
}

double operator * (const vector &a, const vector &b)
{
	int na = a.size();
	#ifdef INDEX_CHECK
	int nb = b.size();
	if (na != nb) vector_error();
	#endif
	double dot_product = 0.0;
	for(int i=0; i<na; ++i) dot_product += a[i]*b[i];

	return dot_product;
}

void add(const vector &a, const vector &b, vector &c)
{
	int na = a.size();
	#ifdef INDEX_CHECK
	int nb = b.size();
	if (na != nb) vector_error();
	#endif
	c.resize(na);

	for(int i=0; i<na; ++i) c[i] = a[i] + b[i];
}

void sub(const vector &a, const vector &b, vector &c)
{
	int na = a.size();
	#ifdef INDEX_CHECK
	int nb = b.size();
	if (na != nb) vector_error();
	#endif
	c.resize(na);

	for(int i=0; i<na; ++i) c[i] = a[i] - b[i];
}

void mul(double scalar, const vector &a, vector &b)
{
	int n = a.size();
	b.resize(n);    
	for(int i=0; i<n; ++i) b[i] = scalar * a[i];
}

void mul(double scalar, vector &a)
{
	int n = a.size();
	for(int i=0; i<n; ++i) a[i] *= scalar;
}

void div(const vector &a, double scalar, vector &b)
{
	int n = a.size();
	b.resize(n);
	for(int i=0; i<n; ++i) b[i] = a[i] / scalar;
}

void div(vector &a, double scalar)
{
	int n = a.size();
	for(int i=0; i<n; ++i) a[i] /= scalar;
}

double dot(const vector &a, const vector &b)
{
	int na = a.size();
	#ifdef INDEX_CHECK
	int nb = b.size();
	if (na != nb) vector_error();
	#endif
	double sum = 0.0;
	for(int i=0; i<na; ++i) sum += a[i] * b[i];

	return sum;
}

void vector_error()
{
	std::cerr << "vector size mismatch" << '\n'; exit(1);
}


// matrix functions

matrix &matrix::operator = (const matrix &rhs)
{
	if( nrows != rhs.nrows ) {
		delete[] m;
		nrows = rhs.nrows;
		m = new vector[nrows];
	}    
	for(int i=0; i<nrows; ++i) m[i] = rhs.m[i];

	return *this;
}

matrix &matrix::operator = (const double &rhs)
{
	for(int i=0; i<nrows; ++i) m[i] = rhs;
	
	return *this;
}

void matrix::index_error() const
{
	std::cerr << "matrix index out of range" << '\n'; exit(1);
}

std::ostream& operator <<(std::ostream& os, const matrix &a)
{    
	int nrows = a.size_row(), ncols = a.size_col();
	for(int i=0; i<nrows; ++i) {
		for(int j=0; j<ncols; ++j) os << a[i][j] << '\t'; 
		os << '\n';
	}

	return os;
}

std::istream& operator >>(std::istream& is, matrix &a)
{
	int nrows = a.size_row(), ncols = a.size_col();
	for(int i=0; i<nrows; ++i) for(int j=0; j<ncols; ++j) is >> a[i][j];

	return is;
}	

matrix operator + (const matrix &a, const matrix &b)
{
	int nrowsa = a.size_row(), ncolsa = a.size_col();
	#ifdef INDEX_CHECK
	int nrowsb = b.size_row(), ncolsb = b.size_col();
	if(nrowsa != nrowsb || ncolsa != ncolsb) matrix_error();
	#endif	
	matrix c(nrowsa, ncolsa);
	for(int i=0; i<nrowsa; ++i) for(int j=0; j<ncolsa; ++j) c[i][j] = a[i][j] + b[i][j];

	return c;
}

matrix operator - (const matrix &a, const matrix &b)
{
	int nrowsa = a.size_row(), ncolsa = a.size_col();
	#ifdef INDEX_CHECK
	int nrowsb = b.size_row(), ncolsb = b.size_col();
	if(nrowsa != nrowsb || ncolsa != ncolsb) matrix_error();
	#endif
	matrix c(nrowsa, ncolsa);
	for(int i=0; i<nrowsa; ++i) for(int j=0; j<ncolsa; ++j) c[i][j] = a[i][j] - b[i][j];

	return c;
}

matrix operator * (double scalar, const matrix &a)
{
	int nrows = a.size_row(), ncols = a.size_col(); 
	matrix b(nrows, ncols);
	for(int i=0; i<nrows; ++i) for(int j=0; j<ncols; ++j) b[i][j] = scalar * a[i][j];

	return b;
}

matrix operator * (const matrix &a, double scalar)
{
	int nrows = a.size_row(), ncols = a.size_col(); 
	matrix b(nrows, ncols);
	for(int i=0; i<nrows; ++i) for(int j=0; j<ncols; ++j) b[i][j] = scalar * a[i][j];

	return b;
}

vector operator * (const matrix &m, const vector &v)
{
	int nrowm = m.size_row(), ncolm = m.size_col();
	#ifdef INDEX_CHECK
	int nv = v.size();
	if( ncolm != nv ) matrix_error();
	#endif
	vector u(nrowm, 0);
	for(int i=0; i<nrowm; ++i) for(int j=0; j<ncolm; ++j) u[i] += m[i][j] * v[j];
	
	return u;
}

vector operator * (const vector &v, const matrix &m)
{
	int nrowm = m.size_row(), ncolm = m.size_col();
	#ifdef INDEX_CHECK
	int nv = v.size();
	if( nrowm != nv ) matrix_error();
	#endif
	vector u(ncolm, 0);
	for(int j=0; j<ncolm; ++j) for(int i=0; i<nrowm; ++i) u[j] += v[i] * m[i][j];
	
	return u;
}

matrix operator / (const matrix &a, double scalar)
{
	int nrows = a.size_row(); int ncols = a.size_col(); 
	matrix b(nrows, ncols);
	for(int i=0; i<nrows; ++i) for(int j=0; j<ncols; ++j) b[i][j] = a[i][j] / scalar;

	return b;
}

matrix &operator += (matrix &a, const matrix &b)
{
	int nrowsa = a.size_row(); int ncolsa = a.size_col();
	#ifdef INDEX_CHECK
	int nrowsb = b.size_row(); int ncolsb = b.size_col();
	if(nrowsa != nrowsb || ncolsa != ncolsb) matrix_error();
	#endif
	for(int i=0; i<nrowsa; ++i) for(int j=0; j<ncolsa; ++j) a[i][j] += b[i][j];

	return a;
}

matrix &operator -= (matrix &a, const matrix &b)
{
	int nrowsa = a.size_row(); int ncolsa = a.size_col();
	#ifdef INDEX_CHECK
	int nrowsb = b.size_row(); int ncolsb = b.size_col();
	if(nrowsa != nrowsb || ncolsa != ncolsb) matrix_error();
	#endif
	for(int i=0; i<nrowsa; ++i) for(int j=0; j<ncolsa; ++j) a[i][j] -= b[i][j];

	return a;
}

matrix &operator *= (matrix &a, double scalar)
{
	int nrows = a.size_row(); int ncols = a.size_col();
	for(int i=0; i<nrows; ++i) for(int j=0; j<ncols; ++j) a[i][j] *= scalar; 

	return a;
}

matrix &operator /= (matrix &a, double scalar)
{
	int nrows = a.size_row(); int ncols = a.size_col();
	for(int i=0; i<nrows; ++i) for(int j=0; j<ncols; ++j) a[i][j] /= scalar; 

	return a;
}

void add(const matrix &a, const matrix &b, matrix &c)
{
	int nrowsa = a.size_row(); int ncolsa = a.size_col();
	#ifdef INDEX_CHECK
	int nrowsb = b.size_row(); int ncolsb = b.size_col();
	if(nrowsa != nrowsb || ncolsa != ncolsb) matrix_error();
	#endif
	c.resize(nrowsa, ncolsa);

	for(int i=0; i<nrowsa; ++i) for(int j=0; j<ncolsa; ++j) c[i][j] = a[i][j] + b[i][j];
}

void sub(const matrix &a, const matrix &b, matrix &c)
{
	int nrowsa = a.size_row(); int ncolsa = a.size_col();
	#ifdef INDEX_CHECK
	int nrowsb = b.size_row(); int ncolsb = b.size_col();
	if(nrowsa != nrowsb || ncolsa != ncolsb) matrix_error();
	#endif
	c.resize(nrowsa, ncolsa);
	for(int i=0; i<nrowsa; ++i) for(int j=0; j<ncolsa; ++j) c[i][j] = a[i][j] - b[i][j];
}

void mul(double scalar, const matrix &a, matrix &b)
{
	int nrows = a.size_row(); int ncols = a.size_col();
	b.resize(nrows, ncols);    
	for(int i=0; i<nrows; ++i) for(int j=0; j<ncols; ++j) b[i][j] = scalar * a[i][j];
}

void mul(double scalar, matrix &a)
{
	int nrows = a.size_row(); int ncols = a.size_col();
	for(int i=0; i<nrows; ++i) for(int j=0; j<ncols; ++j) a[i][j] *= scalar; 
}

void mul(const matrix &m, const vector &va, vector &vb)
{
	int nrows = m.size_row(); int ncols = m.size_col();
	#ifdef INDEX_CHECK
	if( ncols != va.size() ) matrix_error();
	#endif	
	vb = 0;
	for(int i=0; i<nrows; ++i) for(int j=0; j<ncols; ++j) vb[i] += m[i][j] * va[j]; 
}

void mul(const vector &va, const matrix &m, vector &vb)
{
	int nrows = m.size_row(); int ncols = m.size_col();
	#ifdef INDEX_CHECK
	if( nrows != va.size() ) matrix_error();
	#endif
	vb = 0;
	for(int j=0; j<ncols; ++j) for(int i=0; i<nrows; ++i) vb[j] += m[i][j] * va[i]; 
}

void div(const matrix &a, double scalar, matrix &b)
{
	int nrows = a.size_row(); int ncols = a.size_col();
	b.resize(nrows, ncols);    
	for(int i=0; i<nrows; ++i) for(int j=0; j<ncols; ++j) b[i][j] = a[i][j] / scalar;
}

void div(matrix &a, double scalar)
{
	int nrows = a.size_row(); int ncols = a.size_col();
	for(int i=0; i<nrows; ++i) for(int j=0; j<ncols; ++j) a[i][j] /= scalar; 
}

void transpose(const matrix &a, matrix &aTrans)
{
	int nrowsTrans = a.size_col();
	int ncolsTrans = a.size_row();

	aTrans.resize(nrowsTrans, ncolsTrans);

	for(int i=0; i<nrowsTrans; ++i)
		for(int j=0; j<ncolsTrans; ++j)
			aTrans[i][j] = a[j][i];
}

void transpose(matrix &a)
{
	int nrowsTrans = a.size_col();
	int ncolsTrans = a.size_row();

	matrix aTrans(nrowsTrans, ncolsTrans);
	for(int i=0; i<nrowsTrans; ++i)
		for(int j=0; j<ncolsTrans; ++j)
			aTrans[i][j] = a[j][i];

	a = aTrans;
}

double trace(const matrix &m)
{
	#ifdef INDEX_CHECK
	if( m.size_row() != m.size_col() ) { std::cerr << "matrix trace : not a square matrix\n"; exit(1); }
	#endif
	int n = m.size();
	double tr = 0;

	for(int i=0; i<n; ++i) tr += m[i][i];
	
	return tr;
}

double norm( const matrix &m )
{
	double square_sum = 0.;
	for( int i=0; i<m.size_row(); ++i ) for( int j=0; j<m.size_col(); ++j )
		square_sum += m[i][j]*m[i][j];
	
	return sqrt( square_sum );
}

void matrix_error()
{
	std::cerr << "matrix size mismatch\n"; exit(1);
}
