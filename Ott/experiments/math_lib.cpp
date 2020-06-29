#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <ctime>
#include <algorithm>
#include <limits>
#include "math_lib.h"
#include "cblas.h"
#include "f2c.h"
#include "clapack.h"

const double PI  = 3.14159265358979323846264338327950288419716939937510;

// random number generators

double random_uniform()
// Generate a random number between (0,1)
{
	static int idum = -9;
	//static int idum = -int( time(0) );
	return NR::ran1( idum );
}

double random_gaussian( double sigma )
// Generate a random number from a Gaussian distribution of sigma with 0 mean.
{
	static int idum = -9;
	//static int idum = -int( time(0) );
	double x, y, r2;

	do {
		//choose x,y in a uniform square x:(-1,1), y:(-1,1)
		x = 2.0*NR::ran1( idum ) - 1.0;
		y = 2.0*NR::ran1( idum ) - 1.0;
		//see if it is in the unit circle
		r2 = x * x + y * y;
	}
	while(r2 > 1.0 || r2 == 0);
	
	//Box-Muller transform
	return sigma * x * sqrt (-2.0 * log (r2) / r2);
}

// matrix functions

matrix operator * ( const matrix &a, const matrix &b )
{
	matrix c( a.size_row(), b.size_col() );
	mul( a, b, c );
	return c;
}

#ifdef LAPACK

void mul( const matrix &a, const matrix &b, matrix &c )
{
	int i, j;
	int arow = a.size_row(), acol = a.size_col();
	int brow = b.size_row(), bcol = b.size_col();
	int crow = c.size_row(), ccol = c.size_col();

	#ifdef INDEX_CHECK
	if( acol != brow || crow != arow || ccol != bcol ) { std::cerr << "mul : matrix mismatch"; exit(1); }
	#endif

	double *aL = new double[arow*acol];
	double *bL = new double[brow*bcol];
	double *cL = new double[crow*ccol];
	
	for( j=0; j<acol; ++j ) for( i=0; i<arow; ++i ) aL[i+arow*j] = a(i,j);
	for( j=0; j<bcol; ++j ) for( i=0; i<brow; ++i ) bL[i+brow*j] = b(i,j);

	CBLAS_ORDER Order = CblasColMajor;
	CBLAS_TRANSPOSE TransA = CblasNoTrans;
	CBLAS_TRANSPOSE TransB = CblasNoTrans;
	double alpha = 1., beta = 0.;

	cblas_dgemm( Order, TransA, TransB, arow, bcol, acol, alpha, aL, arow, bL, brow, beta, cL, crow );

	for( i=0; i<crow; ++i ) for( j=0; j<ccol; ++j ) c(i,j) = cL[i+crow*j];

	delete [] aL; delete [] bL; delete[] cL;
}

int inverse( const matrix &a, matrix &aInv )
{
	aInv = a;
	return inverse(aInv);
}

int inverse( matrix &a )
{
	integer n = a.size_row(), info, lwork = -1, *ipiv = new integer[n];
	double lwork_op, *work;
	
	int i, j;
	double *aL = new double[n*n];
	for( j=0; j<n; ++j ) for( i=0; i<n; ++i ) aL[i+n*j] = a(i,j);
	
	dgetrf_( &n, &n, aL, &n, ipiv, &info);
	
	if( info != 0 ) return info;
	
	// inquire about optimal lwork
	dgetri_( &n, aL, &n, ipiv, &lwork_op, &lwork, &info );
	lwork = integer( lwork_op + 0.5 );
	work = new double[lwork];

	dgetri_( &n, aL, &n, ipiv, work, &lwork, &info );

	for( i=0; i<n; ++i ) for( j=0; j<n; ++j ) a(i,j) = aL[i+n*j];

	delete [] work; delete[] ipiv; delete [] aL;
	
	return info;
}

int eigen( const matrix &a, vector &eigenvalues, matrix &eigenvectors )
{
	eigenvectors = a;
	return eigen( eigenvectors, eigenvalues );
}

int eigen( matrix &a, vector &eigenvalues )
{
	integer n = eigenvalues.size(), lwork = -1, info;
	char jobz = 'v', uplo = 'u';
	double lwork_op;

	int i, j;
	double *evL = new double[n];
	double *aL = new double[n*n];
	for( j=0; j<n; ++j ) for( i=0; i<n; ++i ) aL[i+n*j] = a(i,j);
		
	// inquire about optimal lwork
	dsyev_( &jobz, &uplo, &n, aL, &n, evL, &lwork_op, &lwork, &info);
	lwork = integer(lwork_op + 0.5);
	double *work = new double[lwork];

	dsyev_( &jobz, &uplo, &n, aL, &n, evL, work, &lwork, &info);
	
	for( i=0; i<n; ++i ) for( j=0; j<n; ++j ) a(i,j) = aL[i+n*j];
	for( i=0; i<n; ++i ) eigenvalues(i) = evL[i];
	delete [] work; delete[] evL; delete[] aL;
	
	return info;
}

#else

void mul(const matrix &a, const matrix &b, matrix &c)
{
	#ifdef INDEX_CHECK
	if( a.size_col() != b.size_row() ) matrix_error();
	#endif	
	int ni = a.size_row();
	int nj = b.size_col();
	int nk = a.size_col();
	double sum;

	c.resize(ni, nj); 

	for(int i=0; i<ni; ++i) {
		for(int j=0; j<nj; ++j) {
			sum = 0;
			for(int k=0; k<nk; ++k) sum += a[i][k]*b[k][j];
			c[i][j] = sum;
		}
	}
}

int inverse( const matrix &a, matrix &aInv )
{
	aInv = a;
	return inverse(aInv);
}

int inverse(matrix &a)
{	
	matrix b;
	return NR::gaussj(a, b);
}

int eigen( const matrix &a, vector &eigenvalues, matrix &eigenvectors )
{
	eigenvectors = a;
	return eigen( eigenvectors, eigenvalues );
}

int eigen( matrix &a, vector &eigenvalues )
// input: a, output: a, eigenvalues. a is replace by the matrix with columns of normalized 
// eigenvectors of corresponding eigenvalues
{
	int info;
	eigenvalues.resize(a.size_row());
	vector e(a.size_row());
	NR::tred2(a, eigenvalues, e);
	info = NR::tqli(eigenvalues, e, a);

	int nRow = a.size_row(), nCol = a.size_col();
	bool flag;	
	do {
		flag = false;
		for(int j=0; j<nCol-1; ++j)
		{
			if( eigenvalues[j] > eigenvalues[j+1] ) 
			{
				NR::SWAP( eigenvalues[j], eigenvalues[j+1] );
				for(int i=0; i<nRow; ++i)
				{
					NR::SWAP( a[i][j], a[i][j+1] );
				}
				flag = true;
			}
		}
	} while( flag );
	
	return info;
}

#endif

// other

// find a Fourier transform of 1-dim data, save real and imaginary values only for positive frequencies
void DFT( const vector &x, vector &real, vector &imag, int kmin, int kmax )
{
	double theta;
	int N = x.size();
	if ( kmax == 0 ) kmax = N/2 - 1;
	
	real.resize( N/2 );
	imag.resize( N/2 );
	real = imag = 0.;
	
	for( int k = kmin; k <= kmax; ++k )
	{
		for( int n = 0; n < N; ++n )
		{
			theta = 2*PI*k*n/N;
			real[k] += x[n] * cos( theta );
			imag[k] -= x[n] * sin( theta );
		}
	}
}

// find a Fourier transform of 2-dim data, save real and imaginary values only for positive frequencies
void DFT( const matrix &z, matrix &real, matrix &imag, int kxmin, int kxmax, int kymin, int kymax )
{
	double theta;
	int Nx = z.size_row(), Ny = z.size_col();
	if( kxmax == 0 ) kxmax = Nx/2 - 1;
	if( kymax == 0 ) kymax = Ny/2 - 1;
	
	real.resize( Nx/2, Ny/2 );
	imag.resize( Nx/2, Ny/2 );
	real = imag = 0.;
	
	for( int kx = kxmin; kx <= kxmax; ++kx ) for( int ky = kymin; ky <= kymax; ++ky )
	{
		for( int nx = 0; nx < Nx; ++nx ) for( int ny = 0; ny < Ny; ++ny )
		{
			theta = 2*PI*( double(kx*nx)/Nx + double(ky*ny)/Ny );
			real[kx][ky] += z[nx][ny] * cos( theta );
			imag[kx][ky] -= z[nx][ny] * sin( theta );
		}
	}
}

void envelope( const vector &RealFT, const vector &ImagFT, vector &env, int kmin, int kmax )
{
	double theta, xReal, xImag;
	int N = env.size();
	if( kmax == 0 ) kmax = N/2 - 1;
	
	// build amplitude from inverse Fourier transform only with positive k's
	for( int n = 0; n < N; ++n )
	{
		xReal = xImag = 0.;
		for( int k = kmin; k <= kmax; ++k )
		{
			theta = 2*PI*k*n/N;
			xReal += RealFT[k] * cos( theta ) - ImagFT[k] * sin( theta );
			xImag += ImagFT[k] * cos( theta ) + RealFT[k] * sin( theta );
		}
		xReal *= 2; xImag *= 2;
		xReal /= N; xImag /= N;
		env[n] = sqrt( xReal*xReal + xImag*xImag );
	}
}

void envelope( const matrix &RealFT, const matrix &ImagFT, matrix &env, int kxmin, int kxmax, int kymin, int kymax )
{
	double theta, zReal, zImag;
	int Nx = env.size_row(), Ny = env.size_col();
	if( kxmax == 0 ) kxmax = Nx/2 - 1;
	if( kymax == 0 ) kymax = Ny/2 - 1;
	
	// build amplitude from inverse Fourier transform only with positive k's
	for( int nx = 0; nx < Nx; ++nx ) for( int ny = 0; ny < Ny; ++ny )
	{
		zReal = zImag = 0.;
		for( int kx = kxmin; kx <= kxmax; ++kx ) for( int ky = kymin; ky <= kymax; ++ky )
		{
			theta = 2*PI*( double(kx*nx)/Nx + double(ky*ny)/Ny );
			zReal += RealFT[kx][ky] * cos( theta ) - ImagFT[kx][ky] * sin( theta );
			zImag += ImagFT[kx][ky] * cos( theta ) + RealFT[kx][ky] * sin( theta );
		}
		zReal *= 2; zImag *= 2;
		zReal /= Nx*Ny; zImag /= Nx*Ny;
		env[nx][ny] = sqrt( zReal*zReal + zImag*zImag );
	}
}

void choose_regular_numbers( int n, int N, std::vector<int> &numbers )
{
	double interval = double(N)/n;
	numbers.resize(n);
	
	for( int i=0; i<n; ++i ) numbers[i] = int(i * interval + 0.5);
}

void choose_random_numbers( int n, int N, std::vector<int> &numbers )
{
	numbers.resize(n);
	int *location_temp = new int[N];
	for(int i=0; i<N; ++i) location_temp[i] = i;
	for(int j, i=0; i<n; ++i) {
		j = int( (N-i) * random_uniform() );
		numbers[i] = location_temp[j];
		location_temp[j] = location_temp[N-i-1];
	}
	delete [] location_temp;
	std::sort( numbers.begin(), numbers.end() );
}

double find_root( double (*func)( double ), double (*funcd)( double ), double x1, double x2, double accuracy )
// rtsafe of NR. Using a combination of Newton-Raphson and bisection, find the root of a function bracketed between x1 and
// x2. The root, returned as the function value, will be refined until its accrarcy is known within plus minus accuracy.
// func and funcd are user-supplied routines that return the function value and the first derivative of the function.
{
	const int MAXIT = 100; // maximum allowed number of iterations
	int j;
	double df, dx, dxold, f, fh, fl;
	double temp, xh, xl, rts;
	
	fl = (*func)( x1 );
	fh = (*func)( x2 );
	if( fl > 0. && fh > 0. || fl < 0. && fh < 0. ) std::cerr << "Root must be bracketed in rtsafe.";
	if( fl == 0. ) return x1;
	if( fh == 0. ) return x2;
	if( fl < 0. ) // Orient the search so that f(xl) < 0.
	{
		xl = x1; xh = x2;
	}
	else
	{
		xh = x1; xl = x2;
	}
	
	rts = 0.5 * (x1+x2);   // Initialize the guess for root,
	dxold = fabs( x2-x1 ); // the "stepsize before last",
	dx = dxold;			   // and the last step.
	f = (*func)( rts );
	df = (*funcd)( rts );
	
	for( j = 1; j <= MAXIT; ++j ) // Loop over allowed iterations.
	{
		if( ((rts-xh)*df-f )*((rts-xl)*df-f) > 0. || fabs(2.*f) > fabs(dxold*df) ) // Bisect if Newton out of range,
		{																		   // or not decreasing fast enough.
			dxold = dx;
			dx = 0.5 * (xh-xl);
			rts = xl + dx;
			if( xl == rts ) return rts; // Change in root is negligible.
		}
		else // Newton step acceptable. Take it.
		{
			dxold = dx;
			dx = f/df;
			temp = rts;
			rts -= dx;
			if( temp == rts ) return rts;
		}
		if( fabs(dx) < accuracy ) return rts; // Convergence criterion.
		
		f = (*func)( rts );
		df = (*funcd)( rts );
		
		if( f < 0. ) xl = rts; // Maintain the bracket on the root. 
		else xh = rts;
	}
	
	std::cerr << "Maximum number of iterations exceeded in rtsafe";
	return 0.; // Never get here.
}

int downhill_simplex( double (*func)( const vector & ), matrix &vertex, double termination )
{
	int dim = vertex.size_col(), nvertex = vertex.size_row(), nstep;
	int i_low, i_shigh, i_high;
	vector f_vertex( nvertex ), cen( dim ), ref( dim ), exp( dim ), con_side( dim ), con( dim ), mean( dim );
	double f_low, f_shigh, f_high, f_ref, f_exp, f_con_side, f_con, rms;
	
	// calculate function values at vertices
	for( int i = 0; i < nvertex; ++i ) f_vertex[i] = func( vertex[i] );
	
	for( nstep = 0; ; ++nstep )
	{
		// find lowest, second highest, and highest vertices
		if( f_vertex[0] < f_vertex[1] )
		{
			i_low = 0; i_shigh = 0;	i_high = 1;
		}
		else
		{
			i_low = 1; i_shigh = 1; i_high = 0;
		}

		for( int i = 2; i < nvertex; ++i )
		{
			if( f_vertex[i] <= f_vertex[i_low] ) i_low = i;
			
			if( f_vertex[i] > f_vertex[i_high] ) { i_shigh = i_high; i_high = i; }
			else if( f_vertex[i] > f_vertex[i_shigh] ) i_shigh = i;
		}

		f_low   = func( vertex[ i_low   ] );
		f_shigh = func( vertex[ i_shigh ] );
		f_high  = func( vertex[ i_high  ] );

		// check termination condition
		mean = 0.;
		for( int i = 0; i < nvertex; ++i ) mean += vertex[i];
		mean /= nvertex;
		
		rms = 0.;
		for( int i = 0; i < nvertex; ++i ) for( int j = 0; j < dim; ++j )
			rms += ( vertex[i][j] - mean[j] ) * ( vertex[i][j] - mean[j] );
		rms = sqrt( rms / ( nvertex * dim ) );

		if( rms < termination )
		{
			mean = vertex[0];
			vertex[0] = vertex[i_low];
			vertex[i_low] = mean;
			break;
		}
		
		// calculate centroid
		cen = 0.;
		for( int i = 0; i < nvertex; ++i )
			if( i != i_high ) cen += vertex[i];
		cen /= nvertex - 1;
		
		// reflect
		ref = cen + cen - vertex[ i_high ];
		f_ref = func( ref );
		
		if( f_ref >= f_low && f_ref < f_shigh )
		{
			// accept reflection
			vertex[ i_high ] = ref;
			f_vertex[ i_high ] = f_ref;
		}
		else if( f_ref < f_low )
		{
			// expand
			exp = cen + 2 * ( ref - cen );
			f_exp = func( exp );
			if( f_exp < f_ref )
			{	// accept expansion
				vertex[ i_high ] = exp;
				f_vertex[ i_high ] = f_exp;
			}
			else 
			{	// accept reflection
				vertex[ i_high ] = ref;
				f_vertex[ i_high ] = f_ref;
			}
		}
		else
		{	// contract
			con_side = f_ref < f_high ? ref : vertex[ i_high ];
			f_con_side = f_ref < f_high ? f_ref : f_high;
			
			con = cen + ( con_side - cen )/2;
			f_con = func( con );
			
			if( f_con < f_con_side )
			{	
				// accept contraction
				vertex[ i_high ] = con;
				f_vertex[ i_high ] = f_con;
			}
			else
			{	// shrink
				for( int i = 0; i < nvertex; ++i )
				{
					if( i != i_low )
					{
						vertex[i] = vertex[i_low] + ( vertex[i] - vertex[i_low] ) / 2;
						f_vertex[i] = func( vertex[i] );
					}
				}
			}
		}
	}
	
	return nstep;
}

// numerical recipe

inline double NR::pythag(const double &a, const double &b)
//computes (a^2 + b^2)^(1/2) without destructive underflow or overflow
{
	double absa = fabs(a), absb = fabs(b);
	if( absa > absb ) return absa * sqrt( 1.0 + SQR(absb/absa) );
	else return ( absb == 0.0 ? 0.0 : absb * sqrt( 1.0 + SQR(absa/absb) ) );
}

int NR::gaussj(matrix &a, matrix &b)
// Gauss-Jordan elimination with full pivoting (Numerical Recipes in C++, 2nd edition)
// Solves a x = b. a: n by n matrix, x and b: n by m matrix. a is replaced by its inverse matrix.
// b is replaced by the solution x.
{
	int i, icol, irow, j, k, l, ll;
	double big, dum, pivinv;

	int n = a.size_row();
	int m = b.size_col();
	
	//These integer arrays are used for bookkeeping on the pivoting
	int *indxc = new int[n], *indxr = new int[n], *ipiv = new int[n];
	for(j=0; j<n; ++j) ipiv[j] = 0;

	//This is the main loop over the columns to be reduced.
	for(i=0; i<n; ++i) {
		big = 0.0;
		//This is the outer loop of the search for a pivot element
		for(j=0; j<n; ++j) {
			if (ipiv[j] != 1) {
				for(k=0; k<n; ++k) {
					if(ipiv[k] == 0) {
						if(fabs(a(j,k)) >= big) {
							big = fabs(a(j,k));
							irow = j;
							icol = k;
						}
					}
				}
			}
		}
		++(ipiv[icol]);

		// We now have the pivot element, so we interchange rows, if needed, to put the pivot element on the diagonal.
		// The columns are not physically interchanged, only relabeled. 
		// indxc[i], the column of the (i+1)th pivot element, is the (i+1)th column that is reduced, 
		// while indxr[i] is the row in which that pivot element was originally located. 
		// If indxr[i]!=indxc[i] there is an implied column interchange.
		// With this form of bookkeeping, the solution b's will end up in the correct order,
		// and the inverse matrix will be scrambled by columns.

		if(irow != icol) {
			for(l=0; l<n; ++l) SWAP(a(irow,l), a(icol,l));
			for(l=0; l<m; ++l) SWAP(b(irow,l), b(icol,l));
		}

		// We are now ready to divide the pivot row by the pivot element, located at irow and icol.
		indxr[i] = irow;
		indxc[i] = icol;
		if( a(icol,icol) == 0.0 ) { std::cerr << "GaussJordan: Singular matrix" << '\n'; return -2; }
		pivinv = 1.0/a(icol,icol);
		a(icol,icol) = 1.0;
		for(l=0; l<n; ++l) a(icol,l) *= pivinv;
		for(l=0; l<m; ++l) b(icol,l) *= pivinv;
		
		//Next we reduce the rows
		for(ll=0; ll<n; ++ll) {			
			if(ll != icol) { //...except for the pivot one, of course.
				dum = a(ll,icol);
				a(ll,icol) = 0.0;
				for(l=0; l<n; ++l) a(ll,l) -= a(icol,l) * dum;
				for(l=0; l<m; ++l) b(ll,l) -= b(icol,l) * dum;
			}
		}
	}
	
	// This is the end of the main loop over columns of the reduction. It only remains to unscramble the solution 
	// in view of the column interchanges.
	// We do this by interchanging pairs of columns in the reverse order that the permutation was buit up.

	for(l=n-1; l>=0; --l)
		if(indxr[l] != indxc[l])
			for(k=0; k<n; ++k) SWAP(a(k,indxr[l]), a(k,indxc[l]));
			
	delete[] indxc; delete[] indxr; delete[] ipiv;
	
	return 0;
}

void NR::tred2(matrix &a, vector &d, vector &e)
// Householder reduction of a real, symmetric matrix a[0..n-1][0..n-1]. On output, a is replaced by the orthogonal 
// matirx Q effecting the transformation. d[0..n-1] returns the dianonal elements of the tridiagonal matrix, 
// and e[0..n-1] the off-diagonal elements, with e[0]=0.
// Several statements, as noted in comments, can be omitted if only eigenvalues are to be found, in which case 
// a contains no useful information on output. Otherwise they are to be included.
{
	int l, k, j, i;
	double scale, hh, h, g, f;

	int n = d.size();
	for(i=n-1; i>0; --i) {
		l = i-1;
		h = scale = 0.0;
		if( l>0 ) {
			for(k=0; k<l+1; ++k) scale += fabs(a(i,k));
			if( scale == 0.0 ) e[i] = a(i,l); //Skip transformation
			else {
				for(k=0; k<l+1; ++k) {
					a(i,k) /= scale; //Use scaled a's for transformation
					h += a(i,k) * a(i,k); //Form sigma in h
				}
				f = a(i,l);
				g = ( f>=0.0 ? -sqrt(h) : sqrt(h) );
				e[i] = scale * g;
				h -= f * g; //Now h is in equation (11.2.4)
				a(i,l) = f - g; //Store u in row i of a
				f = 0.0;
				for(j=0; j<l+1; ++j) {
					//Next statement can be omitted if eigenvectors are not wanted
					a(j,i) = a(i,j)/h; //Store u/H in column i of a
					g = 0.0; //Form an element of A.u in g
					for(k=0; k<j+1; ++k) g += a(j,k) * a(i,k);
					for(k=j+1; k<l+1; ++k) g += a(k,j) * a(i,k);
					e[j] = g/h; //Form element of p in temporarily unused element of e
					f += e[j] * a(i,j);
				}
				hh = f/(h+h); //Form K, equaion (11.2.11)
				for(j=0; j<l+1; ++j) { //Form q and store in e overwriting p
					f = a(i,j);
					e[j] = g = e[j] - hh*f;
					for(k=0; k<j+1; ++k) //Reduce a, equation (11.2.13)
						a(j,k) -= (f*e[k] + g*a(i,k));
				}
			}
		} else e[i] = a(i,l);
		d[i] = h;
	}
	//Next statement can be omitted if eigenvectors are not wanted
	d[0] = 0.0;
	e[0] = 0.0;
	//Contents of the loop can be omitted if eigenvectors are not wanted 
	//except for statement d[i] = a[i][i];
	for(i=0; i<n; ++i) { //Begin accumulation of transformation matrices.
		l = i;
		if( d[i] != 0.0 ) { //This block skipped when i=0
			for(j=0; j<l; ++j) {
				g = 0.0;
				//Use u and u/H stored in a to form P.Q
				for(k=0; k<l; ++k)  g += a(i,k)*a(k,j);
				for(k=0; k<l; ++k) a(k,j) -= g*a(k,i);
			}
		}
		d[i] = a(i,i); //This statement remains
		//Reset row and column of a to identity matix for next iteration
		a(i,i) = 1.0;
		for(j=0; j<l; ++j) a(j,i) = a(i,j) = 0.0;
	}
}

int NR::tqli(vector &d, vector &e, matrix &z)
// QL algorithm with implicit shifts, to determine the eigenvalues and eigenvectors of a real, symmetric, tridiagonal
// matrix, or of a real, symmetric matrix previously reduced by tred2 11.2. On input, d[0..n-1] contains the diagonal
// elements of the tridiagonal matrix. On output, it returns the eigenvalues. The vector e[0..n-1] inputs the subdiagonal
// elements of the tridiagonal matirx, with e[0] arbitrary. On output e is destroyed. When finding only the eigenvalues,
// several lined may be omitted, as noted in the comments. If the eigenvectors of tridiagonal matrix are desired, the matrix
// z[0..n-1][0..n-1] is input as the identity matrix. If the eigenvectors of a matrix that has been reduced by tred2 are
// required, then z is input as the matrix output by tred2. In either case, column k of z returns the normalized eigenvector
// corresponding to d[k].
{
	int m, l, iter, i, k;
	double s, r, p, g, f, dd, c, b;

	int n = d.size();
	for(i=1; i<n; ++i) e[i-1] = e[i]; //Convenient to renumber the elements of e
	e[n-1] = 0.0;
	for(l=0; l<n; ++l) {
		iter = 0;
		do {
			//Look for a single small subdiagonal element to split the matirx
			for(m=l; m<n-1; ++m) {
				dd = fabs(d[m]) + fabs(d[m+1]);
				if(fabs(e[m])+dd == dd) break;
			}
			if( m != l) {
				if( iter++ == 1000 ) { std::cerr << "tqli: too many iterations" << '\n'; return -1; }
				g = (d[l+1]-d[l]) / (2.0*e[l]); //Form shift
				r = pythag(g, 1.0);
				g = d[m] - d[l] + e[l]/(g+SIGN(r,g)); //This is dm-ks
				s = c = 1.0;
				p = 0.0;
				//A plane rotation as in the original QL, followed by Givens rotations 
				//to restore tridianonal form
				for(i=m-1; i>=l; --i) {
					f = s*e[i];
					b = c*e[i];
					e[i+1] = (r=pythag(f, g));
					if( r == 0.0 ) { //Recover from underflow
						d[i+1] -= p;
						e[m] = 0.0;
						break;
					}
					s = f/r;
					c = g/r;
					g = d[i+1] - p;
					r = (d[i]-g)*s + 2.0*c*b;
					d[i+1] = g + (p=s*r);
					g = c*r - b;
					// next loop can be omitted if eigenvectors not wanted
					for(k=0; k<n; ++k) { //Form eigenvectors
						f = z(k,i+1);
						z(k,i+1) = s*z(k,i) + c*f;
						z(k,i) = c*z(k,i) - s*f;
					}
				}
				if( r == 0.0 && i >= l ) continue;
				d[l] -= p;
				e[l] = g;
				e[m] = 0.0;
			}
		} while( m != l );
	}
	
	return 0;
}

double NR::ran0( int &idum )
// Minimal random number generator of Park and Miller. Returns a uniform random deviate between 0.0 and 1.0 Set or reset
// idum to any integer value ( except the unlikely value MASK) to initialize the sequence. idum must not be altered between
// calls for successive deviates in a sequence.
// I think the book uses long for idum because it assumes long is 32 bits, but int is 32 bits on a 64-bit computer.
{
	const int IA = 16807, IM = 2147483647, IQ = 127773, IR = 2836, MASK = 123459876;
	const double AM = 1.0/IM;
	
	int k;
	double ans;
	
	idum ^= MASK; // XORing with MASK allows use of zero and other simple bit patterns for idum
	k = idum / IQ;
	idum = IA * ( idum - k*IQ ) - IR*k; // compute idum=(IA*idum) % IM without overflows by Schrage's method
	if( idum < 0 ) idum += IM;
	ans = AM * idum; // convert idum to a floating result
	idum ^= MASK; // Unmask before return
	return ans;
}

double NR::ran1( int &idum )
// Minimal random number generator of Park and Miller with Bays-Durham shuffle and added safeguards.
// Returns a uniform random deviate between 0.0 and 1.0 (exclusive of the endpoint values). Call with idum a negative
// integer to initialize; thereafter, do not alter idum between successive deviates in a sequence.
// RNMX should approximate the largest floating value that is less than 1.
// I think the book uses long for idum because it assumes long is 32 bits, but int is 32 bits on a 64-bit computer.
{
	const int IA = 16807, IM = 2147483647, IQ = 127773, IR = 2836, NTAB = 32, NDIV = 1 + (IM-1) / NTAB;
	const double EPS = 6.0e-17, AM = 1.0/IM, RNMX = 1.0-EPS;
	
	int j, k;
	static int iy = 0;
	static std::vector<int> iv(NTAB);
	double temp;
	
	if( idum <= 0 || !iy ) { // initialize
		if( -idum < 1 ) idum = 1; // be sure to prevent idum = 0
		else idum = -idum;
		for( j = NTAB+7; j>=0; --j ) { // load the shuffle table ( after 8 warm-ups )
			k = idum / IQ;
			idum = IA * ( idum - k * IQ ) - IR * k;
			if( idum < 0 ) idum += IM;
			if( j < NTAB ) iv[j] = idum;
		}
		iy = iv[0];
	}
	k = idum / IQ; // start here when not initializing
	idum = IA * (idum - k*IQ ) - IR*k; // compute idum = (IA*idum) % IM without over flows by Schrage's method
	if( idum < 0 ) idum += IM;
	j = iy / NDIV; // will be in the range 0 .. NTAB-1
	iy = iv[j]; // output previously stored value and refill the shuffle table
	iv[j] = idum;
	if( ( temp = AM * iy ) > RNMX ) return RNMX; // because users don't expect endpoint values
	else return temp;
}

inline void shft3( double &a, double &b, double &c, const double d ) { a = b; b = c; c = d; }

void NR::mnbrak( double &ax, double &bx, double &cx, double &fa, double &fb, double &fc, double func( const double ) )
// Given a function func, and given distinct initial points ax and bx, this routine searches in the downhill direction
// (defined by the function as evaluated at the initial points) and returns new points ax, bx, cx that bracket a minimum
// of the function. Also returned are the function values at the three points, fa, fb, and fc.
{
	const double GOLD = 1.618034, GLIMIT = 100., TINY = 1.0e-20;
	// Here GOLD is the default ratio by which successive intervals are magnified; GLIMIT is the maximum magnification
	// allowed for a parabolic-fit step.
	double ulim, u, r, q, fu;
	
	fa = func( ax );
	fb = func( bx );
	if( fb > fa ) { SWAP( ax, bx ); SWAP( fb, fa ); }
	// Switch roles of a and b so that we can go downhill in the direction from a to b.
	
	cx = bx + GOLD * ( bx - ax ); // first guess for c
	fc = func( cx );
	while( fb > fc ) // Keep returning here until we bracket.
	{
		r = ( bx - ax ) * ( fb - fc ); // Compute u by parabolic extrapolation from a, b, c.
		q = ( bx - cx ) * ( fb - fa ); // TINY is used to prevent any possible division by zero.
		u = bx - ( (bx-cx)*q-(bx-ax)*r ) / ( 2. * SIGN( MAX( fabs(q-r), TINY ), q-r ) );
		ulim = bx + GLIMIT * ( cx - bx );
		// We won't go farther than this, Test various possibilities:
		if( (bx-u)*(u-cx) > 0. ) // Parabolic u is between b and c: try it.
		{
			fu = func( u );
			if( fu < fc ) // Got a minimum between b and c.
			{
				ax = bx; bx = u; fa = fb; fb = fu; return;
			}
			else if( fu > fb ) // Got a minimum between a and u.
			{
				cx = u; fc = fu; return;
			}
			u = cx + GOLD * (cx-bx); // Parabolic fit was no use. Use default magnification.
			fu = func( u );
		}
		else if( (cx-u)*(u-ulim) > 0. ) // Parabolic fit is between c and its allowed limit.
		{
			fu = func( u );
			if( fu < fc )
			{
				shft3( bx, cx, u, cx + GOLD * (cx-bx) );
				shft3( fb, fc, fu, func(u) );
			}
		}
		else if( (u-ulim)*(ulim-cx) >= 0. ) // Limit parabolic u to maximum allowed value.
		{
			u = ulim;
			fu = func( u );
		}
		else // Reject parabolic u, use default magnification.
		{
			u = cx + GOLD * (cx-bx);
			fu = func( u );
		}
		shft3( ax, bx, cx, u ); // Eliminate oldest point and continue.
		shft3( fa, fb, fc, fu );
	}
}

double NR::brent
( const double ax, const double bx, const double cx, double f( const double ), const double tol, double &xmin )
// Given a function f, and given a bracketing triplet of abscissas ax, bx, cx (such that bx is between ax and cx, and
// f(bx) is less than both f(ax) and f(cx), this routine isolates the minimum to a fractional precision of about tol
// using Brent's method. tol should generally be no smaller than the square root of the machine's floating-point precision.
// The abscissa of the minimum is returned as xmin, and the minimum function value is returned as brent, the returned 
// function value.
{
	const int ITMAX = 1000;
	const double CGOLD = 0.3819660;
	const double ZEPS = std::numeric_limits<double>::epsilon() * 1.0e-3;
	// Here ITMAX is the maximum allowed number of iterations; CGOLD is the golden ratio; ZEPS is a small number that
	// protects against trying to achieve fractional accuracy for a minimum that happens to be exactly zero.
	int iter;
	double a, b, d = 0., etemp, fu, fv, fw, fx;
	double p, q, r, tol1, tol2, u, v, w, x, xm;
	double e = 0.; // This will be the distance moved on the step before last.
	
	a = ( ax < cx ? ax : cx ); // a and b must be in ascending order,
	b = ( ax > cx ? ax : cx ); // but input abscissas need not be.
	x = w = v = bx; // initializations...
	fw = fv = fx = f(x);
	
	for( iter = 0; iter < ITMAX; ++iter ) // main program loop
	{
		xm = 0.5 * (a+b);
		tol2 = 2.0 * ( tol1 = tol*fabs(x) + ZEPS );

		if( fabs(x-xm) <= (tol2-0.5*(b-a)) ) // Test for done here.
		{
			xmin = x;
			return fx;
		}

		if( fabs(e) > tol1 ) // Construct a trial parabolic fit.
		{
			r = (x-w)*(fx-fv);
			q = (x-v)*(fx-fw);
			p = (x-v)*q - (x-w)*r;
			q = 2. * (q-r);
			if( q > 0. ) p = -p;
			q = fabs(q);
			etemp = e;
			e = d;
			if( fabs(p) >= fabs( 0.5*q*etemp ) || p <= q*(a-x) || p >= q*(b-x) )
				d = CGOLD * ( e = ( x >= xm ? a-x : b-x ) );
			// The above conditions determine the acceptability of the parabolic fit. Here we take the golden section step
			// into the larger of the segments.
			else
			{
				d = p/q; // Take the parabolic step.
				u = x+d;
				if( u-a < tol2 || b-u < tol2 ) d = SIGN( tol1, xm-x );
			}
		}
		else d = CGOLD * ( e = ( x >= xm ? a-x : b-x ) );
		
		u = ( fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d) );
		fu = f(u);
		// This is the one funciton evalution per iteration.
		if( fu <= fx ) // Now decide what to do with our function evaluation.
		{
			if( u >= x ) a=x; else b=x;
			shft3( v, w, x, u ); // Housekeeping follows:
			shft3( fv, fw, fx, fu );
		}
		else
		{
			if( u < x ) a=u; else b=u;
			if( fu <= fw || w == x )
			{
				v = w; w = u; fv = fw; fw = fu;
			}
			else if( fu <= fv || v == x || v == w )
			{
				v = u; fv = fu;
			}
		} // Done with housekeeping. Back for another iteration.
	}
	
	std::cerr << "Too many iterations in brent\n";
	xmin = x; // Never get here.
	return fx;
}

//Global variables communicate with f1dim.
int ncom;
double (*nrfunc)( const vector & );
//vector *pcom_p, *xicom_p; // This leads to memory leaks.
vector pcom, xicom;

void NR::linmin( vector &p, vector &xi, double &fret, double func( const vector & ) )
// Given an n-dimensional point p[0..n-1] and an n-dimensional direction xi[0..n-1], moves and resets p to where the
// function func(p) takes on a minimum along the direction xi from p, and replaces xi by the actual vector displacement
// that p was moved. Also returns as fret the value of func at the returned location p. This is actually all accomplished
// by calling the routines mnbrak and brent.
// My experience : Try to make the minimum value close to 0 by adding or subtracting a constant for a better accuracy.
// Then multiply the function by a big constant. This has the effect of making NR::mnbrak of better quality, because
// gradient gets too small near the minimum, and this seems to cause poor performance of mnbrak.
{
	int j;
	const double TOL = 1.0e-8; // tolerance passed to brent
	double xx, xmin, fx, fb, fa, bx, ax;
	
	int n = p.size();
	ncom = n; // Define the global variables.
	//pcom_p = new vector( n ); // This leads to memory leaks.
	//xicom_p = new vector( n ); // This leads to memory leaks.
	pcom.resize( n );
	xicom.resize( n );
	
	nrfunc = func;
	//vector &pcom = *pcom_p, &xicom = *xicom_p; // Make aliases to simplify coding. // This leads to memory leaks.
	for( j=0; j<n; ++j )
	{
		pcom[j] = p[j];
		xicom[j] = xi[j];
	}
	// begin : modified initial guess for brackets
	double xi_length = 0., epsilon = 1./n/5.; // epsilon is an expected distance in each dimension
	for( int i = 0; i < n; ++i ) xi_length += xi[i]*xi[i];
	xi_length = sqrt( xi_length );
	ax = 0.; xx = sqrt(n) * epsilon / xi_length;
	// end
	//ax = 0.; xx = 1.; // initial guess for brackets
	mnbrak( ax, xx, bx, fa, fx, fb, f1dim );
	fret = brent( ax, xx, bx, f1dim, TOL, xmin );
	for( j=0; j<n; ++j ) // Construct the vector results to return.
	{
		xi[j] *= xmin;
		p[j] += xi[j];
	}
	//delete xicom_p, pcom_p; // This leads to memory leaks.
}

double NR::f1dim( const double x ) // Must accompany linmin
{
	int j;

	vector xt( ncom );
	//vector &pcom = *pcom_p, &xicom = *xicom_p; // Make aliases to simplify coding. // This leads to memory leaks.
	for( j = 0; j < ncom; ++j ) xt[j] = pcom[j] + x*xicom[j];
	return nrfunc( xt );
}

void NR::frprmn( vector &p, const double ftol, int &iter, double &fret, double func( const vector & ), void dfunc( const vector &, vector & ) )
// Given a starting point p[0..n-1], Fletcher-Reeves-Polak-Ribiere minimization is performed on a function func, using its
// gradient as calculated by a routine dfunc. The convergence tolerance on the function value is input as ftol. Returned
// quantities are p (the location of the minimum), iter (the number of iterations that were performed), and fret (the
// minimum value of the function). The routine linmin is called to perform line minimizations.
{
	//const int ITMAX = 200;
	const int ITMAX = 100000; // modified
	const double EPS = 1.0e-18;
	// Here ITMAX is the maximum allowed number of iterations, while EPS is a small number to rectify the special case of
	// converging to exactly zero function value.
	int j, its;
	double gg, gam, fp, dgg;
	
	int n = p.size();
	vector g(n), h(n), xi(n);
	fp = func( p ); // initializations
	dfunc( p, xi );

	for( j = 0; j < n; ++j )
	{
		g[j] = -xi[j];
		xi[j] = h[j] = g[j];
	}
	
	for( its = 0; its < ITMAX; ++its ) // Loop over iterations.
	{
		iter = its;
		linmin( p, xi, fret, func ); // Next statement is the normal return:
		if( 2.0 * fabs(fret-fp) <= ftol * (fabs(fret)+fabs(fp)+EPS) ) return;
		fp = fret;
		dfunc( p, xi );
		dgg = gg = 0.;
		for( j = 0; j < n; ++j )
		{
			gg += g[j] * g[j];
			// dgg += xi[j] * xi[j];  // This statement for Fletcher-Reeves.
			dgg += ( xi[j] + g[j] ) * xi[j]; // This statement for Polak-Ribiere.
		}
		if( gg == 0. ) return; // Unlikely. If gradiant is exactly zero then we are already done.
		gam = dgg / gg;
		for( j = 0; j < n; ++j )
		{
			g[j] = -xi[j];
			xi[j] = h[j] = g[j] + gam*h[j];
		}
	}
	
	std::cerr << "Too many iterations in frprmn\n";
}
