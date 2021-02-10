/*============================================================


Name: cmatrix.h

Author: Louis Horner

Date: February 9, 2021

Description:
	Here we define the class cmatrix to deal with 2 dimentional
matricies with complex-valued entries. Needed to run fourier::dft().

cmatrix
	Dependencies: cmatrix, cmath

	Constructors:
			matrix (n_rows, n_cols)

	Public members:
		conjugate()
			takes the complex conjugate of each entry
		print()


============================================================*/
class cmatrix
{
protected:
	int n_rows, n_cols;
	double *** M;

public:
	friend class fourier;

	cmatrix (int tn_rows, int tn_cols) : n_rows(tn_rows), n_cols(tn_cols)
	{
		M = new double ** [tn_rows];
		for (int r=0; r<tn_rows; r++)
		{
			M[r] = new double * [tn_cols];
			for (int c=0; c<tn_cols; c++)
			{
				M[r][c] = new double [2];
			}
		}
	}
	cmatrix () {}

	// destructor
	~cmatrix ()
	{
		delete[] M;
	}

	// printing
	void print ()
	{
		for (int row=0; row<n_rows; row++)
		{
			for (int col=0; col<n_cols; col++)
			{
				std::cout << M[row][col][0] << " + " << M[row][col][1] << "i, ";
			}
			std::cout << "\n";
		}
	}

	// initialization for testing
	cmatrix& operator= (double rhs [4][4][2])
	{
		for (int row=0; row<4; row++)
		{
			for (int col=0; col<4; col++)
			{
				for (int c=0; c<2; c++)
				{
					M[row][col][c] = rhs[row][col][c];
				}
			}
		}
		return *this;
	}
	cmatrix& operator= (double rhs [4][1][2])
	{
		for (int row=0; row<4; row++)
		{
			for (int col=0; col<1; col++)
			{
				for (int c=0; c<2; c++)
				{
					M[row][col][c] = rhs[row][col][c];
				}
			}
		}
		return *this;
	}

	// cmatrix multiplication
	// C = A * B, A = *this, B = rhs
	double *** operator* (cmatrix& rhs)
	{
		if (n_cols == rhs.n_rows)
		{
			double *** C = new double ** [n_rows];
			for (int r=0; r<n_rows; r++)
			{
				C[r] = new double * [rhs.n_cols];
				for (int c=0; c<rhs.n_cols; c++)
				{
					C[r][c] = new double [2];
				}
			}

			double sum_re = 0.0;
			double sum_im = 0.0;
			for (int row_a = 0; row_a < n_rows; row_a++)
			{
				for (int col_b = 0; col_b < rhs.n_cols; col_b++)
				{
					for (int col_a = 0; col_a < n_cols; col_a++)
					{
						sum_re += M[row_a][col_a][0]*rhs.M[col_a][col_b][0] - M[row_a][col_a][1]*rhs.M[col_a][col_b][1];
						sum_im += M[row_a][col_a][0]*rhs.M[col_a][col_b][1] + M[row_a][col_a][1]*rhs.M[col_a][col_b][0];
					}
					C[row_a][col_b][0] = sum_re;
					C[row_a][col_b][1] = sum_im;
					sum_re = 0.0;
					sum_im = 0.0;
				}
			}
			return C;
		}
		else
		{
			std::cout << "Error: you tried to multiply two matricies where the numer of rows in the second did not equal the numer of colums in the first\n";
			throw "Error";
		}
	}

	// copy assignment
	cmatrix& operator= (double *** rhs)
	{
		// copying rhs.M to this->M
		for (int row_r=0; row_r<n_rows; row_r++)
		{
			for (int col_r=0; col_r<n_cols; col_r++)
			{
				for (int p=0; p<2; p++)
				{
					M[row_r][col_r][p] = rhs[row_r][col_r][p];
				}
			}
		}

		// deleting rhs
		delete rhs;

		return *this;
	}

	// setting the entries as their complex conjugates
	cmatrix& conjugate ()
	{
		for (int row=0; row<n_rows; row++)
		{
			for (int col=0; col<n_cols; col++)
			{
				M[row][col][1] *= -1.0;
			}
		}
		return *this;
	}
};
