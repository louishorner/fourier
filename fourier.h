/*============================================================


Name: fourier.h

Author: Louis Horner

Date: February 9, 2021

Description:
	Here we define the class fourier to perform the fourier
tranform of the input function sampler.

fourier
	Dependencies: cmatrix, cmath

	Constructors:
	 		fourier (sampler, unsigned int N, double T)
	 		fourier (sampler, double 2p, unsigned int N)

	Public members:
		dft()
			prints the frequency spectrum


============================================================*/
class fourier
{
protected:
	// input function
	double (*sampler)(double);
	unsigned int N;
	double T;
	double p;
	double omega;

public:
	fourier (
		double (*cf)(double),
		unsigned int cN, 
		double cT
	) : sampler(cf), N(cN), T(cT)
	{}

	fourier (
		double (*cf)(double),
		double a,
		unsigned int cN
	) : sampler(cf), N(cN)
	{
		omega = 2.0*M_PI/a;
		T = a/((double)N);
		p = a/2.0;
	}

	void dft ()
	// discrete fourier transform
	{
		cmatrix F (N, N);
		cmatrix f (N, 1);
		cmatrix c (N, 1);

		// assembling f
		for (int n=0; n<N; n++)
		{
			f.M[n][0][0] = (*sampler)(((double) n)*T);
		}

		// assembling F
		for (int row=0; row<F.n_rows; row++)
		{
			for (int col=0; col<F.n_cols; col++)
			{
				F.M[row][col][0] = std::cos(1.0/((double)N)*2.0*M_PI*((double)row*col));
				F.M[row][col][1] = std::sin(1.0/((double)N)*2.0*M_PI*((double)row*col));
			}
		}

		// final step: C = 1/N * F^bar^ * f
		F.conjugate();
		// multiplying F by 1/N
		for (int row=0; row<F.n_rows; row++)
		{
			for (int col=0; col<F.n_cols; col++)
			{
				for (int c=0; c<2; c++)
				{
					if (F.M[row][col][c] != 0.0)
					{
						F.M[row][col][c] /= N;
					}
				}
			}
		}
		c = F*f;

		// printing the frequencies and aplitudes:
		for (int n=0; n<N; n++)
		{
			std::cout << "f: " << omega*n << ", A: " << std::pow(std::pow(c.M[n][0][0],2) + std::pow(c.M[n][0][1],2), 0.5) << "\n";
		}
	}
};
