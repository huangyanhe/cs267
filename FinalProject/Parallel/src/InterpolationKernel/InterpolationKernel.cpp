#include "InterpolationKernel.H"

InterpolationKernel::InterpolationKernel():
  m_supportSize{0},
  m_coefficients{},
  m_powers{}
{}
InterpolationKernel::~InterpolationKernel()
{}
InterpolationKernel::InterpolationKernel(int a_Order, int a_Smoothness)
{
  if (a_Order == 2)
    {
      if (a_Smoothness == 0)
	{
	  m_supportSize = 1;
	  //m_coefficients(m_supportSize, vector<double>(0,0));
	  //m_coefficients(m_powers, vector<int>(0,0));
	  vector<double> coefficients0;
	  vector<int> powers0;
	  vector<double> coefficients1;
	  vector<int> powers1;
	  coefficients0.push_back(1.0);
	  powers0.push_back(0.0);
	  coefficients0.push_back(-1.0);
	  powers0.push_back(1);
	  coefficients1.push_back(0.0);
	  powers1.push_back(1);
	  

	  m_coefficients.push_back(coefficients0);
	  m_powers.push_back(powers0);
	  m_coefficients.push_back(coefficients1);
	  m_powers.push_back(powers1);
	  
	}
    }
  else if (a_Order == 3)
    {
    }
  else if (a_Order == 4)
    {
    }
  else if (a_Order == 6)
    {
    }
    
}
double InterpolationKernel::apply(double a_x, int supportRegion)
{
  double out = 0.0;
  for (int j = 0; j < m_powers[supportRegion].size(); j++)
    {
      out += m_coefficients[supportRegion][j]*pow(abs(a_x), m_powers[supportRegion][j]);
    }
  return out;
}
int InterpolationKernel::supportSize()
{
  return m_supportSize;
}
