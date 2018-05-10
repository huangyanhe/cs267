#include "ParticleSet.H"

void ParticleShift::init(const ParticleSet& a_particles)
{
  int sizeOfa_particle = a_particles.m_particles.size();
  m_particles.resize(sizeOfa_particle);
  for (int j = 0; j<m_particles.size(); j++)
    {
      //Sets EField to zero
      m_particles[j] *= 0.0;
    }
}
void ParticleShift::initDelta(const ParticleSet& a_particles)
{
  int sizeOfa_particle = a_particles.m_particles.size();
  m_particles.resize(sizeOfa_particle);
  for (int dir =0; dir<DIM;dir++)
    {
      for (int j = 0; j<m_particles.size(); j++)
	{
	  m_particles[j].m_x[dir] *= 0.0;
	  m_particles[j].m_v[dir] *= 0.0;
	}
    }
}
  /// m_particles[k] += a_rhs.m_particles[k]*a_scale.
void ParticleShift::increment(
                 double a_scale, 
                 const ParticleShift& a_rhs)
{
  for (int j = 0; j<a_rhs.m_particles.size(); j++)
    {
      m_particles[j].increment(a_scale, a_rhs.m_particles[j]);
    }
}
void ParticleShift::incrementVelocity(
				      double a_scale,
				      const ParticleShift& a_rhs)
{
 for (int j = 0; j<a_rhs.m_particles.size(); j++)
    {
      m_particles[j].incrementVelocity(a_scale, a_rhs.m_particles[j]);
    } 
}
void ParticleShift::incrementPositionandVelocity(
						 double a_scaleX,
						 double a_scaleV,
						 const ParticleShift& a_rhs)
{
 for (int j = 0; j<a_rhs.m_particles.size(); j++)
    {
      m_particles[j].incrementPositionandVelocity(a_scaleX,  a_scaleV, a_rhs.m_particles[j]);
    } 
}
  /// m_particles[k] *= a_scale
void ParticleShift::operator*=(double a_scale)
{
  for (int j = 0; j<m_particles.size(); j++)
    {
      m_particles[j] *= a_scale;
    }
}
  /// reinitializes the values m_particles[k] to zero. Not used in RK4.
void ParticleShift::setToZero()
{
  for (int j = 0; j<m_particles.size(); j++)
    {
      m_particles[j] *= 0.0;
    }
}
void ParticleShift::zeroEField()
{
  for (int j = 0; j<m_particles.size(); j++)
    {
      for (int dir=0; dir<DIM; dir++)
	{
	  m_particles[j].EField[dir] *= 0.0;
	}
    }
}
