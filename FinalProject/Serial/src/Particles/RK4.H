
#ifndef _RK4_H_
#define _RK4_H_

///  generic explicit RK4 algorithm
/**  class X is the type of your state data. class F provides an evaluation
     of the X' for a specified state of X and the time.  If your X' is 
     determined entirely from the current state values then F can ignore.
     F implements the function:
     void operator(X& k, double a_time, double a_dt, const X& a_state, const X& a_shift);

     F interprets the real state as a_state+a_shift

     Note that F includes the customary factor of dT in the RK4 algorithm. (i.e. this algorithm assumes k_i = F(...) NOT k_i = dt*F(...)
*/
template <class X, class F, class dX>
class RK4
{
public:
  void advance(double a_time, double a_dt, X& a_state);
  void define(X& a_state);
protected:
  //dX is ParticleShift
  // F is ParticleVelocities
  // X is ParticleSet
  dX m_k; 
  dX m_delta;
  F m_f;
};

//==============================================



template <class X, class F, class dX>
void RK4<X, F, dX>::define( X& a_state)
{
  m_f.define(a_state);
}
// void RK4<X, F, dX>::advance(double a_time, double a_dt, X& a_state)
// {
//   // this looks odd, but it is a common programming error
//   //  you can't write "1/6* m_k"; and get what you want.  the first
//   // term is evaluated as integer division, resulting in a zero.
//   double sixth = 1, third=1, half = 1;
//   sixth/=6; third/=3; half/=2;
  
//   m_delta.init(a_state);
//   m_k.init(a_state);
//   m_f(m_k, a_time, a_dt, a_state);  // compute k1
//   m_delta.increment(sixth, m_k);
//   m_k*=half;
//   m_f(m_k, a_time+half*a_dt, a_dt, a_state); // compute k2
//   m_delta.increment(third, m_k);
//   m_k*=half;
//   m_f(m_k, a_time+half*a_dt, a_dt, a_state);  // conpute k3
//   m_delta.increment(third, m_k);
//   m_f(m_k, a_time+a_dt, a_dt, a_state); // compute k4
//   m_delta.increment(sixth, m_k);
//   a_state.increment(m_delta);
// }

//Implementation of Andrews 4th order RK4 scheme
template <class X, class F, class dX>
void RK4<X, F, dX>::advance(double a_time, double a_dt, X& a_state)
{
  // this looks odd, but it is a common programming error
  //  you can't write "1/6* m_k"; and get what you want.  the first
  // term is evaluated as integer division, resulting in a zero.
  double sixth = 1, third=1, half = 1, eigth = 1.0;
  sixth/=6; third/=3; half/=2; eigth/=8.0;

  //m_k are the intermediate ks
  //m_delta is the final one that needs to be incremented 
  double dt2 = a_dt*a_dt; 
  //Init
  m_delta.init(a_state);
  m_k.init(a_state);
  // Compute k1
  m_f(m_k, a_time, 0.0, a_state); 
  //Increment delta by k1
  m_delta.incrementPositionandVelocity(sixth*a_dt*a_dt, sixth*a_dt, m_k);
  //Update m_k to new position to compute k2
  m_k*=eigth;
  m_k*= dt2;
  // Compute k2
  m_f(m_k, a_time, half*a_dt, a_state);
  //Increment delta by 
  m_delta.incrementPositionandVelocity(  third*a_dt*a_dt, 2.0*third*a_dt, m_k);
  m_k*=half;
  m_k*= dt2;
  // Compute k3
  m_f(m_k, a_time, a_dt, a_state);  
  m_delta.incrementVelocity(sixth*a_dt, m_k);
  // Final increment of both m_x and m_v
  a_state.incrementPositionandVelocity(m_delta);
}
  
#endif 
  
