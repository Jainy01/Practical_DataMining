#ifndef OSC_OSCCALCULATORPMNSL3_H
#define OSC_OSCCALCULATORPMNSL3_H

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// \file   OscCalculatorPMNS.h                                          //
//                                                                      //
/// \brief Implementation of oscillations of neutrinos in matter in a
///        three-neutrino framework. 
///
/// Three 3-flavor oscillations based on the following reference:
///.....................................................................
///
/// PHYSICAL REVIEW D       VOLUME 22, NUMBER 11         1 DECEMBER 1980
///
///             Matter effects on three-neutrino oscillation
///
///                      V. Barger and K. Whisnant
///    Physics Department, U. of Wisconsin, Madison, Wisconsin 53706
///
///                            S. Pakvasa
///   Physics Department, U. of Hawaii at Manoa, Honolulu, Hawaii 96822
///
///                         R.J.N. Phillips
///        Rutherford Laboratory, Chilton, Didcot, Oxon, England
///                    (Received 4 August 1980)
///
///                            22 2718
///                            --
///.....................................................................
///
/// \version $Id: PMNS.h,v 1.2 2012/08/29 01:33:49 bckhouse Exp $
///
/// @author messier@indiana.edu
////////////////////////////////////////////////////////////////////////

#include <list>
#include <complex>

using namespace std;

  // Some useful complex numbers
static std::complex<double> zero(0.0,0.0);
static std::complex<double> one(1.0,0.0);

// Unit conversion constants
static const double kK1     = 2.53386551601e-00; ///< (1/2)*(1000/hbarc)
static const double kK2     = 4.62711492217e-09; ///< mole/GeV^2/cm^3 to eV^2
static const double kGeV2eV = 1.0E9;             ///< GeV to eV

/// Adapt the \ref PMNS calculator to standard interface
class OscCalculatorPMNSL3
{
  
 public:

  OscCalculatorPMNSL3();
  virtual ~OscCalculatorPMNSL3();


  

  double P(int flavBefore, int flavAfter, double E);
  
  void SetL     (double L     ){fPropDirty = true; fL      = L;}
  void SetRho   (double rho   ){fPropDirty = true; fRho    = rho;}
  void SetDmsq21(double dmsq21){fDmDirty   = true; fDmsq21 = dmsq21;}
  void SetDmsq32(double dmsq32){fDmDirty   = true; fDmsq32 = dmsq32;}
  void SetTh12  (double th12  ){fMixDirty  = true; fTh12   = th12;}
  void SetTh13  (double th13  ){fMixDirty  = true; fTh13   = th13;}
  void SetTh23  (double th23  ){fMixDirty  = true; fTh23   = th23;}
  void SetdCP   (double dCP   ){fMixDirty  = true; fdCP    = dCP;}
  
  double GetL     () const { return fL      ; }
  double GetRho   () const { return fRho    ; }
  double GetDmsq21() const { return fDmsq21 ; }
  double GetDmsq32() const { return fDmsq32 ; }
  double GetTh12  () const { return fTh12   ; }
  double GetTh13  () const { return fTh13   ; }
  double GetTh23  () const { return fTh23   ; }
  double GetdCP   () const { return fdCP    ; }
  
  /// Construct the PMNS matrix
  /// @param th12 - the "solar" angle (radians)
  /// @param th23 - the "atmospheric" angle (radians)
  /// @param th13 - the "reactor" angle (radians)
  /// @param deltacp - CPV phase angle (radians)
  /// @param dms12 - m^2_2 - m^2_1 (eV^2)
  /// @param dms23 - m^2_3 - m^2_2 (eV^2)
  //  double th12;
  // double th23;
  //double th13;
  //double deltacp;
  //double dms12;
  //double dms23;
  
  /// Print the oscillation matrix
  void PrintMix() const;
  
  /// Print the matrix of mass-squared differneces
  void PrintDeltaMsqrs() const;
  
  /// Return the oscillation probability from flavor i to flavor j
  /// @param i - initial flavor (0,1,2) = (nue,numu,nutau)
  /// @param j - final   flavor (0,1,2) = (nue,numu,nutau)
  double PFinal(int i, int j) const;
  
  /// Set the parameters of the PMNS matrix
  /// @param th12    - The angle theta_12 in radians
  /// @param th23    - The angle theta_23 in radians
  /// @param th13    - The angle theta_13 in radians
  /// @param deltacp - The CPV phase delta_CP in radians
  void SetMix(double th12, double th23, double th13, double deltacp);
  
  /// Set the mixing matrix using the old convention used in the
  /// original 1980 paper. 
  ///
  /// \warning Useful for testing but should not be used for
  /// calculations. Use SetMix above instead
  ///
  /// @param th1,th2,th3,deltacp - parameters of matrix
  void SetMixBWCP(double th1, double th2, double th3, double deltacp);
    
  /// Set the mass-splittings
  /// @param dm21 - m2^2-m1^2 in eV^2
  /// @param dm32 - m3^2-m2^2 in eV^2
  void SetDeltaMsqrs(double dm21, double dm32);
  
  /// Propagate a neutrino through vacuum
  /// @param L    - flight distance in km
  /// @param E    - neutrino energy in GeV
  /// @param anti - +1 = neutrino case, -1 = anti-neutrino case
  void PropVacuum(double L, double E, int anti);
  
  /// Propagate a neutrino through a slab of matter
  /// @param L - length of slab (km)
  /// @param E - neutrino energy in GeV
  /// @param Ne - electron number density of matter in mole/cm^3
  /// @param anti - +1 = neutrino case, -1 = anti-neutrino case
  void PropMatter(double L, double E, double Ne, int anti);
  void PropMatter(const std::list<double>& L,
		  double                   E,
		  const std::list<double>& Ne,
		  int anti);
  
  /// Erase memory of neutrino propagate and reset transition matrix
  /// to unity. Preserves values of mixing and mass-splittings
  void Reset();

 protected:
  //PMNS fPMNS;
  
  // Set by the user. Generally useful to derived classes
  double fRho; // density (g/cm^3)
  double fL; // baseline (km)
  double fDmsq21;
  double fDmsq32;
  double fTh12;
  double fTh13;
  double fTh23;
  double fdCP;
    
  bool fMixDirty;
  bool fDmDirty;
  bool fPropDirty;
  double fPrevE;
  int fPrevAnti;
    
 private:
  // A shorthand...
  typedef std::complex<double> complex;
  
  /// Multiply two matrices: A = B*C
  /// @param A - output matrix
  /// @param B - input matrix
  /// @param C - input matrix
  void Multi(complex A[][3], const complex B[][3], const complex C[][3]);
  
  /// Find the transition matrix for vacuum oscillations using Eqn.2
  /// @param A     - output transition matrix
  /// @param U     - neutrino mixing matrix
  /// @param Udagg - adjoint of neutrino mixing matrix
  /// @param dmsqr - matrix of mass-splittings m_i^2-m_j^2 in eV^2
  /// @param L     - neutrino flight distance in km
    /// @param E     - neutrino energy in GeV
  void EvalEqn2(complex A[][3],
		const complex U[][3],
		const complex Udagg[][3],
		const double  dmsqr[][3],
		double L,
		double E); 
  
  /// Find the Hamiltonian using Eqn.5
  /// @param twoEH - Hamiltonian matrix scaled by 2E in units of eV
  /// @param U     - neutrino mixing matrix
  /// @param Udagg - adjoint of neutrino mixing matrix
  /// @param dmsqr - Matrix of mass-splittings m_i^2-m_j^2 in eV^2
  /// @param E     - neutrino energy in GeV
  /// @param Gf    - Fermi's constant in units of GeV^-2
  /// @param Ne    - electron number density in mole/cm^3
  void EvalEqn5(complex       twoEH[][3],
		const complex U[][3],
		const complex Udagg[][3],
		const double  dmsqr[][3],
		double        E,
		double        Gf,
		double        Ne);
  
  /// Find the transition matrix for matter using Eqn.10
  /// @param A     - the output transition matrix
  /// @param U     - the neutrino mixing matrix
  /// @param X     - matrix resulting from Eqn.11
  /// @param Udagg - the adjoint of the mixing matrix
  void EvalEqn10(complex       A[][3],
		 const complex U[][3],
		 const complex X[][3],
		 const complex Udagg[][3]);
  
  /// Compute the X matrix from Eqn.11
  /// @param X - the "X" matrix. Sorry, no better name...
  /// @param L - the neutrino flight distance in km
  /// @param E - the neutrino energy in GeV
  /// @param twoEH - the Hamiltonian in units of eV
  /// @param Msqr - the matter eigenvalues in units of eV^2
  /// @param dMsqr - the matrix of eigenvalue differences M_i^2-M_j^2
  void EvalEqn11(complex X[][3],
		 double L, double E, 
		 const complex twoEH[][3],
		 const double  Msqr[],
		 const double  dMsqr[][3]);
    
  /// Compute the matter eigenvalues according to Eqn.21
  /// @param Msqr  - Output eigenvalues in eV^2
  /// @param alpha - See Eqn22
  /// @param beta  - See Eqn22
  /// @param gamma - See Eqn22
  void EvalEqn21(double Msqr[],
		 double alpha,
		 double beta,
		 double gamma);
  
  /// Compute the expressions alpha, beta, gamma which appear in
  /// Eqns.21 and 22
  /// @param alpha - output parameter
  /// @param beta  - output parameter
  /// @param gamma - output parameter
  /// @param E     - neutrino energy in GeV
  /// @param Gf    - Fermi's constant in units of GeV^-2
  /// @param Ne    - electron number denstiy in mole/cm^3
  /// @param dmsqr - matrix of vacuum mass-splittings m_i^2-m_j^2
  /// @param U     - neutrino mixing matrix
  void EvalEqn22(double& alpha,
		 double& beta,
		 double& gamma,
		 double  E,
		 double  Gf,
		 double  Ne,
		 const double dmsqr[][3],
		 const complex U[][3]);
  
  /// Sort the matter eignvalues to that they appear in the same order
  /// as the vacuum eigenvalues.
  /// @param dMsqr   - Output matrix of eigenvalue differences M_i^2-M_j^2
  /// @param dmsqr   - the vacuum eigenvalues
  /// @param MsqrVac - vacuum eigenvalues
  /// @param Msqr    - Eigenvalues in matter. In: unsorted, out: sorted
  void SortEigenvalues(double       dMsqr[][3],
		       const double dmsqr[][3],
		       const double MsqrVac[],
		       double       Msqr[]);
  
  /// Utility to print a genetrix complex matrix M
  void DumpMatrix(const complex M[][3]) const;


  double  fdmsqr[3][3]; ///< m^2_i - m^2_j in vacuum
  complex fU[3][3];     ///< PMNS Matrix in all its forms
  complex fUdagg[3][3];
  complex fUstar[3][3];
  complex fUtran[3][3];
  complex fA[3][3];     ///< The neutrino transition matrix

};


#endif
