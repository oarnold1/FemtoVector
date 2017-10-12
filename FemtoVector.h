//FemtoVector.h
//author: O.Arnold

//code can be used (together with ROOT) to analyse output of simulation models (in this case EPOS)

#ifndef FEMTOVECTOR
#define FEMTOVECTOR

#include "TLorentzVector.h"
#include <iostream>
#include "TRandom3.h"
#include "TLorentzVector.h"
#include "TVector3.h"

class FemtoVector
{
 private:
  Double_t ftDx,ftDy,ftDz,ftDT;
  Double_t ftDPx,ftDPy,ftDPz,ftDE;
  Double_t ftDx_gaus,ftDy_gaus,ftDz_gaus,ftDT_gaus;
  Double_t ftDx_cauchy,ftDy_cauchy,ftDz_cauchy,ftDT_cauchy;
  Double_t frelK,fQinv;
  Double_t ftTstar,ftRout,ftRoutstar,ftRside,ftRlong,ftRinv,ftPhi,ftTheta,ftThetakr;
  Double_t ftQout,ftQoutstar,ftQside,ftQlong;
  Double_t ftRout_gaus,ftRoutstar_gaus,ftRside_gaus,ftRlong_gaus,ftRinv_gaus,ftPhi_gaus,ftTheta_gaus;
  Double_t ftRout_cauchy,ftRoutstar_cauchy,ftRside_cauchy,ftRlong_cauchy,ftRinv_cauchy,ftPhi_cauchy,ftTheta_cauchy;
  Double_t ftRlab;
  Double_t ftEqualTimeApprox;
  Double_t ftsigma;
  Double_t fbeta,fbeta_x,fbeta_y,fbeta_z;
  Bool_t finit;
  TRandom3 *rannum;
  TLorentzVector sumAB,spatialDiffAB,spatialDiffAB_gaus,spatialDiffAB_cauchy;
  
  Double_t relKcalc(TLorentzVector track1,TLorentzVector track2);
  Double_t QinvcalcID(TLorentzVector track1,TLorentzVector track2);
  Double_t QinvcalcNonID(TLorentzVector trackV0,TLorentzVector trackProton);

 public:
  
  FemtoVector();
  virtual ~FemtoVector();
  
  //Functions
  Int_t deltafunction(double value,double width);
  TVector3 PairBooster(TLorentzVector vectorToBoost);
  TLorentzVector PairBooster(TLorentzVector VecA, TLorentzVector VecB, TLorentzVector vectorToBoost);

  //Setter
  void SetPairEmissionPoints(TLorentzVector track_xyzt_A,TLorentzVector track_xyzt_B);
  void SetPairMomentumPoints(TLorentzVector track_pxpypzE_A,TLorentzVector track_pxpypzE_B);
  void Setsourceval(Double_t sigma) {ftsigma = sigma;};
  
  //Getter
  Double_t GetrelK() const {return frelK;};
  Double_t GetQinv() const {return fQinv;};

  Double_t GetTstar() const {return ftTstar;};
  Double_t GetRout() const {return ftRout;};
  Double_t GetRoutstar() const {return ftRoutstar;};
  Double_t GetRside() const {return ftRside;};
  Double_t GetRlong() const {return ftRlong;};
  Double_t GetRinv() const {return ftRinv;};
  Double_t GetEqualTimeApproximation() const {return ftEqualTimeApprox;};
  Double_t GetPhi() const {return ftPhi;};
  Double_t GetTheta() const {return ftTheta;};
  Double_t GetThetakr() const {return ftThetakr;};

  Double_t GetRout_gaus() const {return ftRout_gaus;};
  Double_t GetRoutstar_gaus() const {return ftRoutstar_gaus;};
  Double_t GetRside_gaus() const {return ftRside_gaus;};
  Double_t GetRlong_gaus() const {return ftRlong_gaus;};
  Double_t GetRinv_gaus() const {return ftRinv_gaus;};
  Double_t GetPhi_gaus() const {return ftPhi_gaus;};
  Double_t GetTheta_gaus() const {return ftTheta_gaus;};

  Double_t GetRoutstar_cauchy() const {return ftRoutstar_cauchy;};
  Double_t GetRside_cauchy() const {return ftRside_cauchy;};
  Double_t GetRlong_cauchy() const {return ftRlong_cauchy;};
  Double_t GetRinv_cauchy() const {return ftRinv_cauchy;};
  Double_t GetPhi_cauchy() const {return ftPhi_cauchy;};
  Double_t GetTheta_cauchy() const {return ftTheta_cauchy;};

  Double_t GetBeta() const {return fbeta;};
  Double_t GetBetax() const {return fbeta_x;};
  Double_t GetBetay() const {return fbeta_y;};
  Double_t GetBetaz() const {return fbeta_z;};

  Double_t GetRlab() const {return ftRlab;};
  
 protected:
  
#ifdef __ROOT__
  ClassDef(FemtoVector,1);
#endif    
};
#endif
