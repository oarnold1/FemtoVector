#include "FemtoVector.h"

using namespace std;

FemtoVector::FemtoVector():
  ftDx(0.),
  ftDy(0.),
  ftDz(0.),
  ftDT(0.),
  ftDx_gaus(0.),
  ftDy_gaus(0.),
  ftDz_gaus(0.),
  ftDT_gaus(0.),
  frelK(0.),
  fQinv(0.),
  ftTstar(0.),
  ftRout(0.),
  ftRoutstar(0.),
  ftRside(0.),
  ftRlong(0.),
  ftRinv(0.),
  ftPhi(0.),
  ftTheta(0.),
  ftThetakr(0.),
  ftQout(0.),
  ftQoutstar(0.),
  ftQside(0.),
  ftQlong(0.),
  ftRout_gaus(0.),
  ftRoutstar_gaus(0.),
  ftRside_gaus(0.),
  ftRlong_gaus(0.),
  ftRinv_gaus(0.),
  ftPhi_gaus(0.),
  ftTheta_gaus(0.),
  ftRout_cauchy(0.),
  ftRoutstar_cauchy(0.),
  ftRside_cauchy(0.),
  ftRlong_cauchy(0.),
  ftRinv_cauchy(0.),
  ftPhi_cauchy(0.),
  ftTheta_cauchy(0.),
  ftRlab(0.),
  ftEqualTimeApprox(0.),
  ftsigma(0.),
  fbeta(0.),
  fbeta_x(0.),
  fbeta_y(0.),
  fbeta_z(0.),
  finit(kFALSE)
{
  rannum = new TRandom3();
  rannum->SetSeed(0);
}

FemtoVector::~FemtoVector()
{
  rannum->Delete();
  rannum = 0;
}

void FemtoVector::SetPairEmissionPoints(TLorentzVector track_xyzt_A,TLorentzVector track_xyzt_B)
{
  //Takes the spatial Lorentzvectors and calculates the spatio-temporal differences between the particle vectors
  
  //Since lorentz transformations are linear it is allowed to perform first the difference and then the transformation, it commutes
  ftDx = track_xyzt_A.X() - track_xyzt_B.X();
  ftDy = track_xyzt_A.Y() - track_xyzt_B.Y();
  ftDz = track_xyzt_A.Z() - track_xyzt_B.Z();
  ftDT = track_xyzt_A.T() - track_xyzt_B.T();
  
  spatialDiffAB.SetXYZT(ftDx,ftDy,ftDz,ftDT);

  ftRlab = TMath::Sqrt(pow(ftDx,2.) + pow(ftDy,2.) + pow(ftDz,2.));

  finit = kTRUE;

  //it is often interesting how a gaussian source behaves:
  //single particle emitters are fine with sqrt(2.), for two-particle source this gives needed factor when folding two Gaussians:
  if(ftsigma == 0.) std::cout << "Warning, the Gaussian source radius is zero" << std::endl;
  
  double x1 = rannum->Gaus(0,ftsigma);
  double y1 = rannum->Gaus(0,ftsigma);
  double z1 = rannum->Gaus(0,ftsigma);
  double t1 = 0.;//this corresponds to the deltafunction delta(t)

  double x2 = rannum->Gaus(0,ftsigma);
  double y2 = rannum->Gaus(0,ftsigma);
  double z2 = rannum->Gaus(0,ftsigma);
  double t2 = 0.;//this corresponds to the deltafunction delta(t)
  
  ftDx_gaus = x1 - x2;
  ftDy_gaus = y1 - y2;
  ftDz_gaus = z1 - z2;
  ftDT_gaus = t1 - t2;
  
  spatialDiffAB_gaus.SetXYZT(ftDx_gaus,ftDy_gaus,ftDz_gaus,ftDT_gaus);

  //repeat the same story for Cauchy distributions:

  double x1C = rannum->BreitWigner(0,0.5);
  double y1C = rannum->BreitWigner(0,0.5);
  double z1C = rannum->BreitWigner(0,0.5);

  double x2C = rannum->BreitWigner(0,0.5);
  double y2C = rannum->BreitWigner(0,0.5);
  double z2C = rannum->BreitWigner(0,0.5);
  
  ftDx_cauchy = x1C - x2C;
  ftDy_cauchy = y1C - y2C;
  ftDz_cauchy = z1C - z2C;
  ftDT_cauchy = 0.;
  
  spatialDiffAB_cauchy.SetXYZT(ftDx_cauchy,ftDy_cauchy,ftDz_cauchy,ftDT_cauchy);
}

Double_t FemtoVector::relKcalc(TLorentzVector track1,TLorentzVector track2)
{
  //This function is valid for identical AND non-identical pairs, which is a big advantage to other implementations
  //not sensitive to sign of emission, can be implemented, not needed

  TLorentzVector trackSum, track1_cms, track2_cms;
  trackSum = track1 + track2;

  Double_t beta = trackSum.Beta();
  Double_t beta_x = beta*cos(trackSum.Phi())*sin(trackSum.Theta());
  Double_t beta_y = beta*sin(trackSum.Phi())*sin(trackSum.Theta());
  Double_t beta_z = beta*cos(trackSum.Theta());

  track1_cms = track1;
  track2_cms = track2;
  
  track1_cms.Boost(-beta_x,-beta_y,-beta_z);
  track2_cms.Boost(-beta_x,-beta_y,-beta_z);

  TLorentzVector trackCheck,trackCheck2;
  
  trackCheck = track1_cms + track2_cms;
  trackCheck2 = track1 + track2;
  
  trackCheck2.Boost(-beta_x,-beta_y,-beta_z);

  
  TLorentzVector track_relK;
  
  track_relK = track1_cms - track2_cms;
  Double_t relK = 0.5*track_relK.P();
  
  return relK;
}

Double_t FemtoVector::QinvcalcID(TLorentzVector track1,TLorentzVector track2)
{
  TLorentzVector diff12;
  diff12 = track1 - track2;
  

  double product = -diff12.Dot(diff12);
  

  return TMath::Sqrt(product);
}

Double_t FemtoVector::QinvcalcNonID(TLorentzVector trackV0,TLorentzVector trackProton)
{
  TLorentzVector trackQ, trackP;
  trackQ = trackV0 - trackProton;
  trackP = trackV0 + trackProton;
  
  Double_t qinvL = trackQ.Mag2();
  Double_t qP = trackQ.Dot(trackP);
  Double_t pinv = trackP.Mag2();
  
  Double_t QinvLP = TMath::Sqrt(qP*qP/pinv - qinvL);
  
  return QinvLP;
}

Int_t FemtoVector::deltafunction(double value,double width)
{
  //discrete delta function

  if(fabs(value)<width/2.) return 1;
  else return 0;
}

void FemtoVector::SetPairMomentumPoints(TLorentzVector track_pxpypzE_A,TLorentzVector track_pxpypzE_B)
{
  if(!finit) std::cout << "you have to initialize the emission points" << std::endl;
  //Takes the kinematic Lorentzvectors and calculates the kinematic sum of the particle vectors

  sumAB = track_pxpypzE_A + track_pxpypzE_B;
  TLorentzVector diffAB = track_pxpypzE_A - track_pxpypzE_B;
  frelK = relKcalc(track_pxpypzE_A,track_pxpypzE_B);
  
  if(track_pxpypzE_A.M() == track_pxpypzE_B.M()) fQinv = QinvcalcID(track_pxpypzE_A,track_pxpypzE_B);//gives identical results like method above for identical particles
  else fQinv = QinvcalcNonID(track_pxpypzE_A,track_pxpypzE_B);

  //Before the PairBooster can be called, sumAB must be initialized

  //difference of momenta:
  TVector3 tv3kosl = PairBooster(diffAB);


  //calculate the coordinate differences in the PRF from model output:
  TVector3 tv3rosl = PairBooster(spatialDiffAB);
  ftRout = 0.;//not needed at the moment
  ftRoutstar = tv3rosl.X();
  ftRside = tv3rosl.Y();
  ftRlong = tv3rosl.Z();
  ftPhi = tv3rosl.Phi();
  ftTheta = tv3rosl.Theta();
  
  ftThetakr = tv3kosl.Angle(tv3rosl);

  //The one-dimensional radius measured with femto is then just the sum of these three coodrinates
  ftRinv = TMath::Sqrt(pow(ftRoutstar,2.) + pow(ftRside,2.) + pow(ftRlong,2.));
  
  //Test also possible time dependencies:
  TLorentzVector tlorentzosl = PairBooster(track_pxpypzE_A,track_pxpypzE_B,spatialDiffAB);
  
  ftTstar = tlorentzosl.T();
  
  double avMass = 0.5*(track_pxpypzE_A.M() + track_pxpypzE_B.M());

  ftEqualTimeApprox = fabs(tlorentzosl.T()) / (avMass*pow(ftRinv,2.)) * 0.197327;//|t*|/m*r*² (*hbarc)

  //calculate the coordinate differences in the PRF for Gaussian coordinates:
  TVector3 tv3rosl_gaus = PairBooster(spatialDiffAB_gaus);
  ftRout_gaus = 0.;//not needed at the moment
  ftRoutstar_gaus = tv3rosl_gaus.X();
  ftRside_gaus = tv3rosl_gaus.Y();
  ftRlong_gaus = tv3rosl_gaus.Z();
  ftRinv_gaus = TMath::Sqrt(pow(ftRoutstar_gaus,2.) + pow(ftRside_gaus,2.) + pow(ftRlong_gaus,2.));
  ftPhi_gaus = tv3rosl_gaus.Phi();
  ftTheta_gaus = tv3rosl_gaus.Theta();
  
  //calculate the coordinate differences in the PRF for Gaussian coordinates:
  TVector3 tv3rosl_cauchy = PairBooster(spatialDiffAB_cauchy);
  ftRout_cauchy = 0.;//not needed at the moment
  ftRoutstar_cauchy = tv3rosl_cauchy.X();
  ftRside_cauchy = tv3rosl_cauchy.Y();
  ftRlong_cauchy = tv3rosl_cauchy.Z();
  ftRinv_cauchy = TMath::Sqrt(pow(ftRoutstar_cauchy,2.) + pow(ftRside_cauchy,2.) + pow(ftRlong_cauchy,2.));
  ftPhi_cauchy = tv3rosl_cauchy.Phi();
  ftTheta_cauchy = tv3rosl_cauchy.Theta();
  
  fbeta = sumAB.Beta();
  fbeta_x = fbeta*cos(sumAB.Phi())*sin(sumAB.Theta());
  fbeta_y = fbeta*sin(sumAB.Phi())*sin(sumAB.Theta());
  fbeta_z = fbeta*cos(sumAB.Theta());
}

TVector3 FemtoVector::PairBooster(TLorentzVector vectorToBoost)
{
  //If this Method is called, all relevant quantities are also calculated:
  
  
  //LCMS (From M.Lisa paper p.9 - Femtoscopy in heavy ion collisions, arxiv:0505014v2):
  //P=p1+P2 defined on p.4 below Eq. (1)
  //Any 4-vector V can be transformed to the LCMS (P_z=0) via
  //Vout=(PxVx + PyVy)/P_T
  //Vside=(PxVy - PyVx)/P_T
  //Vlong=(P0Vz - PzV0)/M_T
  //Mt^2=P0^2 - P_z^2
  //Pt^2=Px^2 + P_y^2
  //Minv^2=P^2
  
  double totPx = sumAB.Px();
  double totPy = sumAB.Py();
  double totPz = sumAB.Pz();
  

  //We want to transform the coordinate differences, thus V = Difference of coordinates
  double tDVx = vectorToBoost.X();
  double tDVy = vectorToBoost.Y();
  double tDVz = vectorToBoost.Z();
  double tDVT = vectorToBoost.T();

  double totPt = TMath::Sqrt(pow(totPx,2.) + pow(totPy,2.));//=Pt
  double totE = sumAB.E();
  double totMt = TMath::Sqrt(totE*totE - totPz*totPz);//=P0²-Pz²
  double totMinv = TMath::Sqrt(totE*totE - (pow(totPx,2.) + pow(totPy,2.) + pow(totPz,2.)));//=P²
  double totPdotR = totE*tDVT - (totPx*tDVx + totPy*tDVy + totPz*tDVz);//=V \cdot P
  
  //Apply the formulas for the transformation:
  double tVout = (tDVx*totPx + tDVy*totPy)/totPt;//perpendicular to beam ~Pt
  double tVside = (tDVy*totPx - tDVx*totPy)/totPt;//perpendicular to beam ~Pt
  double tVlong = (tDVz*totE - tDVT*totPz)/totMt;//along beam ~Mt ~ Pz
  
  //to transform into the pair rest frame (PRF) one has to boost along total k_T of the pair which affects only the out component
  //Vout'=Minv/Mt (PxVx + PyVy)/Pt - Pt/(Mt*Minv) P \cdot V
  
  double tVoutstar = totMinv/totMt * tVout - totPt/(totMt*totMinv)*totPdotR;//transforms ftRout (LCMS) -> ftRoutstar (PRF)
  
  TVector3 returnVec;
  returnVec.SetXYZ(tVoutstar,tVside,tVlong);

  return returnVec;
}

TLorentzVector FemtoVector::PairBooster(TLorentzVector VecA, TLorentzVector VecB, TLorentzVector vectorToBoost)
{
  TLorentzVector cmsVec;
  cmsVec = VecA + VecB;
  
  Double_t beta = cmsVec.Beta();
  Double_t beta_x = beta*cos(cmsVec.Phi())*sin(cmsVec.Theta());
  Double_t beta_y = beta*sin(cmsVec.Phi())*sin(cmsVec.Theta());
  Double_t beta_z = beta*cos(cmsVec.Theta());
  

  TLorentzVector vectorToBoostPRF;
  vectorToBoostPRF = vectorToBoost;
  
  vectorToBoostPRF.Boost(-beta_x,-beta_y,-beta_z);
  
  return vectorToBoostPRF;
}
