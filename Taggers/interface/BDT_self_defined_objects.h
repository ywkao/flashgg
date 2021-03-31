#ifndef _BDTSEFLOBJECTS_H_
#define _BDTSEFLOBJECTS_H_
#include <memory>
#include <vector>
#include <string>
#include <utility>
#include <math.h> //for the use of pow
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "TMVA/Reader.h"
#include <iostream>
#include <algorithm>

namespace flashgg {

typedef math::PtEtaPhiMLorentzVectorD BDT_ptvec;
//https://cmsdoxygen.web.cern.ch/CMSSW_5_3_0/doc/html/d6/de0/LorentzVector_8h_source.html
//https://root.cern.ch/doc/master/classROOT_1_1Math_1_1LorentzVector-members.html

class BDT_Photon : public BDT_ptvec{
 public:
  BDT_Photon() : BDT_ptvec(0,0,0,0){};
  BDT_Photon(float pt,float eta, float phi, float mass, float IDMVA)
      : BDT_ptvec(pt,eta,phi,mass), _IDMVA(IDMVA){};
  float _IDMVA = 0;
  float IDMVA() const {return _IDMVA;}
};

class BDT_Jet : public BDT_ptvec{
 public:
  BDT_Jet() : BDT_ptvec(0,0,0,0){};
  BDT_Jet(float pt,float eta, float phi, float mass, float deepcsv)
      : BDT_ptvec(pt,eta,phi,mass), _deepcsv(deepcsv){};
  float _deepcsv = 0;
  float deepcsv() const {return _deepcsv;}
};

class BDT_Lepton : public BDT_ptvec{
 public:
  BDT_Lepton() : BDT_ptvec(0,0,0,0){};
  BDT_Lepton(float pt,float eta, float phi, float mass, float id)
      : BDT_ptvec(pt,eta,phi,mass), _id(id){};
  float _id = 0;
  float id() const {return _id;}
};

}

#endif // _BDTSEFLOBJECTS_H_
