#include "pollackProjectHeaders.h"

#include "pollackLinkDef.h"

#include "pollackProjectDict.cxx"

struct DeleteObjectFunctor {
   template <typename T>
   void operator()(const T *ptr) const {
      delete ptr;
   }
   template <typename T, typename Q>
   void operator()(const std::pair<T,Q> &) const {
      // Do nothing
   }
   template <typename T, typename Q>
   void operator()(const std::pair<T,Q*> &ptr) const {
      delete ptr.second;
   }
   template <typename T, typename Q>
   void operator()(const std::pair<T*,Q> &ptr) const {
      delete ptr.first;
   }
   template <typename T, typename Q>
   void operator()(const std::pair<T*,Q*> &ptr) const {
      delete ptr.first;
      delete ptr.second;
   }
};

#ifndef TCMET_cxx
#define TCMET_cxx
TCMET::TCMET() {
}
TCMET::TCMET(const TCMET & rhs)
   : TVector2(const_cast<TCMET &>( rhs ))
   , _genMET(const_cast<TCMET &>( rhs )._genMET)
   , _sumEt(const_cast<TCMET &>( rhs )._sumEt)
   , _muonFraction(const_cast<TCMET &>( rhs )._muonFraction)
   , _neutralHadronFraction(const_cast<TCMET &>( rhs )._neutralHadronFraction)
   , _neutralEMFraction(const_cast<TCMET &>( rhs )._neutralEMFraction)
   , _chargedHadronFraction(const_cast<TCMET &>( rhs )._chargedHadronFraction)
   , _chargedEMFraction(const_cast<TCMET &>( rhs )._chargedEMFraction)
   , _unCorPhi(const_cast<TCMET &>( rhs )._unCorPhi)
   , _Significance(const_cast<TCMET &>( rhs )._Significance)
   , _SigmaX2(const_cast<TCMET &>( rhs )._SigmaX2)
{
   // This is NOT a copy constructor. This is actually a move constructor (for stl container's sake).
   // Use at your own risk!
   if (&rhs) {} // avoid warning about unused parameter
}
TCMET::~TCMET() {
}
#endif // TCMET_cxx

#ifndef TCJet_cxx
#define TCJet_cxx
TCJet::TCJet() {
}
TCJet::TCJet(const TCJet & rhs)
   : TCPhysObject(const_cast<TCJet &>( rhs ))
   , _vtxSumPtFrac(const_cast<TCJet &>( rhs )._vtxSumPtFrac)
   , _vtxSumPt(const_cast<TCJet &>( rhs )._vtxSumPt)
   , _vtxTrackFrac(const_cast<TCJet &>( rhs )._vtxTrackFrac)
   , _vtxNTracks(const_cast<TCJet &>( rhs )._vtxNTracks)
   , _vtxSumPtIndex(const_cast<TCJet &>( rhs )._vtxSumPtIndex)
   , _vtxCountIndex(const_cast<TCJet &>( rhs )._vtxCountIndex)
   , _jesUncertainty(const_cast<TCJet &>( rhs )._jesUncertainty)
   , _chHadFrac(const_cast<TCJet &>( rhs )._chHadFrac)
   , _neuHadFrac(const_cast<TCJet &>( rhs )._neuHadFrac)
   , _chEmFrac(const_cast<TCJet &>( rhs )._chEmFrac)
   , _neuEmFrac(const_cast<TCJet &>( rhs )._neuEmFrac)
   , _numConstit(const_cast<TCJet &>( rhs )._numConstit)
   , _numChPart(const_cast<TCJet &>( rhs )._numChPart)
   , _bDiscrMap(const_cast<TCJet &>( rhs )._bDiscrMap)
   , _jetFlavor(const_cast<TCJet &>( rhs )._jetFlavor)
{
   // This is NOT a copy constructor. This is actually a move constructor (for stl container's sake).
   // Use at your own risk!
   if (&rhs) {} // avoid warning about unused parameter
   TCJet &modrhs = const_cast<TCJet &>( rhs );
   modrhs._bDiscrMap.clear();
}
TCJet::~TCJet() {
}
#endif // TCJet_cxx

#ifndef TCPhysObject_cxx
#define TCPhysObject_cxx
TCPhysObject::TCPhysObject() {
}
TCPhysObject::TCPhysObject(const TCPhysObject & rhs)
   : TLorentzVector(const_cast<TCPhysObject &>( rhs ))
   , _vtx(const_cast<TCPhysObject &>( rhs )._vtx)
   , _IdMap(const_cast<TCPhysObject &>( rhs )._IdMap)
   , _IsoMap(const_cast<TCPhysObject &>( rhs )._IsoMap)
   , _charge(const_cast<TCPhysObject &>( rhs )._charge)
   , _type(const_cast<TCPhysObject &>( rhs )._type)
   , _isPF(const_cast<TCPhysObject &>( rhs )._isPF)
   , _isReco(const_cast<TCPhysObject &>( rhs )._isReco)
{
   // This is NOT a copy constructor. This is actually a move constructor (for stl container's sake).
   // Use at your own risk!
   if (&rhs) {} // avoid warning about unused parameter
   TCPhysObject &modrhs = const_cast<TCPhysObject &>( rhs );
   modrhs._IdMap.clear();
   modrhs._IsoMap.clear();
   modrhs._type.clear();
}
TCPhysObject::~TCPhysObject() {
}
#endif // TCPhysObject_cxx

#ifndef TCElectron_cxx
#define TCElectron_cxx
TCElectron::TCElectron() {
}
TCElectron::TCElectron(const TCElectron & rhs)
   : TCPhysObject(const_cast<TCElectron &>( rhs ))
   , _ptError(const_cast<TCElectron &>( rhs )._ptError)
   , _hadOverEm(const_cast<TCElectron &>( rhs )._hadOverEm)
   , _dPhiSuperCluster(const_cast<TCElectron &>( rhs )._dPhiSuperCluster)
   , _dEtaSuperCluster(const_cast<TCElectron &>( rhs )._dEtaSuperCluster)
   , _sigmaIetaIeta(const_cast<TCElectron &>( rhs )._sigmaIetaIeta)
   , _eOverP(const_cast<TCElectron &>( rhs )._eOverP)
   , _fBrem(const_cast<TCElectron &>( rhs )._fBrem)
   , _r9(const_cast<TCElectron &>( rhs )._r9)
   , _scEta(const_cast<TCElectron &>( rhs )._scEta)
   , _convVeto(const_cast<TCElectron &>( rhs )._convVeto)
   , _convMissHits(const_cast<TCElectron &>( rhs )._convMissHits)
   , _isEB(const_cast<TCElectron &>( rhs )._isEB)
   , _isEE(const_cast<TCElectron &>( rhs )._isEE)
   , _isInGap(const_cast<TCElectron &>( rhs )._isInGap)
   , _normalizedChi2(const_cast<TCElectron &>( rhs )._normalizedChi2)
   , _numberOfValidPixelHits(const_cast<TCElectron &>( rhs )._numberOfValidPixelHits)
   , _numberOfValidTrackerHits(const_cast<TCElectron &>( rhs )._numberOfValidTrackerHits)
   , _numberOfLostPixelHits(const_cast<TCElectron &>( rhs )._numberOfLostPixelHits)
   , _numberOfLostTrackerHits(const_cast<TCElectron &>( rhs )._numberOfLostTrackerHits)
   , _cut95(const_cast<TCElectron &>( rhs )._cut95)
   , _cut90(const_cast<TCElectron &>( rhs )._cut90)
   , _cut85(const_cast<TCElectron &>( rhs )._cut85)
   , _cut80(const_cast<TCElectron &>( rhs )._cut80)
   , _cut70(const_cast<TCElectron &>( rhs )._cut70)
   , _cut60(const_cast<TCElectron &>( rhs )._cut60)
   , _regressionMomCombP4(const_cast<TCElectron &>( rhs )._regressionMomCombP4)
{
   // This is NOT a copy constructor. This is actually a move constructor (for stl container's sake).
   // Use at your own risk!
   if (&rhs) {} // avoid warning about unused parameter
}
TCElectron::~TCElectron() {
}
#endif // TCElectron_cxx

#ifndef TCMuon_cxx
#define TCMuon_cxx
TCMuon::TCMuon() {
}
TCMuon::TCMuon(const TCMuon & rhs)
   : TCPhysObject(const_cast<TCMuon &>( rhs ))
   , _ptError(const_cast<TCMuon &>( rhs )._ptError)
   , _isPF(const_cast<TCMuon &>( rhs )._isPF)
   , _isTRK(const_cast<TCMuon &>( rhs )._isTRK)
   , _isGLB(const_cast<TCMuon &>( rhs )._isGLB)
   , _caloComp(const_cast<TCMuon &>( rhs )._caloComp)
   , _segComp(const_cast<TCMuon &>( rhs )._segComp)
   , _numberOfMatches(const_cast<TCMuon &>( rhs )._numberOfMatches)
   , _numberOfMatchedStations(const_cast<TCMuon &>( rhs )._numberOfMatchedStations)
   , _numberOfValidPixelHits(const_cast<TCMuon &>( rhs )._numberOfValidPixelHits)
   , _numberOfValidTrackerHits(const_cast<TCMuon &>( rhs )._numberOfValidTrackerHits)
   , _numberOfLostPixelHits(const_cast<TCMuon &>( rhs )._numberOfLostPixelHits)
   , _numberOfLostTrackerHits(const_cast<TCMuon &>( rhs )._numberOfLostTrackerHits)
   , _numberOfValidMuonHits(const_cast<TCMuon &>( rhs )._numberOfValidMuonHits)
   , _trackLayersWithMeasurement(const_cast<TCMuon &>( rhs )._trackLayersWithMeasurement)
   , _normalizedChi2(const_cast<TCMuon &>( rhs )._normalizedChi2)
{
   // This is NOT a copy constructor. This is actually a move constructor (for stl container's sake).
   // Use at your own risk!
   if (&rhs) {} // avoid warning about unused parameter
}
TCMuon::~TCMuon() {
}
#endif // TCMuon_cxx

#ifndef TCPhoton__CrystalInfo_cxx
#define TCPhoton__CrystalInfo_cxx
TCPhoton::CrystalInfo::CrystalInfo() {
}
TCPhoton::CrystalInfo::CrystalInfo(const CrystalInfo & rhs)
   : rawId(const_cast<CrystalInfo &>( rhs ).rawId)
   , ieta(const_cast<CrystalInfo &>( rhs ).ieta)
   , iphi(const_cast<CrystalInfo &>( rhs ).iphi)
   , ix(const_cast<CrystalInfo &>( rhs ).ix)
   , iy(const_cast<CrystalInfo &>( rhs ).iy)
   , energy(const_cast<CrystalInfo &>( rhs ).energy)
   , time(const_cast<CrystalInfo &>( rhs ).time)
   , timeErr(const_cast<CrystalInfo &>( rhs ).timeErr)
   , recoFlag(const_cast<CrystalInfo &>( rhs ).recoFlag)
{
   // This is NOT a copy constructor. This is actually a move constructor (for stl container's sake).
   // Use at your own risk!
   if (&rhs) {} // avoid warning about unused parameter
}
TCPhoton::CrystalInfo::~CrystalInfo() {
}
#endif // TCPhoton__CrystalInfo_cxx

#ifndef TCPhoton_cxx
#define TCPhoton_cxx
TCPhoton::TCPhoton() {
   _crysArray = 0;
}
TCPhoton::TCPhoton(const TCPhoton & rhs)
   : TCPhysObject(const_cast<TCPhoton &>( rhs ))
   , _normChi2(const_cast<TCPhoton &>( rhs )._normChi2)
   , _hadOverEm(const_cast<TCPhoton &>( rhs )._hadOverEm)
   , _sigmaIEtaIEta(const_cast<TCPhoton &>( rhs )._sigmaIEtaIEta)
   , _r9(const_cast<TCPhoton &>( rhs )._r9)
   , _sigmaIPhiIPhi(const_cast<TCPhoton &>( rhs )._sigmaIPhiIPhi)
   , _e2OverE9(const_cast<TCPhoton &>( rhs )._e2OverE9)
   , _trackVeto(const_cast<TCPhoton &>( rhs )._trackVeto)
   , _SCdPhi(const_cast<TCPhoton &>( rhs )._SCdPhi)
   , _SCdEta(const_cast<TCPhoton &>( rhs )._SCdEta)
   , _SCeta(const_cast<TCPhoton &>( rhs )._SCeta)
   , _SCphi(const_cast<TCPhoton &>( rhs )._SCphi)
   , _SCenergy(const_cast<TCPhoton &>( rhs )._SCenergy)
   , _convVeto(const_cast<TCPhoton &>( rhs )._convVeto)
   , _crysArray(const_cast<TCPhoton &>( rhs )._crysArray)
   , _nCrystals(const_cast<TCPhoton &>( rhs )._nCrystals)
{
   // This is NOT a copy constructor. This is actually a move constructor (for stl container's sake).
   // Use at your own risk!
   if (&rhs) {} // avoid warning about unused parameter
   TCPhoton &modrhs = const_cast<TCPhoton &>( rhs );
   modrhs._crysArray = 0;
}
TCPhoton::~TCPhoton() {
   delete _crysArray;   _crysArray = 0;
}
#endif // TCPhoton_cxx

#ifndef TCTriggerObject_cxx
#define TCTriggerObject_cxx
TCTriggerObject::TCTriggerObject() {
}
TCTriggerObject::TCTriggerObject(const TCTriggerObject & rhs)
   : TLorentzVector(const_cast<TCTriggerObject &>( rhs ))
   , _id(const_cast<TCTriggerObject &>( rhs )._id)
   , _HLTName(const_cast<TCTriggerObject &>( rhs )._HLTName)
   , _moduleName(const_cast<TCTriggerObject &>( rhs )._moduleName)
{
   // This is NOT a copy constructor. This is actually a move constructor (for stl container's sake).
   // Use at your own risk!
   if (&rhs) {} // avoid warning about unused parameter
   TCTriggerObject &modrhs = const_cast<TCTriggerObject &>( rhs );
   modrhs._HLTName.clear();
   modrhs._moduleName.clear();
}
TCTriggerObject::~TCTriggerObject() {
}
#endif // TCTriggerObject_cxx

#ifndef TCGenJet_cxx
#define TCGenJet_cxx
TCGenJet::TCGenJet() {
}
TCGenJet::TCGenJet(const TCGenJet & rhs)
   : TCPhysObject(const_cast<TCGenJet &>( rhs ))
   , _progenitorP4(const_cast<TCGenJet &>( rhs )._progenitorP4)
   , _hadEnergy(const_cast<TCGenJet &>( rhs )._hadEnergy)
   , _emEnergy(const_cast<TCGenJet &>( rhs )._emEnergy)
   , _invEnergy(const_cast<TCGenJet &>( rhs )._invEnergy)
   , _auxEnergy(const_cast<TCGenJet &>( rhs )._auxEnergy)
   , _numConstit(const_cast<TCGenJet &>( rhs )._numConstit)
   , _numChPart(const_cast<TCGenJet &>( rhs )._numChPart)
   , _jetFlavor(const_cast<TCGenJet &>( rhs )._jetFlavor)
{
   // This is NOT a copy constructor. This is actually a move constructor (for stl container's sake).
   // Use at your own risk!
   if (&rhs) {} // avoid warning about unused parameter
}
TCGenJet::~TCGenJet() {
}
#endif // TCGenJet_cxx

#ifndef TCGenParticle_cxx
#define TCGenParticle_cxx
TCGenParticle::TCGenParticle() {
}
TCGenParticle::TCGenParticle(const TCGenParticle & rhs)
   : TCPhysObject(const_cast<TCGenParticle &>( rhs ))
   , mother(const_cast<TCGenParticle &>( rhs ).mother)
   , grandmother(const_cast<TCGenParticle &>( rhs ).grandmother)
   , PDGID(const_cast<TCGenParticle &>( rhs ).PDGID)
   , status(const_cast<TCGenParticle &>( rhs ).status)
   , isParton_(const_cast<TCGenParticle &>( rhs ).isParton_)
{
   // This is NOT a copy constructor. This is actually a move constructor (for stl container's sake).
   // Use at your own risk!
   if (&rhs) {} // avoid warning about unused parameter
}
TCGenParticle::~TCGenParticle() {
}
#endif // TCGenParticle_cxx

#ifndef TCPrimaryVtx_cxx
#define TCPrimaryVtx_cxx
TCPrimaryVtx::TCPrimaryVtx() {
}
TCPrimaryVtx::TCPrimaryVtx(const TCPrimaryVtx & rhs)
   : TVector3(const_cast<TCPrimaryVtx &>( rhs ))
   , _nDof(const_cast<TCPrimaryVtx &>( rhs )._nDof)
   , _chi2(const_cast<TCPrimaryVtx &>( rhs )._chi2)
   , _isFake(const_cast<TCPrimaryVtx &>( rhs )._isFake)
   , _nTracks(const_cast<TCPrimaryVtx &>( rhs )._nTracks)
   , _sumPt2Trks(const_cast<TCPrimaryVtx &>( rhs )._sumPt2Trks)
{
   // This is NOT a copy constructor. This is actually a move constructor (for stl container's sake).
   // Use at your own risk!
   if (&rhs) {} // avoid warning about unused parameter
}
TCPrimaryVtx::~TCPrimaryVtx() {
}
#endif // TCPrimaryVtx_cxx

