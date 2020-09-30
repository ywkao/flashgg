#include "flashgg/DataFormats/interface/THQHadronicTag.h"
#include <algorithm>

using namespace flashgg;

THQHadronicTag::THQHadronicTag() : DiPhotonTagBase::DiPhotonTagBase()
{}

THQHadronicTag::~THQHadronicTag()
{}

// N.B. Other attributes are set using methods in header file
THQHadronicTag::THQHadronicTag( edm::Ptr<DiPhotonCandidate> diPho, edm::Ptr<DiPhotonMVAResult> mvares ) : THQHadronicTag::THQHadronicTag( diPho, *mvares ) {}
THQHadronicTag::THQHadronicTag( edm::Ptr<DiPhotonCandidate> dipho, DiPhotonMVAResult mvares ) : DiPhotonTagBase::DiPhotonTagBase( dipho, mvares ) {}
// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
