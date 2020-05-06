#include <cmath>
#include <TLorentzVector.h>
#include <TVector3.h>

std::pair<TVector3, TVector3> csBoostedProtons(TLorentzVector dilepton) {
    float protonMass = 0.938272;
    float energy = 6500;
    int zsign = dilepton.Z() > 0 ? 1 : -1;
    TLorentzVector proton1(0., 0., zsign*energy, hypot(energy, protonMass));
    TLorentzVector proton2(0., 0., -1*zsign*energy, hypot(energy, protonMass));
    proton1.Boost(-1*dilepton.BoostVector());
    proton2.Boost(-1*dilepton.BoostVector());
    return std::make_pair<TVector3, TVector3>(proton1.Vect(), proton2.Vect());
}

const TVector3 csframe(TLorentzVector dilepton) {
    std::pair<TVector3, TVector3> protons = csBoostedProtons(dilepton);
    TVector3 csAxis = (protons.first.Unit() - protons.second.Unit()).Unit();
    return csAxis;
}

const TVector3 csframeY(TLorentzVector dilepton) {
    std::pair<TVector3, TVector3> protons = csBoostedProtons(dilepton);
    TVector3 csYAxis = protons.first.Unit().Cross(protons.second.Unit());
    return csYAxis.Unit();
}

const TVector3 csframeX(TLorentzVector dilepton) {
    TVector3 csAxis = csframe(dilepton);
    TVector3 csYAxis = csframeY(dilepton);
    TVector3 csXAxis = csYAxis.Cross(csAxis);
    return csXAxis.Unit();
}

float cosThetaCS(ROOT::Math::PtEtaPhiMVector lplus, ROOT::Math::PtEtaPhiMVector lminus) {
    ROOT::Math::PtEtaPhiMVector dilepton = lplus + lminus;
    TLorentzVector dilep(dilepton.X(), dilepton.Y(), dilepton.Z(), dilepton.T());
    TLorentzVector boostedLep(lplus.X(), lplus.Y(), lplus.Z(), lplus.T());
    boostedLep.Boost(-1*dilep.BoostVector());
    const TVector3 csFrame = csframe(dilep);
    return cos(boostedLep.Angle(csFrame));
}

float phiCS(ROOT::Math::PtEtaPhiMVector lplus, ROOT::Math::PtEtaPhiMVector lminus) {
    ROOT::Math::PtEtaPhiMVector dilepton = lplus + lminus;
    TLorentzVector dilep(dilepton.X(), dilepton.Y(), dilepton.Z(), dilepton.T());
    TLorentzVector boostedLep(lplus.X(), lplus.Y(), lplus.Z(), lplus.T());
    boostedLep.Boost(-1*dilep.BoostVector());
    const TVector3 csFrameX = csframeX(dilep);
    const TVector3 csFrameY = csframeY(dilep);
    float phi = atan2(boostedLep.Vect()*csFrameY, boostedLep.Vect()*csFrameX);
    return phi >= 0 ? phi : phi + 2*M_PI;
}

void CSFrame() {
    TLorentzVector test;
    test.SetPtEtaPhiM(20., 1.5, 1., 90.);

    TVector3 cs = csframe(test);
    std::cout << "cs frame x,y,z " << cs.X() << " " 
              << cs.Y() << " " << cs.Z() << std::endl;
    TVector3 csx = csframeX(test);
    std::cout << "cs frame x: x,y,z " << csx.X() << " " 
              << csx.Y() << " " << csx.Z() << std::endl;
    TVector3 csy = csframeY(test);
    std::cout << "cs frame y: y,y,z " << csy.y() << " " 
              << csy.Y() << " " << csy.Z() << std::endl;

    ROOT::Math::PtEtaPhiMVector test1;
    test1.SetPt(28.5);
    test1.SetEta(2.7);
    test1.SetPhi(-2.7);
    test1.SetM(0.106);
    ROOT::Math::PtEtaPhiMVector test2;
    test1.SetPt(55.);
    test1.SetEta(1.68);
    test1.SetPhi(0.20);
    test1.SetM(.106);
    std::cout << "cosThetaCS = " << cosThetaCS(test1, test2) << std::endl;
    std::cout << "phiCS = " << phiCS(test1, test2) << std::endl;
}

