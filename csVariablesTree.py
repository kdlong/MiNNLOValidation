import ROOT
df_read = ROOT.RDataFrame("Events", "/eos/cms/store/user/kelong/ZJToMuMu_mWPilot_TuneCP5_13TeV-powheg-MiNNLO-pythia8-photos/NanoGen/200429_162601/0000/*")

tmpdf = df_read.Define("pmu1", "ROOT::Math::PtEtaPhiMVector(LHEPart_pt[0], LHEPart_eta[0], LHEPart_phi[0], LHEPart_mass[0])")
tmpdf = tmpdf.Define("pmu2", "ROOT::Math::PtEtaPhiMVector(LHEPart_pt[1], LHEPart_eta[1], LHEPart_phi[1], LHEPart_mass[1])")
newdf = tmpdf.Define("pZ", "pmu1+pmu2")
dfz = newdf.Define("ptZ", "pZ.Pt()")
dfz = dfz.Define("yZ", "pZ.Rapidity()")
dfz = dfz.Define("absyZ", "fabs(pZ.Rapidity())")
dfz = dfz.Define("phiZ", "pZ.Phi()")
dfz = dfz.Define("mZ", "pZ.M()")
dfz = dfz.Filter("mZ < 101 && mZ > 81")

ROOT.gInterpreter.Declare("""
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
""")

tmpdf2 = dfz.Define("pmup", "LHEPart_pdgId[0] < 0 ? pmu1 : pmu2")
tmpdf2 = tmpdf2.Define("pmum", "LHEPart_pdgId[0] > 0 ? pmu1 : pmu2")
tmpdf2 = tmpdf2.Define("pmumpt", "pmum.Pt()")

dfcs = tmpdf2.Define("costcs", "cosThetaCS(pmup, pmum)")
dfcs = dfcs.Define("phics", "phiCS(pmup, pmum)")

cols  = ROOT.RDFDetail.ColumnNames_t()
for c in ["costcs", "phics", "yZ", "ptZ", "mZ", "phiZ", "genWeight"]:
    cols.push_back(c)
dfcs.Snapshot("Events", "newtest.root", cols)

