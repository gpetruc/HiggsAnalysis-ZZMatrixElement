#include "ZZMatrixElement/PythonWrapper/interface/MEMCalculatorsWrapper.h"

#include "ZZMatrixElement/MEMCalculators/interface/MEMCalculators.h"
#include "ZZMatrixElement/MELA/src/computeAngles.h"
#include "ZZMatrixElement/MELA/interface/Mela.h"

MEMCalculatorsWrapper::MEMCalculatorsWrapper(double collisionEnergy, double sKD_mass) {
    mem_ = new MEMs(collisionEnergy,sKD_mass,"CTEQ6L");
    mela_ = mem_->m_MELA;
}

MEMCalculatorsWrapper::~MEMCalculatorsWrapper() {
    if(mem_ !=0) delete mem_;
}

MEMCalculatorsWrapper::Angles 
MEMCalculatorsWrapper::computeAngles(TLorentzVector Z1_lept1, int Z1_lept1Id,
		   TLorentzVector Z1_lept2, int Z1_lept2Id,
		   TLorentzVector Z2_lept1, int Z2_lept1Id,
		   TLorentzVector Z2_lept2, int Z2_lept2Id) {
    Angles ret;
    mela::computeAngles(Z1_lept1,Z1_lept1Id,
                        Z1_lept2,Z1_lept2Id,
                        Z2_lept1,Z2_lept1Id,
                        Z2_lept2,Z2_lept2Id,
                    ret.costhetastar,ret.costheta1,ret.costheta2,ret.phi,ret.phistar1);
    return ret;
  }

MEMCalculatorsWrapper::Angles 
MEMCalculatorsWrapper::computeAngles(const math::XYZTLorentzVector & Z1_lept1, int Z1_lept1Id,
		   const math::XYZTLorentzVector & Z1_lept2, int Z1_lept2Id,
		   const math::XYZTLorentzVector & Z2_lept1, int Z2_lept1Id,
		   const math::XYZTLorentzVector & Z2_lept2, int Z2_lept2Id) {
    return computeAngles(TLorentzVector(Z1_lept1.Px(),Z1_lept1.Py(),Z1_lept1.Pz(),Z1_lept1.E()), Z1_lept1Id,
                         TLorentzVector(Z1_lept2.Px(),Z1_lept2.Py(),Z1_lept2.Pz(),Z1_lept2.E()), Z1_lept2Id,
                         TLorentzVector(Z2_lept1.Px(),Z2_lept1.Py(),Z2_lept1.Pz(),Z2_lept1.E()), Z2_lept1Id,
                         TLorentzVector(Z2_lept2.Px(),Z2_lept2.Py(),Z2_lept2.Pz(),Z2_lept2.E()), Z2_lept2Id);
  }


std::vector<std::pair<std::string,float>>  
MEMCalculatorsWrapper::computeNew(
        const math::XYZTLorentzVector & Z1_lept1, int Z1_lept1Id,
        const math::XYZTLorentzVector & Z1_lept2, int Z1_lept2Id,
        const math::XYZTLorentzVector & Z2_lept1, int Z2_lept1Id,
        const math::XYZTLorentzVector & Z2_lept2, int Z2_lept2Id,
        const std::vector<math::XYZTLorentzVector> & jets)
{
    std::vector<TLorentzVector> partP;
    partP.emplace_back(Z1_lept1.Px(),Z1_lept1.Py(),Z1_lept1.Pz(),Z1_lept1.E());
    partP.emplace_back(Z1_lept2.Px(),Z1_lept2.Py(),Z1_lept2.Pz(),Z1_lept2.E());
    partP.emplace_back(Z2_lept1.Px(),Z2_lept1.Py(),Z2_lept1.Pz(),Z2_lept1.E());
    partP.emplace_back(Z2_lept2.Px(),Z2_lept2.Py(),Z2_lept2.Pz(),Z2_lept2.E());

    std::vector<int> partId;
    partId.push_back(Z1_lept1Id);
    partId.push_back(Z1_lept2Id);
    partId.push_back(Z2_lept1Id);
    partId.push_back(Z2_lept2Id);

    double p0plus_VAJHU, bkg_VAMCFM, p0plus_m4l, bkg_m4l;
    double p0minus_VAJHU, Dgg10_VAMCFM;

    mem_->computeME(MEMNames::kSMHiggs, MEMNames::kJHUGen, partP, partId, p0plus_VAJHU); // Calculation of SM gg->H->4l JHUGen ME
    mem_->computeME(MEMNames::k0minus, MEMNames::kJHUGen, partP, partId, p0minus_VAJHU); // Calculation of PS (0-, fa3=1) gg->H->4l JHUGen ME 
    mem_->computeME(MEMNames::kggHZZ_10, MEMNames::kMCFM, partP, partId, Dgg10_VAMCFM); // Direct calculation of Dgg (D^kin for off-shell) from MCFM MEs
    mem_->computeME(MEMNames::kqqZZ, MEMNames::kMCFM, partP, partId, bkg_VAMCFM); // qq->4l background calculation from MCFM
    mem_->computePm4l(partP,partId, MEMNames::kNone, p0plus_m4l, bkg_m4l); // m4l probabilities for signal and background, nominal resolution

    double D_bkg_kin = p0plus_VAJHU / ( p0plus_VAJHU + bkg_VAMCFM ); // D^kin_bkg
    double D_bkg = p0plus_VAJHU * p0plus_m4l / ( p0plus_VAJHU * p0plus_m4l + bkg_VAMCFM * bkg_m4l ); // D^kin including superMELA
    double D_g4 = p0plus_VAJHU / ( p0plus_VAJHU + p0minus_VAJHU ); // D_0-
    std::vector<std::pair<std::string,float>> ret;
    ret.emplace_back("D_bkg^kin",D_bkg_kin);
    ret.emplace_back("D_bkg",D_bkg);
    ret.emplace_back("D_gg",Dgg10_VAMCFM);
    ret.emplace_back("D_0-", D_g4);

    TLorentzVector higgs_undec = partP[0]+partP[1]+partP[2]+partP[3];
    std::vector<TLorentzVector> partPprod(partP);
    std::vector<int>            partIdprod(partId);
    for (const auto &p4 : jets) {
        partPprod.emplace_back(p4.Px(), p4.Py(), p4.Pz(), p4.E());
        partIdprod.push_back(0);
        if (partPprod.size() == 6) break;
    }
    if (jets.size() >= 2) {
        double phjj_VAJHU, pvbf_VAJHU;
        float pwh_hadronic_VAJHU, pzh_hadronic_VAJHU;
        mem_->computeME(MEMNames::kJJ_SMHiggs_GG, MEMNames::kJHUGen,  partPprod, partIdprod, phjj_VAJHU); // SM gg->H+2j
        mem_->computeME(MEMNames::kJJ_SMHiggs_VBF, MEMNames::kJHUGen, partPprod, partIdprod, pvbf_VAJHU);  // SM VBF->H

        mela_->setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::WH);
        mela_->computeProdP(partPprod[4], partIdprod[4], partPprod[5], partIdprod[5], higgs_undec, pwh_hadronic_VAJHU);  // SM W(->2j)H
        mela_->setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::ZH);
        mela_->computeProdP(partPprod[4], partIdprod[4], partPprod[5], partIdprod[5], higgs_undec, pzh_hadronic_VAJHU); // SM Z(->2j)H

        double Djet_VAJHU = pvbf_VAJHU / ( pvbf_VAJHU + phjj_VAJHU ); // D^VBF_HJJ
        double D_WHh_VAJHU = pwh_hadronic_VAJHU / ( pwh_hadronic_VAJHU + 100000.*phjj_VAJHU ); // W(->2j)H vs. gg->H+2j
        double D_ZHh_VAJHU = pzh_hadronic_VAJHU / ( pzh_hadronic_VAJHU + 10000.*phjj_VAJHU ); // Z(->2j)H vs. gg->H+2j
        ret.emplace_back("D_HJJ^VBF", Djet_VAJHU);
        ret.emplace_back("D_HJJ^WH", D_WHh_VAJHU);
        ret.emplace_back("D_HJJ^ZH", D_ZHh_VAJHU);
    } else if (jets.size() == 1) {
        TLorentzVector nullFourVector(0, 0, 0, 0);
        float pvbf_VAJHU, pAux_vbf_VAJHU, phj_VAJHU;
        mela_->setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::JH);
        mela_->computeProdP(partPprod[4], partIdprod[4], nullFourVector, 0, higgs_undec, phj_VAJHU); // SM gg->H+1j
        mela_->setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::JJVBF);
        mela_->computeProdP(partPprod[4], partIdprod[4], nullFourVector, 0, higgs_undec, pvbf_VAJHU); // Un-integrated ME
        mela_->get_PAux(pAux_vbf_VAJHU); // = Integrated / un-integrated
        double D_VBF1j_VAJHU = pvbf_VAJHU*pAux_vbf_VAJHU / ( pvbf_VAJHU*pAux_vbf_VAJHU + 5.*phj_VAJHU ); // VBF(1j) vs. gg->H+1j
        ret.emplace_back("D_HJJ^VBF", D_VBF1j_VAJHU);
        ret.emplace_back("D_HJJ^WH", -1.0);
        ret.emplace_back("D_HJJ^ZH", -1.0);
    } else {
        ret.emplace_back("D_HJJ^VBF", -1.0);
        ret.emplace_back("D_HJJ^WH", -1.0);
        ret.emplace_back("D_HJJ^ZH", -1.0);
    }

    return ret;
}

