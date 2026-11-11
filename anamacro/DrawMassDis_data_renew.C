using namespace std;

#include <optional>
#include <cmath>

struct MatchedTrkCalo
{
    int evt;//event id
    int itrk;//trk index
    int icalo;//calo index
    // track info
    float trk_pt;
    float trk_px;
    float trk_py;
    float trk_pz;
    float trk_eta;
    float trk_phi;
    int trk_charge;
    float trk_x;
    float trk_y;
    float trk_z;
    int nmaps;
    int nintt;
    float chi2ndf;
    // track projection at EMCal
    float trk_proj_phi;
    float trk_proj_z;
    // calo cluster info
    float calo_phi;
    float calo_z;
    float calo_E;

    // for calo matching before correction
    float dphi; // track_phi_emc - calo_phi
    float dz;   // z_emc_trk   - calo_z
    // E/p
    float eop;
    float score;
    //vertex
    float zvtx;

};

struct PairMass
{
    int evt;
    int i_1; // index in matched_this_evt
    int i_2;

    float reconst_mass;

    int charge_1, charge_2;
    float pt_1, pt_2;
    float px_1, px_2;
    float py_1, py_2;
    float pz_1, pz_2;
    float trk_z_1, trk_z_2;
    float eta_1, eta_2;
    float phi_1, phi_2;
    float dphi_1, dphi_2;
    float dz_1,dz_2;
    int nmaps_1, nmaps_2;
    int nintt_1, nintt_2;
    float chi2ndf_1, chi2ndf_2;
    float zvtx_1, zvtx_2;
    float eop_1, eop_2;

};

struct PassFlags
{
    bool pass_pt;
    bool pass_z0;
    bool pass_hits;
    bool pass_chi2ndf;
    bool pass_dz;
    bool pass_dphi;
    bool pass_trk_z;
    bool pass_zvtx;
    bool pass_dphidz;
    bool pass_all;
};
//koko
PassFlags cutFlags(
            const MatchedTrkCalo &mt,
            double PtCut,
            double trk_zMinCut,
            double trk_zMaxCut,
            double Ztrk_zvtxCut,
            int NmapsCut,
            int NinttCut,
            double Chi2ndfCut,
            double ZvtxCut,
            double DphiCut,
            double DzCut
            )
{
    PassFlags fl;
    fl.pass_pt = (mt.trk_pt >= PtCut);
    fl.pass_z0 = fabs(mt.trk_z - mt.zvtx) < Ztrk_zvtxCut;
    fl.pass_zvtx = fabs(mt.zvtx) <= ZvtxCut;
    fl.pass_hits = (mt.nmaps == NmapsCut && mt.nintt >= NinttCut);
    fl.pass_dphi = fabs(mt.dphi) < DphiCut;
    fl.pass_dz = fabs(mt.dz) < DzCut;
    fl.pass_dphidz = (fl.pass_dz && fl.pass_dphi);
    fl.pass_chi2ndf = mt.chi2ndf <= Chi2ndfCut;
    fl.pass_trk_z = (mt.trk_z > trk_zMinCut && mt.trk_z < trk_zMaxCut);
    fl.pass_all = (fl.pass_pt && fl.pass_z0 && fl.pass_zvtx && fl.pass_hits && fl.pass_dphi && fl.pass_dz && fl.pass_chi2ndf && fl.pass_trk_z);

    return fl;
}

float wrapPi(float dphi)
{
    while (dphi >= M_PI)
        dphi -= 2 * M_PI;
    while (dphi < -M_PI)
        dphi += 2 * M_PI;
    return dphi;
}

std::optional<MatchedTrkCalo> matchTrackCalo(
                                float trk_proj_x,
                                float trk_proj_y,
                                float trk_proj_z,
                                vector<float>& calo_phi_,
                                vector<float>& calo_z_,
                                double Wphi,
                                double r_emc
                                )
{   
    float min_score = 1e9;
    int min_icalo = -1;
    float min_dphi = 0.0;
    float min_dz = 0.0;
    double trk_proj_phi = atan2(trk_proj_y, trk_proj_x);

    for(size_t icalo = 0; icalo < calo_phi_.size(); icalo++)
    {   
        float dphi = wrapPi(trk_proj_phi - calo_phi_[icalo]);
        float dphi_cm = dphi * r_emc;
        float dz = trk_proj_z - calo_z_[icalo];

        double score = sqrt(dz * dz + (dphi_cm * dphi_cm) / Wphi);

        if( score < min_score)
        {
            min_score = score;
            min_icalo = icalo;
            min_dphi = dphi;
            min_dz = dz;
        }
    }

    if (min_icalo < 0)return std::nullopt;//cannot find matching
    {
        MatchedTrkCalo result;
        result.icalo = min_icalo;
        result.trk_proj_phi = trk_proj_phi;
        result.trk_proj_z   = trk_proj_z;

        result.calo_phi = calo_phi_[min_icalo];
        result.calo_z   = calo_z_[min_icalo];

        result.dphi = min_dphi;
        result.dz   = min_dz;
        result.score = min_score;

        return result;
    }; 
}

std::vector<MatchedTrkCalo> AllmatchTrackCalo(
                                    float trk_proj_x,
                                    float trk_proj_y,
                                    float trk_proj_z,
                                    const std::vector<float>& calo_phi_,
                                    const std::vector<float>& calo_z_,
                                    const std::vector<float>& calo_E_,
                                    double Wphi,
                                    double r_emc
                                )
{
    std::vector<MatchedTrkCalo> results;

    double trk_proj_phi = atan2(trk_proj_y, trk_proj_x);

    for (size_t icalo = 0;icalo < calo_phi_.size();icalo++)
    {
        float dphi = wrapPi(trk_proj_phi - calo_phi_[icalo]);
        float dphi_cm = dphi * r_emc;
        float dz = trk_proj_z - calo_z_[icalo];
        float score = sqrt( dz * dz + (dphi_cm * dphi_cm) / Wphi);

        MatchedTrkCalo result;

        result.icalo    = icalo;
        result.calo_phi = calo_phi_[icalo];
        result.calo_z   = calo_z_[icalo];
        result.calo_E   = calo_E_[icalo];

        result.trk_proj_phi = trk_proj_phi;
        result.trk_proj_z   = trk_proj_z;

        result.dphi  = dphi;
        result.dz    = dz;
        result.score = score;

        results.push_back(result);
    }

    return results;
}

std::vector<PairMass> mass_reconstruction( 
                        vector<MatchedTrkCalo>& matched_this_evt,
                        float emass 
                        )
{
    vector<PairMass> results;
    size_t n_matched = matched_this_evt.size();
    if(n_matched < 2)return results;//track数が２本未満

    for(size_t i_1 = 0; i_1 < n_matched; i_1++)
    {
        MatchedTrkCalo &trk_1 = matched_this_evt[i_1];
        float charge_1 = trk_1.trk_charge;
        float px_1 = trk_1.trk_px;
        float py_1 = trk_1.trk_py;
        float pz_1 = trk_1.trk_pz;
        float p_1 = sqrt(px_1*px_1 + py_1*py_1 + pz_1*pz_1);
        float E_1 = sqrt(p_1*p_1 + emass*emass);

        for(size_t i_2 = i_1 +1; i_2 < n_matched; i_2++)
        {
            MatchedTrkCalo &trk_2 = matched_this_evt[i_2];
            float charge_2 = trk_2.trk_charge;
            if(charge_1*charge_2 >= 0)continue;//++ or --

            float px_2 = trk_2.trk_px;
            float py_2 = trk_2.trk_py;
            float pz_2 = trk_2.trk_pz;
            float p_2 = sqrt(px_2*px_2 + py_2*py_2 + pz_2*pz_2);
            float E_2 = sqrt(p_2*p_2 + emass*emass);

            float px = px_1 + px_2;
            float py = py_1 + py_2;
            float pz = pz_1 + pz_2;
            float p_squared = px*px + py*py + pz*pz;

            float reconst_mass = sqrt((E_1 + E_2)*(E_1 + E_2) - p_squared);
            //cout<<reconst_mass<<endl;
            PairMass result;

            result.i_1 = i_1;
            result.i_2 = i_2;
            result.charge_1 = trk_1.trk_charge;
            result.charge_2 = trk_2.trk_charge;
            result.reconst_mass = reconst_mass;

            result.pt_1 = trk_1.trk_pt;
            result.pt_2 = trk_2.trk_pt;
            result.px_1 = trk_1.trk_px;
            result.px_2 = trk_2.trk_px;
            result.py_1 = trk_1.trk_py;
            result.py_2 = trk_2.trk_py;
            result.pz_1 = trk_1.trk_pz;
            result.pz_2 = trk_2.trk_pz;
            result.trk_z_1 = trk_1.trk_z;
            result.trk_z_2 = trk_2.trk_z;

            result.nmaps_1 = trk_1.nmaps;
            result.nmaps_2 = trk_2.nmaps;
            result.nintt_1 = trk_1.nintt;
            result.nintt_2 = trk_2.nintt;
            result.eta_1 = trk_1.trk_eta;
            result.eta_2 = trk_2.trk_eta;
            result.phi_1 = trk_1.trk_phi;
            result.phi_2 = trk_2.trk_phi;
            result.trk_z_1 = trk_1.trk_z;
            result.trk_z_2 = trk_2.trk_z;
            result.chi2ndf_1 = trk_1.chi2ndf;
            result.chi2ndf_2 = trk_2.chi2ndf;

            result.dphi_1 = trk_1.dphi;
            result.dphi_2 = trk_2.dphi;
            result.dz_1 = trk_1.dz;
            result.dz_2 = trk_2.dz;
            result.zvtx_1 = trk_1.zvtx;
            result.zvtx_2 = trk_2.zvtx;
            result.eop_1 = trk_1.calo_E / trk_1.trk_pt * std::cosh(trk_1.trk_eta);
            result.eop_2 = trk_2.calo_E / trk_2.trk_pt * std::cosh(trk_2.trk_eta);

            results.push_back(result);
        }
    }
    return results;//+-のpairがない
}

void DrawMassDis_data_renew(
    //const string &inputfile = "/sphenix/u/mizukishigihara/tgmizuki/Jpsi/RunQA/cpdata/ana_53877_00001.root"
    //const string &inputfile = "/sphenix/tg/tg01/commissioning/INTT/work/hachiya/SiCaloTrack/data/run2pp/53741/ana_53741_00000.root"
    const string &inputfile = "/sphenix/tg/tg01/commissioning/INTT/work/hachiya/SiCaloTrack/data/run3pp/79508/ana_79508_00000.root"
)
{
    //path
    //const string outdir = "/sphenix/tg/tg01/commissioning/INTT/work/mizuki/Jpsi_renew/anamacro/result/data_run3";
    //const string outdir = "/sphenix/tg/tg01/commissioning/INTT/work/mizuki/Jpsi_renew/anamacro/result/data";
    const string outdir_allmatched = "/sphenix/tg/tg01/commissioning/INTT/work/mizuki/Jpsi_renew/anamacro/result/data_run3_allmatched";
    //const string outdir_allmatched = "/sphenix/tg/tg01/commissioning/INTT/work/mizuki/Jpsi_renew/anamacro/result/data_allmatched";
    //変数
    double Wphi = 3.31;
    double r_emc =93.5;
    float electron_mass = 0.000511;//[Gev/c^2]
    //cut値
    const double PtCut = 0.5;
    const double DphiCut = 0.178;
    const double DzCut = 5.016;
    const double trk_zMinCut = 0.8;
    const double trk_zMaxCut = 1.2;
    const double Ztrk_zvtxCut = 1;
    const int NmapsCut = 3;
    const int NinttCut = 1;
    const double Chi2ndfCut = 22.781;
    const double ZvtxCut = 10.0;
    
    //input file
    TFile *f = TFile::Open(inputfile.c_str(),"READ");
    if((!f ||f->IsZombie())){
        std::cerr <<"Error: cannot open file: "<< inputfile.c_str() << std::endl;
        return;
    }
    // //outputfile
    // TString fname_out = gSystem->BaseName(inputfile.c_str());
    // fname_out.ReplaceAll("ana", "mass_reco");
    // TFile *fout = new TFile(Form("%s/%s", outdir.c_str(), fname_out.Data()),"RECREATE");
    //outputfile_allmached
    TString fname_out_allmatched = gSystem->BaseName(inputfile.c_str());
    fname_out_allmatched.ReplaceAll("ana", "allmatched");
    TFile *fout_allmatched = new TFile(Form("%s/%s", outdir_allmatched.c_str(), fname_out_allmatched.Data()),"RECREATE");
    //TFile *fout = new TFile(Form("%s/test.root", outdir.c_str()),"RECREATE");


    //output TTree MatchedSiEMC
    // TTree *MatchedSiEMC = nullptr;
    // MatchedSiEMC = new TTree("MatchedSiEMC", "Matched track-emc info");
    // MatchedTrkCalo mt;

    // //evt
    // MatchedSiEMC->Branch("evt", &mt.evt, "evt/I");
    // MatchedSiEMC->Branch("itrk", &mt.itrk, "itrk/I");
    // MatchedSiEMC->Branch("icalo", &mt.icalo, "icalo/I");
    // // track info
    // MatchedSiEMC->Branch("trk_pt", &mt.trk_pt, "trk_pt/F");
    // MatchedSiEMC->Branch("trk_px", &mt.trk_px, "trk_px/F");
    // MatchedSiEMC->Branch("trk_py", &mt.trk_py, "trk_py/F");
    // MatchedSiEMC->Branch("trk_pz", &mt.trk_pz, "trk_pz/F");
    // MatchedSiEMC->Branch("trk_eta", &mt.trk_eta, "trk_eta/F");
    // MatchedSiEMC->Branch("trk_phi", &mt.trk_phi, "trk_phi/F");
    // MatchedSiEMC->Branch("trk_charge", &mt.trk_charge, "trk_charge/I");
    // MatchedSiEMC->Branch("trk_z", &mt.trk_z, "trk_z/F");
    // MatchedSiEMC->Branch("nmaps", &mt.nmaps, "nmaps/I");
    // MatchedSiEMC->Branch("nintt", &mt.nintt, "nintt/I");
    // MatchedSiEMC->Branch("chi2ndf", &mt.chi2ndf, "chi2ndf/F");
    // // track projection si-emc info
    // MatchedSiEMC->Branch("trk_proj_phi", &mt.trk_proj_phi, "trk_proj_phi/F");
    // MatchedSiEMC->Branch("trk_proj_z", &mt.trk_proj_z, "trk_proj_z/F");
    // //calo cluster info
    // MatchedSiEMC->Branch("calo_phi", &mt.calo_phi, "calo_phi/F");
    // MatchedSiEMC->Branch("calo_z", &mt.calo_z, "calo_z/F");
    // MatchedSiEMC->Branch("calo_E", &mt.calo_E, "calo_E/F");
    // //matching dphi and dz
    // MatchedSiEMC->Branch("dphi", &mt.dphi, "dphi/F");
    // MatchedSiEMC->Branch("dz", &mt.dz, "dz/F");
    // //E/p
    // MatchedSiEMC->Branch("eop", &mt.eop, "eop/F");
    // MatchedSiEMC->Branch("score",        &mt.score,        "score/F");
    // //vertex
    // MatchedSiEMC->Branch("zvtx", &mt.zvtx, "zvtx/F");

     //output TTree alldata
    TTree *AllMatchedSiEMC = new TTree("AllMatchedSiEMC", "All track-calo matched data");
    MatchedTrkCalo ac;

    //evt
    AllMatchedSiEMC->Branch("evt", &ac.evt, "evt/I");
    AllMatchedSiEMC->Branch("itrk", &ac.itrk, "itrk/I");
    AllMatchedSiEMC->Branch("icalo", &ac.icalo, "icalo/I");
    // track info
    AllMatchedSiEMC->Branch("trk_pt", &ac.trk_pt, "trk_pt/F");
    AllMatchedSiEMC->Branch("trk_px", &ac.trk_px, "trk_px/F");
    AllMatchedSiEMC->Branch("trk_py", &ac.trk_py, "trk_py/F");
    AllMatchedSiEMC->Branch("trk_pz", &ac.trk_pz, "trk_pz/F");
    AllMatchedSiEMC->Branch("trk_eta", &ac.trk_eta, "trk_eta/F");
    AllMatchedSiEMC->Branch("trk_phi", &ac.trk_phi, "trk_phi/F");
    AllMatchedSiEMC->Branch("trk_charge", &ac.trk_charge, "trk_charge/I");
    AllMatchedSiEMC->Branch("trk_z", &ac.trk_z, "trk_z/F");
    AllMatchedSiEMC->Branch("nmaps", &ac.nmaps, "nmaps/I");
    AllMatchedSiEMC->Branch("nintt", &ac.nintt, "nintt/I");
    AllMatchedSiEMC->Branch("chi2ndf", &ac.chi2ndf, "chi2ndf/F");
    // track projection si-emc info
    AllMatchedSiEMC->Branch("trk_proj_phi", &ac.trk_proj_phi, "trk_proj_phi/F");
    AllMatchedSiEMC->Branch("trk_proj_z", &ac.trk_proj_z, "trk_proj_z/F");
    //calo cluster info
    AllMatchedSiEMC->Branch("calo_phi", &ac.calo_phi, "calo_phi/F");
    AllMatchedSiEMC->Branch("calo_z", &ac.calo_z, "calo_z/F");
    AllMatchedSiEMC->Branch("calo_E", &ac.calo_E, "calo_E/F");
    //matching dphi and dz
    AllMatchedSiEMC->Branch("dphi", &ac.dphi, "dphi/F");
    AllMatchedSiEMC->Branch("dz", &ac.dz, "dz/F");
    //E/p
    AllMatchedSiEMC->Branch("eop", &ac.eop, "eop/F");
    AllMatchedSiEMC->Branch("score", &ac.score,"score/F");
    //vertex
    AllMatchedSiEMC->Branch("zvtx", &ac.zvtx, "zvtx/F");

    //output TTree recons    
    // TTree *ReconstMass = new TTree("ReconstMass", "e+e- pair mass reconstruction");
    // PairMass pm;

    // ReconstMass->Branch("evt", &pm.evt, "evt/I");
    // ReconstMass->Branch("i_1", &pm.i_1, "i_1/I");
    // ReconstMass->Branch("i_2", &pm.i_2, "i_2/I");
    // ReconstMass->Branch("reconst_mass", &pm.reconst_mass, "reconst_mass/F");

    // ReconstMass->Branch("charge_1", &pm.charge_1, "charge_1/I");
    // ReconstMass->Branch("charge_2", &pm.charge_2, "charge_2/I");
    // ReconstMass->Branch("pt_1", &pm.pt_1, "pt_1/F");
    // ReconstMass->Branch("pt_2", &pm.pt_2, "pt_2/F");
    // ReconstMass->Branch("px_1", &pm.px_1, "px_1/F");
    // ReconstMass->Branch("px_2", &pm.px_2, "px_2/F");
    // ReconstMass->Branch("py_1", &pm.py_1, "py_1/F");
    // ReconstMass->Branch("py_2", &pm.py_2, "py_2/F");
    // ReconstMass->Branch("pz_1", &pm.pz_1, "pz_1/F");
    // ReconstMass->Branch("pz_2", &pm.pz_2, "pz_2/F");
    // ReconstMass->Branch("trk_z_1", &pm.trk_z_1, "trk_z_1/F");
    // ReconstMass->Branch("trk_z_2", &pm.trk_z_2, "trk_z_2/F");
    // ReconstMass->Branch("eta_1", &pm.eta_1, "eta_1/F");
    // ReconstMass->Branch("eta_2", &pm.eta_2, "eta_2/F");
    // ReconstMass->Branch("phi_1", &pm.phi_1, "phi_1/F");
    // ReconstMass->Branch("phi_2", &pm.phi_2, "phi_2/F");
    // ReconstMass->Branch("trk_z_1", &pm.trk_z_1, "trk_z_1/F");
    // ReconstMass->Branch("trk_z_2", &pm.trk_z_2, "trk_z_2/F");
    // ReconstMass->Branch("dphi_1", &pm.dphi_1, "dphi_1/F");
    // ReconstMass->Branch("dphi_2", &pm.dphi_2, "dphi_2/F");
    // ReconstMass->Branch("dz_1", &pm.dz_1, "dz_1/F");
    // ReconstMass->Branch("dz_2", &pm.dz_2, "dz_2/F");
    // ReconstMass->Branch("nmaps_1", &pm.nmaps_1, "nmaps_1/I");
    // ReconstMass->Branch("nmaps_2", &pm.nmaps_2, "nmaps_2/I");
    // ReconstMass->Branch("nintt_1", &pm.nintt_1, "nintt_1/I");
    // ReconstMass->Branch("nintt_2", &pm.nintt_2, "nintt_2/I");
    // ReconstMass->Branch("chi2ndf_1", &pm.chi2ndf_1, "chi2ndf_1/F");
    // ReconstMass->Branch("chi2ndf_2", &pm.chi2ndf_2, "chi2ndf_2/F");
    // ReconstMass->Branch("zvtx_1", &pm.zvtx_1, "zvtx_1/F");
    // ReconstMass->Branch("zvtx_2", &pm.zvtx_2, "zvtx_2/F");
    // ReconstMass->Branch("eop_1", &pm.eop_1, "eop_1/F");
    // ReconstMass->Branch("eop_2", &pm.eop_2, "eop_2/F");

    //input TTree
    TTree *trackTree = (TTree *)f->Get("trackTree");
    TTree *caloTree = (TTree *)f->Get("caloTree");
    TTree *evtTree = (TTree *)f->Get("evtTree");

    //input branches
    int trk_evt = 0;
    vector<float> *trk_proj_x_ = nullptr;
    vector<float> *trk_proj_y_ = nullptr;
    vector<float> *trk_proj_z_ = nullptr;
    vector<float> *trk_pt_ = nullptr;
    vector<float> *trk_px_ = nullptr;
    vector<float> *trk_py_ = nullptr;
    vector<float> *trk_pz_ = nullptr;
    vector<float> *trk_eta_ = nullptr;
    vector<int> *trk_charge_ = nullptr;
    vector<float> *trk_x_ =  nullptr;
    vector<float> *trk_y_ =  nullptr;
    vector<float> *trk_z_ =  nullptr;

    vector<int> *trk_nmaps_ = nullptr;
    vector<int> *trk_nintt_ = nullptr;
    vector<float> *trk_chi2ndf_ = nullptr;
    trackTree->SetBranchAddress("evt", &trk_evt);
    trackTree->SetBranchAddress("x_proj_emc", &trk_proj_x_);
    trackTree->SetBranchAddress("y_proj_emc", &trk_proj_y_);
    trackTree->SetBranchAddress("z_proj_emc", &trk_proj_z_);
    trackTree->SetBranchAddress("pt0", &trk_pt_);
    trackTree->SetBranchAddress("px0", &trk_px_);
    trackTree->SetBranchAddress("py0", &trk_py_);
    trackTree->SetBranchAddress("pz0", &trk_pz_);
    trackTree->SetBranchAddress("eta0", &trk_eta_);
    trackTree->SetBranchAddress("charge", &trk_charge_);
    trackTree->SetBranchAddress("x0", &trk_x_);
    trackTree->SetBranchAddress("y0", &trk_y_);
    trackTree->SetBranchAddress("z0", &trk_z_);
    trackTree->SetBranchAddress("nmaps", &trk_nmaps_);
    trackTree->SetBranchAddress("nintt", &trk_nintt_);
    trackTree->SetBranchAddress("chi2ndf", &trk_chi2ndf_);

    int calo_evt = 0;
    vector<float> *calo_phi_ = nullptr;
    vector<float> *calo_z_ = nullptr;
    vector<float> *calo_E_ = nullptr;
    caloTree->SetBranchAddress("calo_evt", &calo_evt);
    caloTree->SetBranchAddress("phi", &calo_phi_);
    caloTree->SetBranchAddress("z", &calo_z_);
    caloTree->SetBranchAddress("energy", &calo_E_);

    float_t zvtx = 0.0;
    evtTree->SetBranchAddress("zvtx", &zvtx);

    //hist
    TH1F *h_mass = new TH1F("h_mass","e^{+}e^{-} invariant mass;Mass [GeV/c^{2}];Counts",100,0.0,5.0);
    TH2F *h2_trkcalo_phi = new TH2F("h2_trkcalo_phi","proj_phi-calo_phi; trk_proj[rad]; calo[rad]",100,-3.2,3.2,100,-3.2,3.2);
    TH2F *h2_trkcalo_z = new TH2F("h2_trkcalo_z","proj_z-calo_z; trk_proj[cm]; calo[cm]",100,-150,150,100,-150,150);

    //Si-EMC matching
    Long64_t nentries = trackTree->GetEntries();
    for(Long64_t ievt = 0; ievt < nentries; ievt++)
    {
        vector<MatchedTrkCalo> matched_this_evt;
        //vector<PassFlags> flags_this_evt;
        trackTree->GetEntry(ievt);
        caloTree->GetEntry(ievt);
        evtTree->GetEntry(ievt);

        for(size_t itrk = 0; itrk < trk_proj_x_->size(); itrk++){

            //double trk_proj_phi = atan2((*trk_proj_y_)[itrk], (*trk_proj_x_)[itrk]);
            float trk_proj_x = (*trk_proj_x_)[itrk];
            float trk_proj_y = (*trk_proj_y_)[itrk];
            float trk_proj_z = (*trk_proj_z_)[itrk];

            auto Allmatched_trk = AllmatchTrackCalo(trk_proj_x, trk_proj_y, trk_proj_z, *calo_phi_, *calo_z_, *calo_E_, Wphi, r_emc);
             for (const auto& result : Allmatched_trk)
            {   
                ac = result;
                ac.evt       = trk_evt;
                ac.itrk      = itrk;
                ac.trk_pt    = (*trk_pt_)[itrk];
                ac.trk_px    = (*trk_px_)[itrk];
                ac.trk_py    = (*trk_py_)[itrk];
                ac.trk_pz    = (*trk_pz_)[itrk];
                ac.trk_eta   = (*trk_eta_)[itrk];
                ac.trk_charge = (*trk_charge_)[itrk];
                ac.trk_x = (*trk_x_)[itrk];
                ac.trk_y = (*trk_y_)[itrk];
                ac.trk_z = (*trk_z_)[itrk];
                ac.trk_phi = atan2((*trk_y_)[itrk], (*trk_x_)[itrk]);
                ac.eop = result.calo_E / (*trk_pt_)[itrk] * std::cosh((*trk_eta_)[itrk]);
                ac.zvtx = zvtx;
                ac.nmaps = (*trk_nmaps_)[itrk];
                ac.nintt = (*trk_nintt_)[itrk];
                ac.chi2ndf = (*trk_chi2ndf_)[itrk];
                AllMatchedSiEMC->Fill();

                h2_trkcalo_phi->Fill(ac.trk_proj_phi,ac.calo_phi);
                h2_trkcalo_z->Fill(ac.trk_proj_z,ac.calo_z);
            }
        
            // auto matched_trk = matchTrackCalo(trk_proj_x, trk_proj_y, trk_proj_z, *calo_phi_, *calo_z_, Wphi,r_emc);
            // if(!matched_trk)continue;
            // mt = *matched_trk;
           
            // mt.evt       = trk_evt;
            // mt.itrk      = itrk;
            // mt.trk_pt    = (*trk_pt_)[itrk];
            // mt.trk_px    = (*trk_px_)[itrk];
            // mt.trk_py    = (*trk_py_)[itrk];
            // mt.trk_pz    = (*trk_pz_)[itrk];
            // mt.trk_eta   = (*trk_eta_)[itrk];
            // mt.trk_charge = (*trk_charge_)[itrk];
            // mt.trk_x = (*trk_x_)[itrk];
            // mt.trk_y = (*trk_y_)[itrk];
            // mt.trk_z = (*trk_z_)[itrk];
            // mt.trk_phi = atan2((*trk_y_)[itrk], (*trk_x_)[itrk]);
            // mt.eop = (*calo_E_)[mt.icalo] / (*trk_pt_)[itrk] * std::cosh((*trk_eta_)[itrk]);
            // mt.zvtx = zvtx;
            // mt.nmaps = (*trk_nmaps_)[itrk];
            // mt.nintt = (*trk_nintt_)[itrk];
            // mt.chi2ndf = (*trk_chi2ndf_)[itrk];
            // AllMatchedSiEMC->Fill();
            // matched_this_evt.push_back(mt);

            // MatchedSiEMC->Fill();

            //PassFlags fl = cutFlags(mt,PtCut,trk_zMinCut,trk_zMaxCut, Ztrk_zvtxCut, NmapsCut, NinttCut, Chi2ndfCut, ZvtxCut, DphiCut, DzCut);
            //flags_this_evt.push_back(fl);
        }

       //mass reconstruction
    //    vector<PairMass> reconstructed_mass = mass_reconstruction(matched_this_evt, electron_mass);
    //    if(reconstructed_mass.empty())continue;

    //    for (const auto &result : reconstructed_mass)
    //     {
    //         pm = result;
    //         ReconstMass->Fill();
    //         h_mass->Fill(pm.reconst_mass);
    //     }
    //    //cout << "mass = " << pm.reconst_mass << endl;
    }
        fout_allmatched->cd();
        AllMatchedSiEMC->Write();
        fout_allmatched->Close();

        // fout->cd();
        // MatchedSiEMC->Write();
        // // ReconstMass->Write();
        // // h_mass->Write();
        // h2_trkcalo_phi->Write();
        // h2_trkcalo_z->Write();
        // fout->Close();

        f->Close();
}
