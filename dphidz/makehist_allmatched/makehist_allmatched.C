
// Return false for noisy calorimeter regions.
// Histograms without "_cut" are filled before this cut;
// histograms with "_cut" are filled only if this function returns true.
bool PassCaloCut(float calo_phi, float calo_z)
{
    if (calo_phi > -3.5 && calo_phi < -3.0 && calo_z > 130) return false;
    if (calo_phi > -1.7 && calo_phi < -1.4 && calo_z > 50)  return false;
    if (calo_phi > -1.4 && calo_phi < -1.2 && calo_z > 130) return false;
    if (calo_phi > -0.2 && calo_phi < 0.2 && calo_z > 50)   return false;
    if (calo_phi > 1.4 && calo_phi < 1.8 && calo_z > 50)    return false;
    if (calo_phi > 3.0 && calo_phi < 3.5 && calo_z > 50)    return false;
    return true;
}
bool PassDphiCut(float dphi)
{
    if(dphi > -0.4 && dphi < 0.4 )return true;
    return false;
}
bool PassDzCut(float dz)
{
    if(dz > -30 && dz < 30 )return true;
    return false;
}

void makehist_allmatched(
    const string inputfile =
    "/sphenix/u/mizukishigihara/tgmizuki/Jpsi_renew/anamacro/result/data_merged/53741_allmatched_merged.root")
{
    TFile *f = TFile::Open(inputfile.c_str(), "READ");
    if (!f || f->IsZombie()) {
        cerr << "Error: cannot open input file: " << inputfile << endl;
        if (f) { f->Close(); delete f; }
        return;
    }

    TTree *tree = dynamic_cast<TTree*>(f->Get("AllMatchedSiEMC"));
    if (!tree) {
        cerr << "Error: TTree AllMatchedSiEMC was not found in " << inputfile << endl;
        f->Close();
        delete f;
        return;
    }

    const string outdir =
        "/sphenix/tg/tg01/commissioning/INTT/work/mizuki/Jpsi_renew/dphidz/makehist_allmatched/result";
    gSystem->mkdir(outdir.c_str(), true);

    TString outname = gSystem->BaseName(inputfile.c_str());
    outname.ReplaceAll("allmatched_merged", "makehist_3cut");
    TString outfile = Form("%s/%s", outdir.c_str(), outname.Data());

    TFile *fout = TFile::Open(outfile, "RECREATE");
    if (!fout || fout->IsZombie()) {
        cerr << "Error: cannot create output file: " << outfile << endl;
        if (fout) { fout->Close(); delete fout; }
        f->Close();
        delete f;
        return;
    }

    float trk_pt = 0;
    int trk_charge = 0;
    float dphi = 0, dz = 0;
    float calo_phi = 0, calo_z = 0;
    float trk_proj_phi = 0, trk_proj_z = 0;

    tree->SetBranchAddress("trk_pt", &trk_pt);
    tree->SetBranchAddress("trk_charge", &trk_charge);
    tree->SetBranchAddress("dphi", &dphi);
    tree->SetBranchAddress("dz", &dz);
    tree->SetBranchAddress("calo_phi", &calo_phi);
    tree->SetBranchAddress("calo_z", &calo_z);
    tree->SetBranchAddress("trk_proj_phi", &trk_proj_phi);
    tree->SetBranchAddress("trk_proj_z", &trk_proj_z);

    // Every histogram is created as a pair:
    // base = no PassCaloCut requirement; _cut = PassCaloCut passed.
    // phi axes: 300 bins, -3.5 to 3.5 rad.
    // z axes:   400 bins, -200 to 200 cm.
    #define MAKE2(name, title, nx, xmin, xmax, ny, ymin, ymax) \
        TH2F *name = new TH2F(#name, title, nx, xmin, xmax, ny, ymin, ymax)

    MAKE2(h2_calohit, "Calo hit;calo_phi [rad];calo_z [cm]", 300, -3.5, 3.5, 400, -200, 200);
    MAKE2(h2_calohit_cut, "Calo hit after PassCaloCut;calo_phi [rad];calo_z [cm]", 300, -3.5, 3.5, 400, -200, 200);

    MAKE2(h2_Alldphidz, "All matched;dphi [rad];dz [cm]", 300, -3.5, 3.5, 400, -200, 200);
    MAKE2(h2_Alldphidz_cut, "All matched after PassCaloCut;dphi [rad];dz [cm]", 300, -3.5, 3.5, 400, -200, 200);

    MAKE2(h2_dphi_calo_phi, "dphi vs calo_phi;calo_phi [rad];dphi [rad]", 300, -3.5, 3.5, 300, -3.5, 3.5);
    MAKE2(h2_dphi_calo_phi_cut, "dphi vs calo_phi after PassCaloCut;calo_phi [rad];dphi [rad]", 300, -3.5, 3.5, 300, -3.5, 3.5);

    MAKE2(h2_dz_calo_z, "dz vs calo_z;calo_z [cm];dz [cm]", 400, -200, 200, 400, -200, 200);
    MAKE2(h2_dz_calo_z_cut, "dz vs calo_z after PassCaloCut;calo_z [cm];dz [cm]", 400, -200, 200, 400, -200, 200);

    MAKE2(h2_dphi_proj_phi, "dphi vs trk_proj_phi;trk_proj_phi [rad];dphi [rad]", 300, -3.5, 3.5, 300, -3.5, 3.5);
    MAKE2(h2_dphi_proj_phi_cut, "dphi vs trk_proj_phi after PassCaloCut;trk_proj_phi [rad];dphi [rad]", 300, -3.5, 3.5, 300, -3.5, 3.5);

    MAKE2(h2_dz_proj_z, "dz vs trk_proj_z;trk_proj_z [cm];dz [cm]", 400, -200, 200, 400, -200, 200);
    MAKE2(h2_dz_proj_z_cut, "dz vs trk_proj_z after PassCaloCut;trk_proj_z [cm];dz [cm]", 400, -200, 200, 400, -200, 200);

    MAKE2(h2_phi_trk_calo, "trk_proj_phi vs calo_phi;trk_proj_phi [rad];calo_phi [rad]", 300, -3.5, 3.5, 300, -3.5, 3.5);
    MAKE2(h2_phi_trk_calo_cut, "trk_proj_phi vs calo_phi after PassCaloCut;trk_proj_phi [rad];calo_phi [rad]", 300, -3.5, 3.5, 300, -3.5, 3.5);

    MAKE2(h2_z_trk_calo, "trk_proj_z vs calo_z;trk_proj_z [cm];calo_z [cm]", 400, -200, 200, 400, -200, 200);
    MAKE2(h2_z_trk_calo_cut, "trk_proj_z vs calo_z after PassCaloCut;trk_proj_z [cm];calo_z [cm]", 400, -200, 200, 400, -200, 200);

    // Charge/pT-split histograms. Each base histogram has a matching _cut histogram.
    const char* categories[] = {
        "charge_plus_highpt", "charge_minus_highpt",
        "charge_plus_lowpt",  "charge_minus_lowpt"
    };

    // Declare the category histograms with unique ROOT object names.
    #define MAKE_CAT_PAIR(var, title, nx, xmin, xmax, ny, ymin, ymax) \
        TH2F *var = new TH2F(#var, title, nx, xmin, xmax, ny, ymin, ymax); \
        TH2F *var##_cut = new TH2F(#var "_cut", title " after PassCaloCut", nx, xmin, xmax, ny, ymin, ymax)

    // trk projection phi vs calorimeter phi
    MAKE_CAT_PAIR(h2_phi_trk_calo_charge_plus_highpt,  "trk_proj_phi vs calo_phi, charge>0, pT>=3.5;trk_proj_phi [rad];calo_phi [rad]", 300, -3.5, 3.5, 300, -3.5, 3.5);
    MAKE_CAT_PAIR(h2_phi_trk_calo_charge_minus_highpt, "trk_proj_phi vs calo_phi, charge<0, pT>=3.5;trk_proj_phi [rad];calo_phi [rad]", 300, -3.5, 3.5, 300, -3.5, 3.5);
    MAKE_CAT_PAIR(h2_phi_trk_calo_charge_plus_lowpt,   "trk_proj_phi vs calo_phi, charge>0, pT<3.5;trk_proj_phi [rad];calo_phi [rad]", 300, -3.5, 3.5, 300, -3.5, 3.5);
    MAKE_CAT_PAIR(h2_phi_trk_calo_charge_minus_lowpt,  "trk_proj_phi vs calo_phi, charge<0, pT<3.5;trk_proj_phi [rad];calo_phi [rad]", 300, -3.5, 3.5, 300, -3.5, 3.5);

    // trk projection z vs calorimeter z
    MAKE_CAT_PAIR(h2_z_trk_calo_charge_plus_highpt,  "trk_proj_z vs calo_z, charge>0, pT>=3.5;trk_proj_z [cm];calo_z [cm]", 400, -200, 200, 400, -200, 200);
    MAKE_CAT_PAIR(h2_z_trk_calo_charge_minus_highpt, "trk_proj_z vs calo_z, charge<0, pT>=3.5;trk_proj_z [cm];calo_z [cm]", 400, -200, 200, 400, -200, 200);
    MAKE_CAT_PAIR(h2_z_trk_calo_charge_plus_lowpt,   "trk_proj_z vs calo_z, charge>0, pT<3.5;trk_proj_z [cm];calo_z [cm]", 400, -200, 200, 400, -200, 200);
    MAKE_CAT_PAIR(h2_z_trk_calo_charge_minus_lowpt,  "trk_proj_z vs calo_z, charge<0, pT<3.5;trk_proj_z [cm];calo_z [cm]", 400, -200, 200, 400, -200, 200);

    // dphi vs calorimeter phi
    MAKE_CAT_PAIR(h2_dphi_calo_phi_charge_plus_highpt,  "dphi vs calo_phi, charge>0, pT>=3.5;calo_phi [rad];dphi [rad]", 300, -3.5, 3.5, 300, -3.5, 3.5);
    MAKE_CAT_PAIR(h2_dphi_calo_phi_charge_minus_highpt, "dphi vs calo_phi, charge<0, pT>=3.5;calo_phi [rad];dphi [rad]", 300, -3.5, 3.5, 300, -3.5, 3.5);
    MAKE_CAT_PAIR(h2_dphi_calo_phi_charge_plus_lowpt,   "dphi vs calo_phi, charge>0, pT<3.5;calo_phi [rad];dphi [rad]", 300, -3.5, 3.5, 300, -3.5, 3.5);
    MAKE_CAT_PAIR(h2_dphi_calo_phi_charge_minus_lowpt,  "dphi vs calo_phi, charge<0, pT<3.5;calo_phi [rad];dphi [rad]", 300, -3.5, 3.5, 300, -3.5, 3.5);

    // dz vs calorimeter z
    MAKE_CAT_PAIR(h2_dz_calo_z_charge_plus_highpt,  "dz vs calo_z, charge>0, pT>=3.5;calo_z [cm];dz [cm]", 400, -200, 200, 400, -200, 200);
    MAKE_CAT_PAIR(h2_dz_calo_z_charge_minus_highpt, "dz vs calo_z, charge<0, pT>=3.5;calo_z [cm];dz [cm]", 400, -200, 200, 400, -200, 200);
    MAKE_CAT_PAIR(h2_dz_calo_z_charge_plus_lowpt,   "dz vs calo_z, charge>0, pT<3.5;calo_z [cm];dz [cm]", 400, -200, 200, 400, -200, 200);
    MAKE_CAT_PAIR(h2_dz_calo_z_charge_minus_lowpt,  "dz vs calo_z, charge<0, pT<3.5;calo_z [cm];dz [cm]", 400, -200, 200, 400, -200, 200);

    // dphi vs track projection phi
    MAKE_CAT_PAIR(h2_dphi_proj_phi_charge_plus_highpt,  "dphi vs trk_proj_phi, charge>0, pT>=3.5;trk_proj_phi [rad];dphi [rad]", 300, -3.5, 3.5, 300, -3.5, 3.5);
    MAKE_CAT_PAIR(h2_dphi_proj_phi_charge_minus_highpt, "dphi vs trk_proj_phi, charge<0, pT>=3.5;trk_proj_phi [rad];dphi [rad]", 300, -3.5, 3.5, 300, -3.5, 3.5);
    MAKE_CAT_PAIR(h2_dphi_proj_phi_charge_plus_lowpt,   "dphi vs trk_proj_phi, charge>0, pT<3.5;trk_proj_phi [rad];dphi [rad]", 300, -3.5, 3.5, 300, -3.5, 3.5);
    MAKE_CAT_PAIR(h2_dphi_proj_phi_charge_minus_lowpt,  "dphi vs trk_proj_phi, charge<0, pT<3.5;trk_proj_phi [rad];dphi [rad]", 300, -3.5, 3.5, 300, -3.5, 3.5);

    // dz vs track projection z
    MAKE_CAT_PAIR(h2_dz_proj_z_charge_plus_highpt,  "dz vs trk_proj_z, charge>0, pT>=3.5;trk_proj_z [cm];dz [cm]", 400, -200, 200, 400, -200, 200);
    MAKE_CAT_PAIR(h2_dz_proj_z_charge_minus_highpt, "dz vs trk_proj_z, charge<0, pT>=3.5;trk_proj_z [cm];dz [cm]", 400, -200, 200, 400, -200, 200);
    MAKE_CAT_PAIR(h2_dz_proj_z_charge_plus_lowpt,   "dz vs trk_proj_z, charge>0, pT<3.5;trk_proj_z [cm];dz [cm]", 400, -200, 200, 400, -200, 200);
    MAKE_CAT_PAIR(h2_dz_proj_z_charge_minus_lowpt,  "dz vs trk_proj_z, charge<0, pT<3.5;trk_proj_z [cm];dz [cm]", 400, -200, 200, 400, -200, 200);

    Long64_t nentries = tree->GetEntries();
    for (Long64_t i = 0; i < nentries; ++i) {
        tree->GetEntry(i);

        // Fill all non-_cut histograms with no calorimeter cut.
        h2_calohit->Fill(calo_phi, calo_z);
        h2_Alldphidz->Fill(dphi, dz);
        h2_dphi_calo_phi->Fill(calo_phi, dphi);
        h2_dz_calo_z->Fill(calo_z, dz);
        h2_dphi_proj_phi->Fill(trk_proj_phi, dphi);
        h2_dz_proj_z->Fill(trk_proj_z, dz);
        h2_phi_trk_calo->Fill(trk_proj_phi, calo_phi);
        h2_z_trk_calo->Fill(trk_proj_z, calo_z);

        const bool highpt = (trk_pt >= 3.5f);
        const bool lowpt  = (trk_pt < 3.5f);
        const bool plus   = (trk_charge > 0);
        const bool minus  = (trk_charge < 0);

        // Fill the charge/pT split base histograms before applying PassCaloCut.
        if (plus && highpt) {
            h2_phi_trk_calo_charge_plus_highpt->Fill(trk_proj_phi, calo_phi);
            h2_z_trk_calo_charge_plus_highpt->Fill(trk_proj_z, calo_z);
            h2_dphi_calo_phi_charge_plus_highpt->Fill(calo_phi, dphi);
            h2_dz_calo_z_charge_plus_highpt->Fill(calo_z, dz);
            h2_dphi_proj_phi_charge_plus_highpt->Fill(trk_proj_phi, dphi);
            h2_dz_proj_z_charge_plus_highpt->Fill(trk_proj_z, dz);
        } else if (minus && highpt) {
            h2_phi_trk_calo_charge_minus_highpt->Fill(trk_proj_phi, calo_phi);
            h2_z_trk_calo_charge_minus_highpt->Fill(trk_proj_z, calo_z);
            h2_dphi_calo_phi_charge_minus_highpt->Fill(calo_phi, dphi);
            h2_dz_calo_z_charge_minus_highpt->Fill(calo_z, dz);
            h2_dphi_proj_phi_charge_minus_highpt->Fill(trk_proj_phi, dphi);
            h2_dz_proj_z_charge_minus_highpt->Fill(trk_proj_z, dz);
        } else if (plus && lowpt) {
            h2_phi_trk_calo_charge_plus_lowpt->Fill(trk_proj_phi, calo_phi);
            h2_z_trk_calo_charge_plus_lowpt->Fill(trk_proj_z, calo_z);
            h2_dphi_calo_phi_charge_plus_lowpt->Fill(calo_phi, dphi);
            h2_dz_calo_z_charge_plus_lowpt->Fill(calo_z, dz);
            h2_dphi_proj_phi_charge_plus_lowpt->Fill(trk_proj_phi, dphi);
            h2_dz_proj_z_charge_plus_lowpt->Fill(trk_proj_z, dz);
        } else if (minus && lowpt) {
            h2_phi_trk_calo_charge_minus_lowpt->Fill(trk_proj_phi, calo_phi);
            h2_z_trk_calo_charge_minus_lowpt->Fill(trk_proj_z, calo_z);
            h2_dphi_calo_phi_charge_minus_lowpt->Fill(calo_phi, dphi);
            h2_dz_calo_z_charge_minus_lowpt->Fill(calo_z, dz);
            h2_dphi_proj_phi_charge_minus_lowpt->Fill(trk_proj_phi, dphi);
            h2_dz_proj_z_charge_minus_lowpt->Fill(trk_proj_z, dz);
        }

        // From this point onward, only entries passing PassCaloCut are filled.
        if (!PassCaloCut(calo_phi, calo_z)) continue;
        bool bool_dphi = PassDphiCut(dphi);
        bool bool_dz = PassDzCut(dz);

        h2_calohit_cut->Fill(calo_phi, calo_z);
        h2_Alldphidz_cut->Fill(dphi, dz);
        if(bool_dphi)
        {
            h2_dz_calo_z_cut->Fill(calo_z, dz);
            h2_dz_proj_z_cut->Fill(trk_proj_z, dz);
            h2_z_trk_calo_cut->Fill(trk_proj_z, calo_z);
        }
        if(bool_dz)
        {
            h2_dphi_calo_phi_cut->Fill(calo_phi, dphi);
            h2_dphi_proj_phi_cut->Fill(trk_proj_phi, dphi);
            h2_phi_trk_calo_cut->Fill(trk_proj_phi, calo_phi);
        }
        

        if (plus && highpt) {
            if(bool_dphi)
            {
                h2_z_trk_calo_charge_plus_highpt_cut->Fill(trk_proj_z, calo_z);
                h2_dz_proj_z_charge_plus_highpt_cut->Fill(trk_proj_z, dz);
                h2_dz_calo_z_charge_plus_highpt_cut->Fill(calo_z, dz);
            }
            if(bool_dz)
            {
                h2_phi_trk_calo_charge_plus_highpt_cut->Fill(trk_proj_phi, calo_phi);
                h2_dphi_calo_phi_charge_plus_highpt_cut->Fill(calo_phi, dphi);
                h2_dphi_proj_phi_charge_plus_highpt_cut->Fill(trk_proj_phi, dphi);
            }            
        } else if (minus && highpt) {
            if(bool_dphi)
            {
                h2_z_trk_calo_charge_minus_highpt_cut->Fill(trk_proj_z, calo_z);
                h2_dz_calo_z_charge_minus_highpt_cut->Fill(calo_z, dz);
                h2_dz_proj_z_charge_minus_highpt_cut->Fill(trk_proj_z, dz);
            }
            if(bool_dz)
            {
                h2_phi_trk_calo_charge_minus_highpt_cut->Fill(trk_proj_phi, calo_phi);            
                h2_dphi_calo_phi_charge_minus_highpt_cut->Fill(calo_phi, dphi);            
                h2_dphi_proj_phi_charge_minus_highpt_cut->Fill(trk_proj_phi, dphi);
            }
            
        } else if (plus && lowpt) {
            if(bool_dphi){
                h2_z_trk_calo_charge_plus_lowpt_cut->Fill(trk_proj_z, calo_z);
                h2_dz_calo_z_charge_plus_lowpt_cut->Fill(calo_z, dz);
                h2_dz_proj_z_charge_plus_lowpt_cut->Fill(trk_proj_z, dz);

            }
            if(bool_dz){
                h2_phi_trk_calo_charge_plus_lowpt_cut->Fill(trk_proj_phi, calo_phi);
                h2_dphi_calo_phi_charge_plus_lowpt_cut->Fill(calo_phi, dphi);
                h2_dphi_proj_phi_charge_plus_lowpt_cut->Fill(trk_proj_phi, dphi);
            }
        } else if (minus && lowpt) {
            if(bool_dphi){
                h2_z_trk_calo_charge_minus_lowpt_cut->Fill(trk_proj_z, calo_z);
                h2_dz_calo_z_charge_minus_lowpt_cut->Fill(calo_z, dz);
                h2_dz_proj_z_charge_minus_lowpt_cut->Fill(trk_proj_z, dz);
            }
            if(bool_dz){
                h2_phi_trk_calo_charge_minus_lowpt_cut->Fill(trk_proj_phi, calo_phi);
                h2_dphi_calo_phi_charge_minus_lowpt_cut->Fill(calo_phi, dphi);
                h2_dphi_proj_phi_charge_minus_lowpt_cut->Fill(trk_proj_phi, dphi);
            }
        }
    }

    fout->cd();

    // Write all histogram pairs.
    h2_calohit->Write(); h2_calohit_cut->Write();
    h2_Alldphidz->Write(); h2_Alldphidz_cut->Write();
    h2_dphi_calo_phi->Write(); h2_dphi_calo_phi_cut->Write();
    h2_dz_calo_z->Write(); h2_dz_calo_z_cut->Write();
    h2_dphi_proj_phi->Write(); h2_dphi_proj_phi_cut->Write();
    h2_dz_proj_z->Write(); h2_dz_proj_z_cut->Write();
    h2_phi_trk_calo->Write(); h2_phi_trk_calo_cut->Write();
    h2_z_trk_calo->Write(); h2_z_trk_calo_cut->Write();

    #define WRITE_CAT_PAIR(var) var->Write(); var##_cut->Write()
    WRITE_CAT_PAIR(h2_phi_trk_calo_charge_plus_highpt);
    WRITE_CAT_PAIR(h2_phi_trk_calo_charge_minus_highpt);
    WRITE_CAT_PAIR(h2_phi_trk_calo_charge_plus_lowpt);
    WRITE_CAT_PAIR(h2_phi_trk_calo_charge_minus_lowpt);

    WRITE_CAT_PAIR(h2_z_trk_calo_charge_plus_highpt);
    WRITE_CAT_PAIR(h2_z_trk_calo_charge_minus_highpt);
    WRITE_CAT_PAIR(h2_z_trk_calo_charge_plus_lowpt);
    WRITE_CAT_PAIR(h2_z_trk_calo_charge_minus_lowpt);

    WRITE_CAT_PAIR(h2_dphi_calo_phi_charge_plus_highpt);
    WRITE_CAT_PAIR(h2_dphi_calo_phi_charge_minus_highpt);
    WRITE_CAT_PAIR(h2_dphi_calo_phi_charge_plus_lowpt);
    WRITE_CAT_PAIR(h2_dphi_calo_phi_charge_minus_lowpt);

    WRITE_CAT_PAIR(h2_dz_calo_z_charge_plus_highpt);
    WRITE_CAT_PAIR(h2_dz_calo_z_charge_minus_highpt);
    WRITE_CAT_PAIR(h2_dz_calo_z_charge_plus_lowpt);
    WRITE_CAT_PAIR(h2_dz_calo_z_charge_minus_lowpt);

    WRITE_CAT_PAIR(h2_dphi_proj_phi_charge_plus_highpt);
    WRITE_CAT_PAIR(h2_dphi_proj_phi_charge_minus_highpt);
    WRITE_CAT_PAIR(h2_dphi_proj_phi_charge_plus_lowpt);
    WRITE_CAT_PAIR(h2_dphi_proj_phi_charge_minus_lowpt);

    WRITE_CAT_PAIR(h2_dz_proj_z_charge_plus_highpt);
    WRITE_CAT_PAIR(h2_dz_proj_z_charge_minus_highpt);
    WRITE_CAT_PAIR(h2_dz_proj_z_charge_plus_lowpt);
    WRITE_CAT_PAIR(h2_dz_proj_z_charge_minus_lowpt);

    fout->Close();
    f->Close();
    delete fout;
    delete f;

    #undef WRITE_CAT_PAIR
    #undef MAKE_CAT_PAIR
    #undef MAKE2
}
