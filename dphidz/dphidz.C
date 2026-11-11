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

void dphidz(const string inputfile = "/sphenix/u/mizukishigihara/tgmizuki/Jpsi_renew/anamacro/result/data_merged/53018_allmatched_merged.root")
{
    TFile *f = TFile::Open(inputfile.c_str(),"READ");
    if((!f ||f->IsZombie())){
        std::cerr <<"Error: cannot open file: "<< inputfile.c_str() << std::endl;
        return;
    }

    //outputfile
    string outdir = "/sphenix/u/mizukishigihara/tgmizuki/Jpsi_renew/dphidz/result/data";
    TString fname_out = gSystem->BaseName(inputfile.c_str());
    fname_out.ReplaceAll("allmatched_merged", "Alldphidz");
    TFile *fout0 = new TFile(Form("%s/200_1000%s", outdir.c_str(), fname_out.Data()),"RECREATE");
    TFile *fout1 = new TFile(Form("%s/1000_2000%s", outdir.c_str(), fname_out.Data()),"RECREATE");
    TFile *fout2 = new TFile(Form("%s/2000_3000%s", outdir.c_str(), fname_out.Data()),"RECREATE");
    TFile *fout3 = new TFile(Form("%s/3000_4000%s", outdir.c_str(), fname_out.Data()),"RECREATE");
    TFile *fout4 = new TFile(Form("%s/4000_5000%s", outdir.c_str(), fname_out.Data()),"RECREATE");


    TTree *AllMatchedSiEMC = (TTree*)f->Get("AllMatchedSiEMC");

    float trk_pt;
    int trk_charge;
    float dphi;
    float dz;
    float calo_phi;
    float calo_z;

    AllMatchedSiEMC->SetBranchAddress("trk_pt", &trk_pt);
    AllMatchedSiEMC->SetBranchAddress("dphi", &dphi);
    AllMatchedSiEMC->SetBranchAddress("dz", &dz);
    AllMatchedSiEMC->SetBranchAddress("trk_charge", &trk_charge);
    AllMatchedSiEMC->SetBranchAddress("calo_phi", &calo_phi);
    AllMatchedSiEMC->SetBranchAddress("calo_z", &calo_z);

    TH2F *h2_Alldphidz = new TH2F("h2_Alldphipz","All_dphi_dz; dphi[rad]; dz[cm]",300,-3.5,3.5,400,-400,400);
    TH1F *h_Alldphi = new TH1F("h_Alldphi","Alldphi;dphi [rad];Counts",300,-3.5,3.5);
    TH1F *h_Alldz = new TH1F("h_Alldz","Alldz;dz [cm];Counts",400,-400,400);
    TH1F *h_Alldphi_charge_plas = new TH1F("h_Alldphi_charge_plas","Alldphi_charge_plas;dphi [rad];Counts",300,-3.5,3.5);
    TH1F *h_Alldz_charge_plas = new TH1F("h_Alldz_charge_plas","Alldz_charge_plas;dz [cm];Counts",400,-400,400);
    TH1F *h_Alldphi_charge_minus = new TH1F("h_Alldphi_charge_minus","Alldphi_charge_minus;dphi [rad];Counts",300,-3.5,3.5);
    TH1F *h_Alldz_charge_minus = new TH1F("h_Alldz_charge_minus","Alldz_charge_minus;dz [cm];Counts",400,-400,400);
    //200_1000
    TH1F *h_Alldphi_200_400MeV = new TH1F("h_Alldphi_200_400MeV","Alldphi_200_400MeV;dphi [rad];Counts",300,-3.5,3.5);
    TH1F *h_Alldz_200_400MeV = new TH1F("h_Alldz_200_400MeV","Alldz_200_400MeV;dz [cm];Counts",400,-400,400);
    TH1F *h_Alldphi_400_600MeV = new TH1F("h_Alldphi_400_600MeV","Alldphi_400_600MeV;dphi [rad];Counts",300,-3.5,3.5);
    TH1F *h_Alldz_400_600MeV = new TH1F("h_Alldz_400_600MeV","Alldz_400_600MeV;dz [cm];Counts",400,-400,400);
    TH1F *h_Alldphi_600_800MeV = new TH1F("h_Alldphi_600_800MeV","Alldphi_600_800MeV;dphi [rad];Counts",300,-3.5,3.5);
    TH1F *h_Alldz_600_800MeV = new TH1F("h_Alldz_600_800MeV","Alldz_600_800MeV;dz [cm];Counts",400,-400,400);
    TH1F *h_Alldphi_800_1000MeV = new TH1F("h_Alldphi_800_1000MeV","Alldphi_800_1000MeV;dphi [rad];Counts",300,-3.5,3.5);
    TH1F *h_Alldz_800_1000MeV = new TH1F("h_Alldz_800_1000MeV","Alldz_800_1000MeV;dz [cm];Counts",400,-400,400);
    //1000_2000MeV
    TH1F *h_Alldphi_1000_1200MeV = new TH1F("h_Alldphi_1000_1200MeV","Alldphi_1000_1200MeV;dphi [rad];Counts",300,-3.5,3.5);
    TH1F *h_Alldz_1000_1200MeV = new TH1F("h_Alldz_1000_1200MeV","Alldz_1000_1200MeV;dz [cm];Counts",400,-400,400);
    TH1F *h_Alldphi_1200_1400MeV = new TH1F("h_Alldphi_1200_1400MeV","Alldphi_1200_1400MeV;dphi [rad];Counts",300,-3.5,3.5);
    TH1F *h_Alldz_1200_1400MeV = new TH1F("h_Alldz_1200_1400MeV","Alldz_1200_1400MeV;dz [cm];Counts",400,-400,400);
    TH1F *h_Alldphi_1400_1600MeV = new TH1F("h_Alldphi_1400_1600MeV","Alldphi_1400_1600MeV;dphi [rad];Counts",300,-3.5,3.5);
    TH1F *h_Alldz_1400_1600MeV = new TH1F("h_Alldz_1400_1600MeV","Alldz_1400_1600MeV;dz [cm];Counts",400,-400,400);
    TH1F *h_Alldphi_1600_1800MeV = new TH1F("h_Alldphi_1600_1800MeV","Alldphi_1600_1800MeV;dphi [rad];Counts",300,-3.5,3.5);
    TH1F *h_Alldz_1600_1800MeV = new TH1F("h_Alldz_1600_1800MeV","Alldz_1600_1800MeV;dz [cm];Counts",400,-400,400);
    TH1F *h_Alldphi_1800_2000MeV = new TH1F("h_Alldphi_1800_2000MeV","Alldphi_1800_2000MeV;dphi [rad];Counts",300,-3.5,3.5);
    TH1F *h_Alldz_1800_2000MeV = new TH1F("h_Alldz_1800_2000MeV","Alldz_1800_2000MeV;dz [cm];Counts",400,-400,400);
    //2000_3000
    TH1F *h_Alldphi_2000_2200MeV = new TH1F("h_Alldphi_2000_2200MeV","Alldphi_2000_2200MeV;dphi [rad];Counts",300,-3.5,3.5);
    TH1F *h_Alldz_2000_2200MeV = new TH1F("h_Alldz_2000_2200MeV","Alldz_2000_2200MeV;dz [cm];Counts",400,-400,400);
    TH1F *h_Alldphi_2200_2400MeV = new TH1F("h_Alldphi_2200_2400MeV","Alldphi_2200_2400MeV;dphi [rad];Counts",300,-3.5,3.5);
    TH1F *h_Alldz_2200_2400MeV = new TH1F("h_Alldz_2200_2400MeV","Alldz_2200_2400MeV;dz [cm];Counts",400,-400,400);
    TH1F *h_Alldphi_2400_2600MeV = new TH1F("h_Alldphi_2400_2600MeV","Alldphi_2400_2600MeV;dphi [rad];Counts",300,-3.5,3.5);
    TH1F *h_Alldz_2400_2600MeV = new TH1F("h_Alldz_2400_2600MeV","Alldz_2400_2600MeV;dz [cm];Counts",400,-400,400);
    TH1F *h_Alldphi_2600_2800MeV = new TH1F("h_Alldphi_2600_2800MeV","Alldphi_2600_2800MeV;dphi [rad];Counts",300,-3.5,3.5);
    TH1F *h_Alldz_2600_2800MeV = new TH1F("h_Alldz_2600_2800MeV","Alldz_2600_2800MeV;dz [cm];Counts",400,-400,400);
    TH1F *h_Alldphi_2800_3000MeV = new TH1F("h_Alldphi_2800_3000MeV","Alldphi_2800_3000MeV;dphi [rad];Counts",300,-3.5,3.5);
    TH1F *h_Alldz_2800_3000MeV = new TH1F("h_Alldz_2800_3000MeV","Alldz_2800_3000MeV;dz [cm];Counts",400,-400,400);
    //3000_4000
    TH1F *h_Alldphi_3000_3200MeV = new TH1F("h_Alldphi_3000_3200MeV","Alldphi_3000_3200MeV;dphi [rad];Counts",300,-3.5,3.5);
    TH1F *h_Alldz_3000_3200MeV = new TH1F("h_Alldz_3000_3200MeV","Alldz_3000_3200MeV;dz [cm];Counts",400,-400,400);
    TH1F *h_Alldphi_3200_3400MeV = new TH1F("h_Alldphi_3200_3400MeV","Alldphi_3200_3400MeV;dphi [rad];Counts",300,-3.5,3.5);
    TH1F *h_Alldz_3200_3400MeV = new TH1F("h_Alldz_3200_3400MeV","Alldz_3200_3400MeV;dz [cm];Counts",400,-400,400);
    TH1F *h_Alldphi_3400_3600MeV = new TH1F("h_Alldphi_3400_3600MeV","Alldphi_3400_3600MeV;dphi [rad];Counts",300,-3.5,3.5);
    TH1F *h_Alldz_3400_3600MeV = new TH1F("h_Alldz_3400_3600MeV","Alldz_3400_3600MeV;dz [cm];Counts",400,-400,400);
    TH1F *h_Alldphi_3600_3800MeV = new TH1F("h_Alldphi_3600_3800MeV","Alldphi_3600_3800MeV;dphi [rad];Counts",300,-3.5,3.5);
    TH1F *h_Alldz_3600_3800MeV = new TH1F("h_Alldz_3600_3800MeV","Alldz_3600_3800MeV;dz [cm];Counts",400,-400,400);
    TH1F *h_Alldphi_3800_4000MeV = new TH1F("h_Alldphi_3800_4000MeV","Alldphi_3800_4000MeV;dphi [rad];Counts",300,-3.5,3.5);
    TH1F *h_Alldz_3800_4000MeV = new TH1F("h_Alldz_3800_4000MeV","Alldz_3800_4000MeV;dz [cm];Counts",400,-400,400);
    //4000_5000
    TH1F *h_Alldphi_4000_4200MeV = new TH1F("h_Alldphi_4000_4200MeV","Alldphi_4000_4200MeV;dphi [rad];Counts",300,-3.5,3.5);
    TH1F *h_Alldz_4000_4200MeV = new TH1F("h_Alldz_4000_4200MeV","Alldz_4000_4200MeV;dz [cm];Counts",400,-400,400);
    TH1F *h_Alldphi_4200_4400MeV = new TH1F("h_Alldphi_4200_4400MeV","Alldphi_4200_4400MeV;dphi [rad];Counts",300,-3.5,3.5);
    TH1F *h_Alldz_4200_4400MeV = new TH1F("h_Alldz_4200_4400MeV","Alldz_4200_4400MeV;dz [cm];Counts",400,-400,400);
    TH1F *h_Alldphi_4400_4600MeV = new TH1F("h_Alldphi_4400_4600MeV","Alldphi_4400_4600MeV;dphi [rad];Counts",300,-3.5,3.5);
    TH1F *h_Alldz_4400_4600MeV = new TH1F("h_Alldz_4400_4600MeV","Alldz_4400_4600MeV;dz [cm];Counts",400,-400,400);
    TH1F *h_Alldphi_4600_4800MeV = new TH1F("h_Alldphi_4600_4800MeV","Alldphi_4600_4800MeV;dphi [rad];Counts",300,-3.5,3.5);
    TH1F *h_Alldz_4600_4800MeV = new TH1F("h_Alldz_4600_4800MeV","Alldz_4600_4800MeV;dz [cm];Counts",400,-400,400);
    TH1F *h_Alldphi_4800_5000MeV = new TH1F("h_Alldphi_4800_5000MeV","Alldphi_4800_5000MeV;dphi [rad];Counts",300,-3.5,3.5);
    TH1F *h_Alldz_4800_5000MeV = new TH1F("h_Alldz_4800_5000MeV","Alldz_4800_5000MeV;dz [cm];Counts",400,-400,400);

    
    //dphi
    Long64_t n = AllMatchedSiEMC->GetEntries();
    for (Long64_t i = 0; i < n; i++) 
    {
        AllMatchedSiEMC->GetEntry(i);

        if(!PassCaloCut(calo_phi, calo_z))continue;
        h2_Alldphidz->Fill(dphi,dz);
        if(!PassDzCut(dz))continue;
 
        h_Alldphi->Fill(dphi);

        if(trk_charge < 0)
        {
            h_Alldphi_charge_minus->Fill(dphi);
        }
        if(trk_charge > 0)
        {
            h_Alldphi_charge_plas->Fill(dphi);
        }

        if(0.2 <= trk_pt && trk_pt< 0.4)
        {
            h_Alldphi_200_400MeV->Fill(dphi);
        }
        if(0.4 <= trk_pt && trk_pt< 0.6)
        {
            h_Alldphi_400_600MeV->Fill(dphi);
        }
         if(0.6 <= trk_pt && trk_pt< 0.8)
        {
            h_Alldphi_600_800MeV->Fill(dphi);
        }
         if(0.8 <= trk_pt && trk_pt< 1.0)
        {
            h_Alldphi_800_1000MeV->Fill(dphi);
        }
        if(1.0 <= trk_pt && trk_pt< 1.2)
        {
            h_Alldphi_1000_1200MeV->Fill(dphi);
        }
        if(1.2 <= trk_pt && trk_pt< 1.4)
        {
            h_Alldphi_1200_1400MeV->Fill(dphi);
        }
        if(1.4 <= trk_pt && trk_pt< 1.6)
        {
            h_Alldphi_1400_1600MeV->Fill(dphi);
        }
        if(1.6 <= trk_pt && trk_pt< 1.8)
        {
            h_Alldphi_1600_1800MeV->Fill(dphi);
        }
        if(1.8 <= trk_pt && trk_pt< 2.0)
        {
            h_Alldphi_1800_2000MeV->Fill(dphi);
        }
        if(2.0 <= trk_pt && trk_pt< 2.2)
        {
            h_Alldphi_2000_2200MeV->Fill(dphi);
        }
        if(2.2 <= trk_pt && trk_pt< 2.4)
        {
            h_Alldphi_2200_2400MeV->Fill(dphi);
        }
        if(2.4 <= trk_pt && trk_pt< 2.6)
        {
            h_Alldphi_2400_2600MeV->Fill(dphi);
        }
        if(2.6 <= trk_pt && trk_pt< 2.8)
        {
            h_Alldphi_2600_2800MeV->Fill(dphi);
        }
        if(2.8 <= trk_pt && trk_pt< 3.0)
        {
            h_Alldphi_2800_3000MeV->Fill(dphi);
        }
        if(3.0 <= trk_pt && trk_pt< 3.2)
        {
            h_Alldphi_3000_3200MeV->Fill(dphi);
        }
        if(3.2 <= trk_pt && trk_pt< 3.4)
        {
            h_Alldphi_3200_3400MeV->Fill(dphi);
        }
        if(3.4 <= trk_pt && trk_pt< 3.6)
        {
            h_Alldphi_3400_3600MeV->Fill(dphi);
        }
        if(3.6 <= trk_pt && trk_pt< 3.8)
        {
            h_Alldphi_3600_3800MeV->Fill(dphi);
        }
        if(3.8 <= trk_pt && trk_pt< 4.0)
        {
            h_Alldphi_3800_4000MeV->Fill(dphi);
        }
        if(4.0 <= trk_pt && trk_pt< 4.2)
        {
            h_Alldphi_4000_4200MeV->Fill(dphi);
        }
        if(4.2 <= trk_pt && trk_pt< 4.4)
        {
            h_Alldphi_4200_4400MeV->Fill(dphi);
        }
        if(4.4 <= trk_pt && trk_pt< 4.6)
        {
            h_Alldphi_4400_4600MeV->Fill(dphi);
        }
        if(4.6 <= trk_pt && trk_pt< 4.8)
        {
            h_Alldphi_4600_4800MeV->Fill(dphi);
        }
        if(4.8 <= trk_pt && trk_pt< 5.0)
        {
            h_Alldphi_4800_5000MeV->Fill(dphi);
        }

    }

    //dz
    for (Long64_t i = 0; i < n; i++) 
    {
        AllMatchedSiEMC->GetEntry(i);

        if(!PassCaloCut(calo_phi, calo_z))continue;
        if(!PassDphiCut(dphi))continue;

        h_Alldz->Fill(dz);

        if(trk_charge < 0)
        {
            h_Alldz_charge_minus->Fill(dz);
        }
        if(trk_charge > 0)
        {
            h_Alldz_charge_plas->Fill(dz);
        }

        if(0.2 <= trk_pt && trk_pt< 0.4)
        {
            h_Alldz_200_400MeV->Fill(dz);
        }
        if(0.4 <= trk_pt && trk_pt< 0.6)
        {
            h_Alldz_400_600MeV->Fill(dz);
        }
         if(0.6 <= trk_pt && trk_pt< 0.8)
        {
            h_Alldz_600_800MeV->Fill(dz);
        }
         if(0.8 <= trk_pt && trk_pt< 1.0)
        {
            h_Alldz_800_1000MeV->Fill(dz);
        }
        if(1.0 <= trk_pt && trk_pt< 1.2)
        {
            h_Alldz_1000_1200MeV->Fill(dz);
        }
        if(1.2 <= trk_pt && trk_pt< 1.4)
        {
            h_Alldz_1200_1400MeV->Fill(dz);
        }
        if(1.4 <= trk_pt && trk_pt< 1.6)
        {
            h_Alldz_1400_1600MeV->Fill(dz);
        }
        if(1.6 <= trk_pt && trk_pt< 1.8)
        {
            h_Alldz_1600_1800MeV->Fill(dz);
        }
        if(1.8 <= trk_pt && trk_pt< 2.0)
        {
            h_Alldz_1800_2000MeV->Fill(dz);
        }
        if(2.0 <= trk_pt && trk_pt< 2.2)
        {
            h_Alldz_2000_2200MeV->Fill(dz);
        }
        if(2.2 <= trk_pt && trk_pt< 2.4)
        {
            h_Alldz_2200_2400MeV->Fill(dz);
        }
        if(2.4 <= trk_pt && trk_pt< 2.6)
        {
            h_Alldz_2400_2600MeV->Fill(dz);
        }
        if(2.6 <= trk_pt && trk_pt< 2.8)
        {
            h_Alldz_2600_2800MeV->Fill(dz);
        }
        if(2.8 <= trk_pt && trk_pt< 3.0)
        {
            h_Alldz_2800_3000MeV->Fill(dz);
        }
        if(3.0 <= trk_pt && trk_pt< 3.2)
        {
            h_Alldz_3000_3200MeV->Fill(dz);
        }
        if(3.2 <= trk_pt && trk_pt< 3.4)
        {
            h_Alldz_3200_3400MeV->Fill(dz);
        }
        if(3.4 <= trk_pt && trk_pt< 3.6)
        {
            h_Alldz_3400_3600MeV->Fill(dz);
        }
        if(3.6 <= trk_pt && trk_pt< 3.8)
        {
            h_Alldz_3600_3800MeV->Fill(dz);
        }
        if(3.8 <= trk_pt && trk_pt< 4.0)
        {
            h_Alldz_3800_4000MeV->Fill(dz);
        }
        if(4.0 <= trk_pt && trk_pt< 4.2)
        {
            h_Alldz_4000_4200MeV->Fill(dz);
        }
        if(4.2 <= trk_pt && trk_pt< 4.4)
        {
            h_Alldz_4200_4400MeV->Fill(dz);
        }
        if(4.4 <= trk_pt && trk_pt< 4.6)
        {
            h_Alldz_4400_4600MeV->Fill(dz);
        }
        if(4.6 <= trk_pt && trk_pt< 4.8)
        {
            h_Alldz_4600_4800MeV->Fill(dz);
        }
        if(4.8 <= trk_pt && trk_pt< 5.0)
        {
            h_Alldz_4800_5000MeV->Fill(dz);
        }

    }

    fout0->cd();
    h2_Alldphidz->Write();
    h_Alldphi->Write();
    h_Alldz->Write();
    h_Alldphi_charge_minus->Write();
    h_Alldz_charge_minus->Write();
    h_Alldphi_charge_plas->Write();
    h_Alldz_charge_plas->Write();

    h_Alldphi_200_400MeV->Write();
    h_Alldz_200_400MeV->Write();
    h_Alldphi_400_600MeV->Write();
    h_Alldz_400_600MeV->Write();
    h_Alldphi_600_800MeV->Write();
    h_Alldz_600_800MeV->Write();
    h_Alldphi_800_1000MeV->Write();
    h_Alldz_800_1000MeV->Write();

    fout1->cd();
    h2_Alldphidz->Write();
    h_Alldphi->Write();
    h_Alldz->Write();
    h_Alldphi_charge_minus->Write();
    h_Alldz_charge_minus->Write();
    h_Alldphi_charge_plas->Write();
    h_Alldz_charge_plas->Write();

    h_Alldphi_1000_1200MeV->Write();
    h_Alldz_1000_1200MeV->Write();
    h_Alldphi_1200_1400MeV->Write();
    h_Alldz_1200_1400MeV->Write();
    h_Alldphi_1400_1600MeV->Write();
    h_Alldz_1400_1600MeV->Write();
    h_Alldphi_1600_1800MeV->Write();
    h_Alldz_1600_1800MeV->Write();
    h_Alldphi_1800_2000MeV->Write();
    h_Alldz_1800_2000MeV->Write();

    fout2->cd();
    h2_Alldphidz->Write();
    h_Alldphi->Write();
    h_Alldz->Write();
    h_Alldphi_charge_minus->Write();
    h_Alldz_charge_minus->Write();
    h_Alldphi_charge_plas->Write();
    h_Alldz_charge_plas->Write();

    h_Alldphi_2000_2200MeV->Write();
    h_Alldz_2000_2200MeV->Write();
    h_Alldphi_2200_2400MeV->Write();
    h_Alldz_2200_2400MeV->Write();
    h_Alldphi_2400_2600MeV->Write();
    h_Alldz_2400_2600MeV->Write();
    h_Alldphi_2600_2800MeV->Write();
    h_Alldz_2600_2800MeV->Write();
    h_Alldphi_2800_3000MeV->Write();
    h_Alldz_2800_3000MeV->Write();

    fout3->cd();
    h2_Alldphidz->Write();
    h_Alldphi->Write();
    h_Alldz->Write();
    h_Alldphi_charge_minus->Write();
    h_Alldz_charge_minus->Write();
    h_Alldphi_charge_plas->Write();
    h_Alldz_charge_plas->Write();

    h_Alldphi_3000_3200MeV->Write();
    h_Alldz_3000_3200MeV->Write();
    h_Alldphi_3200_3400MeV->Write();
    h_Alldz_3200_3400MeV->Write();
    h_Alldphi_3400_3600MeV->Write();
    h_Alldz_3400_3600MeV->Write();
    h_Alldphi_3600_3800MeV->Write();
    h_Alldz_3600_3800MeV->Write();
    h_Alldphi_3800_4000MeV->Write();
    h_Alldz_3800_4000MeV->Write();

    fout4->cd();
    h2_Alldphidz->Write();
    h_Alldphi->Write();
    h_Alldz->Write();
    h_Alldphi_charge_minus->Write();
    h_Alldz_charge_minus->Write();
    h_Alldphi_charge_plas->Write();
    h_Alldz_charge_plas->Write();

    h_Alldphi_4000_4200MeV->Write();
    h_Alldz_4000_4200MeV->Write();
    h_Alldphi_4200_4400MeV->Write();
    h_Alldz_4200_4400MeV->Write();
    h_Alldphi_4400_4600MeV->Write();
    h_Alldz_4400_4600MeV->Write();
    h_Alldphi_4600_4800MeV->Write();
    h_Alldz_4600_4800MeV->Write();
    h_Alldphi_4800_5000MeV->Write();
    h_Alldz_4800_5000MeV->Write();

    fout0->Close();
    fout1->Close();
    fout2->Close();
    fout3->Close();
    fout4->Close();
    f->Close();


}