void titleStyle(TH1* h1);
void plot2hist(TH1* h1, TString leg1, TH1* h2, TString leg2, Bool_t isTitleStyle, TString savename="", Float_t xlegmin=0.4, Float_t ylegmin=0.4 );
void plot3hist(TH1* h1, TString leg1, TH1* h2, TString leg2, TH1* h3, TString leg3, Bool_t isTitleStyle,TString savename="", Float_t xlegmin=0.4, Float_t ylegmin=0.4 );
void test_th23th13(){
    //T2K 2022 data can be downloaded from https://zenodo.org/record/7741399
    TFile *f1 = new TFile("/Users/smithjainy/Documents/GitHub/nguyenthihien.github.io/Internship/T2K_arxiv2303.03222_DataRelease/Bayesian_DataRelease.root");

    TString output = f1->GetName();

        output.ReplaceAll(".root", "");


        // hsinsqth13_pred: Histogram for storing the distribution of predicted values for sinsqth13.
    // hsinsqth13_worc: Histogram for storing the distribution of random values for sinsqth13 with reactor.
    //hsinsqth13_wrc: Histogram for storing the distribution of random values for sinsqth13 without reactor.

    const char *phistname[] = {
    "h1D_th23posterior_woRC_both",//0
    "h1D_dm32posterior_woRC_both",//1
    "h1D_dCPposterior_woRC_both",//2
    "h1D_th13posterior_woRC_both",//3
    "h1D_th23posterior_woRC_IH",//4
    "h1D_dm32posterior_woRC_IH",//5
    "h1D_dCPposterior_woRC_IH",//6
    "h1D_th13posterior_woRC_IH",//7
    "h1D_th23posterior_woRC_NH",//8
    "h1D_dm32posterior_woRC_NH",//9
    "h1D_dCPposterior_woRC_NH",//10
    "h1D_th13posterior_woRC_NH",//11
    "h1D_th23posterior_wRC_both",//12
    "h1D_dm32posterior_wRC_both",//13
    "h1D_dCPposterior_wRC_both",//14
    "h1D_th13posterior_wRC_both",//15
    "h1D_th23posterior_wRC_IH",//16
    "h1D_dm32posterior_wRC_IH",//17
    "h1D_dCPposterior_wRC_IH",//18
    "h1D_th13posterior_wRC_IH",//19
    "h1D_th23posterior_wRC_NH",//20
    "h1D_dm32posterior_wRC_NH",//21
    "h1D_dCPposterior_wRC_NH",//22
    "h1D_th13posterior_wRC_NH"//23
       };
       Int_t NHIST1D =sizeof(phistname)/sizeof(phistname[0]);
       TH1D **hhist1d = new TH1D*[NHIST1D];
    
    for (Int_t ihist=0; ihist<NHIST1D; ++ihist) {
        hhist1d[ihist]=(TH1D*)f1->Get(phistname[ihist]);
    }
    
    /*new TCanvas;
    hhist1d[8]->SetLineColor(2);
    hhist1d[8]->Draw("hist");
    hhist1d[20]->SetLineColor(4);
    hhist1d[20]->Draw("histsame");
    gPad->Print(Form("test_%s_th23_nh_comp.eps",output.Data()));*/
    

    //No. of throw
    //Long_t NTHROW =hhist1d[8]->Integral();
    //cout<<"No. of throw "<<NTHROW<<endl;
    Int_t NTHROW = 1000000;
    double rnd_sinsqth23_worc;
    double rnd_sinsq2th23_worc;
    double rnd_sinsq2th23_mod_worc;
    double rnd_sinsqth13_worc;
    
    double rnd_sinsqth13_wrc;
    double rnd_sinsqth23_wrc;
    double rnd_sinsq2th23_wrc;
    double rnd_sinsq2th23_mod_wrc;
    
    double sinsqth23predval;
    double sinsq2th23predval;
    
    Int_t NBINTH23 = hhist1d[8]->GetNbinsX();
    Double_t* xbinsTH23 = new Double_t[NBINTH23+1];
    TAxis *axisTH23 = hhist1d[8]->GetXaxis();
    for (Int_t ibin=0;ibin<NBINTH23;ibin++)  {
        xbinsTH23[ibin] = axisTH23->GetBinLowEdge(1+ibin);
    }
    xbinsTH23[NBINTH23]=axisTH23->GetXmax();
    
    Int_t NBINTH13 = hhist1d[11]->GetNbinsX();
    Double_t* xbinsTH13 = new Double_t[NBINTH13+1];
    TAxis *axisTH13 = hhist1d[11]->GetXaxis();
    for (Int_t ibin=0;ibin<NBINTH13;ibin++)  {
        xbinsTH13[ibin] = axisTH13->GetBinLowEdge(1+ibin);
    }
    xbinsTH13[NBINTH13]=axisTH13->GetXmax();
    
//  
    TH1D *hsinsqth13_pred = new TH1D("hsinsqth13_pred","",NBINTH13,xbinsTH13); //Change from 23 --> 13
    TH1D *hsinsqth13_worc = new TH1D("hsinsqth13_worc","",NBINTH13,xbinsTH13);
    TH1D *hsinsqth13_wrc = new TH1D("hsinsqth13_wrc","",NBINTH13,xbinsTH13);
    
    TH1D *hsinsqth23_worc = new TH1D("hsinsqth23_worc","",NBINTH23,xbinsTH23);
    TH1D *hsinsqth23_pred = new TH1D("hsinsqth23_pred","",NBINTH23,xbinsTH23);
    TH1D *hsinsqth23_wrc = new TH1D("hsinsqth23_wrc","",NBINTH23,xbinsTH23);
    
    TH1D *hsinsq2th23_worc = new TH1D("hsinsq2th23_worc","",NBINTH23,0.95,1.);
    TH1D *hsinsq2th23_wrc = new TH1D("hsinsq2th23_wrc","",NBINTH23,0.95, 1.);
    
    TH1D *hsinsq2th23_mod_worc = new TH1D("hsinsq2th23_mod_worc","",NBINTH23,0.9,1.);
    TH1D *hsinsq2th23_mod_wrc = new TH1D("hsinsq2th23_mod_wrc","",NBINTH23,0.9, 1.);
//  Monte Carlo Experiment 
    
    for (Int_t ithrow=0; ithrow<NTHROW; ++ithrow) {
        rnd_sinsqth13_worc = hhist1d[11]->GetRandom();
        hsinsqth13_worc->Fill(rnd_sinsqth13_worc);
        
        rnd_sinsqth23_worc = hhist1d[8]->GetRandom();
        hsinsqth23_worc->Fill(rnd_sinsqth23_worc);
        hsinsq2th23_worc->Fill(4*rnd_sinsqth23_worc*(1-rnd_sinsqth23_worc));
        
        rnd_sinsq2th23_worc = 4*rnd_sinsqth23_worc*(1-rnd_sinsqth23_worc);
        hsinsq2th23_worc->Fill(rnd_sinsq2th23_worc);
        
       //%% 
        rnd_sinsq2th23_mod_worc = rnd_sinsq2th23_worc-rnd_sinsqth23_worc*4*rnd_sinsqth13_worc*(1-rnd_sinsqth13_worc)*(1-2*rnd_sinsqth23_worc);
        hsinsq2th23_mod_worc->Fill(rnd_sinsq2th23_mod_worc);
    //%%%
        rnd_sinsqth13_wrc = hhist1d[23]->GetRandom();
        hsinsqth13_wrc->Fill(rnd_sinsqth13_wrc);
        
        rnd_sinsqth23_wrc = hhist1d[20]->GetRandom();
        hsinsqth23_wrc->Fill(rnd_sinsqth23_wrc);
        rnd_sinsq2th23_wrc = 4*rnd_sinsqth23_wrc*(1-rnd_sinsqth23_wrc);
        hsinsq2th23_wrc->Fill(rnd_sinsq2th23_wrc);

        rnd_sinsq2th23_mod_wrc = rnd_sinsq2th23_wrc-rnd_sinsqth23_wrc*4*rnd_sinsqth13_wrc*(1-rnd_sinsqth13_wrc)*(1-2*rnd_sinsqth23_wrc);
        hsinsq2th23_mod_wrc->Fill(rnd_sinsq2th23_mod_wrc);
        
        
        sinsqth23predval = rnd_sinsqth23_worc*rnd_sinsqth13_worc*(1-rnd_sinsqth13_worc)/(rnd_sinsqth13_wrc*(1-rnd_sinsqth13_wrc));
        
        sinsq2th23predval = 4*sinsqth23predval*(1-sinsqth23predval);
        
        if(sinsq2th23predval>0.98)hsinsqth23_pred->Fill(sinsqth23predval);
        if (ithrow<10) {
            cout<<"th13_worc "<<rnd_sinsqth13_worc<<" val "<<sinsqth23predval<<endl;
        }
    }
    hsinsqth13_wrc->GetXaxis()->SetTitle("sin^{2}#theta_{13}");
    hsinsqth13_wrc->GetYaxis()->SetTitle("No. of Throws");
    plot2hist(hsinsqth13_wrc,"w/ reactor ", hsinsqth13_worc,"wo/ reactor",true,Form("%s_th13_rand",output.Data()),0.6,0.6);
    
    hsinsqth23_wrc->GetXaxis()->SetTitle("sin^{2}#theta_{23}");
    hsinsqth23_wrc->GetYaxis()->SetTitle("No. of Throws");
    plot2hist(hsinsqth23_wrc,"w/ reactor ", hsinsqth23_worc,"wo/ reactor",true,Form("%s_th23_rand",output.Data()),0.6,0.6);
    
    plot3hist(hsinsqth23_wrc,"w/ reactor ", hsinsqth23_worc,"wo/ reactor",hsinsqth23_pred,"prediction",false,Form("%s_th23_pred",output.Data()),0.6,0.6);
    
    hsinsq2th23_wrc->GetXaxis()->SetTitle("sin^{2}2#theta_{23}");
    hsinsq2th23_wrc->GetYaxis()->SetTitle("No. of Throws");
    plot2hist(hsinsq2th23_wrc,"w/ reactor ", hsinsq2th23_worc,"wo/ reactor",true, Form("%s_2th23_rand",output.Data()),0.6,0.6);
    
    hsinsq2th23_mod_wrc->GetXaxis()->SetTitle("sin^{2}2#theta_{23}^{eff.}");
    hsinsq2th23_mod_wrc->GetYaxis()->SetTitle("No. of Throws");
    plot2hist(hsinsq2th23_mod_wrc,"w/ reactor ", hsinsq2th23_mod_worc,"wo/ reactor",true, Form("%s_2th23_mod_rand",output.Data()),0.6,0.6);
    
}

void plot3hist(TH1* h1, TString leg1, TH1* h2, TString leg2, TH1* h3, TString leg3, Bool_t isTitleStyle, TString savename="", Float_t xlegmin=0.4, Float_t ylegmin=0.4 ){
    new TCanvas;
    gStyle->SetOptStat(0);
    gPad->SetBottomMargin(gPad->GetBottomMargin()*1.2);
    Double_t ymax = h1->GetMaximum();
    if (ymax<h2->GetMaximum()) {
        ymax = h2->GetMaximum();
    }
    if (ymax<h3->GetMaximum()) {
        ymax = h3->GetMaximum();
    }
    Double_t ymin = 0;
    
    if(isTitleStyle) titleStyle(h1);
    h1->SetLineColor(2);
    h1->SetLineWidth(2);
    h1->GetYaxis()->SetRangeUser(ymin,ymax*1.2);
    h1->DrawCopy("");
    
    h2->SetLineColor(4);
    h2->SetLineWidth(2);
    h2->DrawCopy("same");
    
    h3->SetLineColor(6);
    h3->SetLineWidth(2);
    h3->DrawCopy("same");
    
    TLegend* leg = new TLegend(xlegmin,ylegmin,xlegmin+0.3,ylegmin+0.25);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.06);
    leg->AddEntry(h1, leg1.Data(),"l");
    leg->AddEntry(h2, leg2.Data(),"l");
    leg->AddEntry(h3, leg3.Data(),"l");
    leg->Draw();
    
    gPad->Print(savename+".eps");
    //gPad->Print(savename+".png");
    
}


void plot2hist(TH1* h1, TString leg1, TH1* h2, TString leg2, Bool_t isTitleStyle, TString savename="", Float_t xlegmin=0.4, Float_t ylegmin=0.4 ){
    new TCanvas;
    gStyle->SetOptStat(0);
    gPad->SetBottomMargin(gPad->GetBottomMargin()*1.2);
    Double_t ymax = h1->GetMaximum();
    if (ymax<h2->GetMaximum()) {
        ymax = h2->GetMaximum();
    }
    Double_t ymin = 0;
    
    if(isTitleStyle) titleStyle(h1);
    h1->SetLineColor(2);
    h1->SetLineWidth(2);
    h1->GetYaxis()->SetRangeUser(ymin,ymax*1.2);
    h1->DrawCopy("");
    
    h2->SetLineColor(4);
    h2->SetLineWidth(2);
    h2->DrawCopy("same");
    
    TLegend* leg = new TLegend(xlegmin,ylegmin,xlegmin+0.3,ylegmin+0.25);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.06);
    leg->AddEntry(h1, leg1.Data(),"l");
    leg->AddEntry(h2, leg2.Data(),"l");
    leg->Draw();
    
    gPad->Print(savename+".eps");
    //gPad->Print(savename+".png");
    
}

void titleStyle(TH1* h1){
    h1->SetTitle("");
    h1->GetYaxis()->CenterTitle();
    h1->GetXaxis()->CenterTitle();
    h1->GetXaxis()->SetLabelSize(h1->GetXaxis()->GetTitleSize()*1.2);
    h1->GetYaxis()->SetLabelSize(h1->GetYaxis()->GetTitleSize()*1.2);
    h1->GetXaxis()->SetTitleSize(h1->GetXaxis()->GetLabelSize()*1.2);
    h1->GetYaxis()->SetTitleSize(h1->GetYaxis()->GetLabelSize()*1.2);
    h1->GetYaxis()->SetTitleOffset(1.);
}
