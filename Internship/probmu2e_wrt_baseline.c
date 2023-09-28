#include "TH1.h"
#include "TMath.h"
#include "TF1.h"
#include "TLegend.h"
#include "TCanvas.h"


double total(double Ene, double dmsq_23, double dmsq_21, double baseL, double theta_12, double theta_23, double theta_13, double delta, double density) {
        double probmu2e = 4*pow(sin(theta_13),2)*pow(sin(theta_23),2)*pow(cos(theta_13),2)*pow(sin(1.27*dmsq_23*baseL/Ene),2) - 8*pow(sin(theta_13),2)*pow(sin(theta_23),2)*pow(cos(theta_13),2)*(density*0.6/dmsq_23)*(2*pow(sin(theta_13), 2) - 1)*pow(sin(1.27*dmsq_23*baseL/Ene),2)

+ 8*pow(sin(theta_13),2)*pow(sin(theta_23),2)*pow(cos(theta_13),2)*(density*0.6*baseL/(4*Ene))*(2*pow(sin(theta_13), 2) - 1)*sin(1.27*dmsq_23*baseL/Ene)*cos(1.27*dmsq_23*baseL/Ene)

- 8*sin(theta_12)*sin(theta_13)*sin(theta_23)*cos(theta_12)*pow(cos(theta_13),2)*cos(theta_23)*sin(delta)*pow(sin(1.27*dmsq_23*baseL/Ene),2)*sin(1.27*dmsq_21*baseL/Ene)

+ 8*sin(theta_12)*sin(theta_13)*sin(theta_23)*pow(cos(theta_13),2)*(cos(theta_12)*cos(theta_23)*cos(delta) - sin(theta_12)*sin(theta_13)*sin(theta_23))*sin(1.27*dmsq_21*baseL/Ene)*sin(1.27*dmsq_23*baseL/Ene)*cos(1.27*dmsq_23*baseL/Ene)

+ 4*pow(sin(theta_12),2) *pow(cos(theta_13),2)*(pow(cos(theta_12),2)*pow(cos(theta_23),2) + pow(sin(theta_12),2)*pow(sin(theta_13),2)*pow(sin(theta_23),2) - 2*sin(theta_12)*sin(theta_13)*sin(theta_23)*cos(theta_12)*cos(theta_23)*cos(delta))*pow(sin(1.27*dmsq_21*baseL/Ene),2);

return probmu2e;
    }  
void probmu2e_wrt_baseline(){
    TCanvas *c1 = new TCanvas("c1","",200,10,700,500);
    double myEne;// = 0.6;//T2K peak energy in GeV
    double myDmsq23 = 2.5e-3;//in eV^2
    double myDmsq21 = 7.6e-5;//in eV^2
    double myBaseL = 295;//T2K km
    double myTheta12 = asin(sqrt(0.8704))/2; //in radian 45*Pi/180
    double myTheta23 = asin(sqrt(0.5)); //in radian 45*Pi/180
    double myTheta13 = asin(sqrt(0.085))/2; //in radian 45*Pi/180
    double myDelta = -TMath::Pi()/2; //in radian
    double myDensity = 2.05e-4; //t2k matter effect
    	    

    double myProbmu2ef = total(myEne, myDmsq23,myDmsq21, myBaseL,myTheta12,myTheta23,myTheta13,myDelta, myDensity);
    //cout<<"probably at the peak"<<myProbmu2ef<<endl;
    //set number of steps you want to calculate
    const Int_t nEnStep = 1000;
    // set maximal energy to plot
    double EnergyMax = 2.0;
    //double BaseLMax = 1500;
    //Graph of oscillation prob. as function of energy
    //Note:formula is not work at zero energy
    TGraph *pGraphMu2ef;
    Double_t pBaseL[nEnStep], pProb[nEnStep];

    for (Int_t istep=1; istep<=nEnStep; ++istep){
        //energy value at each step
        //pEne[istep-1] = BaseLMax*istep/nEnStep;
        pBaseL[istep-1] = EnergyMax*istep/nEnStep;
        //calculated osc. prob. at each step
        pProb[istep-1] = total(pBaseL[istep-1],myDmsq23,myDmsq21, myBaseL,myTheta12,myTheta23,myTheta13,myDelta, myDensity);
       // cout<<"step "<<istep<<" Energy "<<pEne[istep-1]<<" Prob. "<<pProb[istep-1]<<endl;
    }
    //create new Graph based on the array
    pGraphMu2ef = new TGraph(nEnStep,pBaseL,pProb);
    //set the width for graph line
    pGraphMu2ef->SetLineWidth(2);  
    //create new Canvas
   // new TCanvas;
    
    pGraphMu2ef->SetTitle("Oscillation probability w.r.t baseline");    //set title for the graph
    pGraphMu2ef->GetXaxis()->SetTitle("Energy[eV]");
    pGraphMu2ef->GetYaxis()->SetTitle("Osillation Probability #nu_{#mu}#rightarrow #nu_e");
    pGraphMu2ef->GetXaxis()->SetRangeUser(0.0, 1500);
    pGraphMu2ef->GetYaxis()->SetRangeUser(0,0.8);
    //pGraphMu2ef->GetXaxis()->CenterTitle();
   // pGraphMu2ef->GetYaxis()->CenterTitle();
    pGraphMu2ef->GetYaxis()->SetTitleOffset(1.2);
    pGraphMu2ef->SetLineColor(2);
  //  pGraphMu2e->SetLineStyle(7);
    pGraphMu2ef->Draw("AL");//draw graph
 /* Draw horizontal lines   */
    c1->Update();
    TLine *l1 = new TLine(c1->GetUxmin(), 0.0, c1->GetUxmax(), 0.0);
    l1->SetLineColor(1);
    l1->Draw("L SAME");

    TLegend * pleg = new TLegend(0.7,0.65,0.85,0.87);
    //pleg -> AddEntry( pGraphMu2ef , "Total" ,"l");
    //pleg -> AddEntry( pGraphMu2e , "Leading" ,"l");
    //pleg -> AddEntry( pGraphMu2eb, "Matter effect" ,"l");
    //pleg -> AddEntry( pGraphMu2ec , "CP-violating" ,"l");
    //pleg -> AddEntry( pGraphMu2ed , "CP-conserving" ,"l");
    //pleg -> AddEntry( pGraphMu2ee , "Solar" ,"l");

    pleg -> SetFillColor(0);
    pleg -> SetBorderSize(0);
    pleg ->SetTextSize(14);
    pleg ->SetTextFont(43);
    pleg -> Draw();


 gStyle->SetOptStat(0);
 pleg = new TLegend(0.6,0.5,0.8,0.6);
 pleg->SetTextSize(0.03);
 pleg->SetTextColor(1);
 pleg->SetFillColor(0);
 pleg->SetBorderSize(0);
 pleg ->SetTextSize(14);
 pleg ->SetTextFont(43);
 pleg->AddEntry(pGraphMu2ef,"L = 295km","");
 pleg->AddEntry(pGraphMu2ef,"#delta_{CP} = -#pi/2","");
 pleg->Draw();
   //save the graph
    gPad->Print("prob_mu2e_wrt_baseline_t2k.pdf");//pdf format
 
}



