#include "RiceStyle.h"
using namespace std;
void getImage()
{
    //input file from analyzing the output, which has t distribution from MC and REC.
    const char* file = "results.root";
    TFile* input = new TFile(file);

    TH1D* hdsigmadt_MC = (TH1D*)input->Get("h_tMC");
    TH1D* hdsigmadt_REC = (TH1D*)input->Get("h_tREC");

    int nbins = hdsigmadt_MC->GetNbinsX();
    double dsigmadt_MC,dsigmadt_REC, tBinWidth,t,b,delta, F_b_MC,F_b_REC, result1=0,result2=0;

    // ************************************
    double t_cut = 0.25; //upper integration limit
    // ************************************

    double bmin= -12;
    double bmax= 12;
    double noOfBins = 300;
    double hbarc = 0.197;

    //define the b histos for MC and REC
    TH1D* hF_b_MC = new TH1D("hF_b_MC", "F_b_MC", noOfBins, bmin, bmax);
    TH1D* hF_b_REC = new TH1D("hF_b_REC", "F_b_MC", noOfBins, bmin, bmax);

    //numerical integration.
    /* reference, https://arxiv.org/pdf/1211.3048, page 8, Eq(22).  
    
    Hints:
      1.  \delta^2 = -t, the measured cross section is in t, but integration is in \delta,
        we need to do a change of variable.
      2.  hbarc is your good friend for converting to dimensionless quantity. 
          what is the unit of hbarc?  
      3. Finding the charge radius of the ion can be approached differently. 
          For this exercise, let's use the `Width at Half Maximum` as a way to 
          obtain the R value of the ion. 

    */
    for (int j=1; j<=noOfBins; j++)
    {
        F_b_MC= 0;F_b_REC=0;
        b = hF_b_MC->GetBinCenter(j);
        double prefactor = 1/6.28;

        for (int i=1; i<=nbins; i++)
        {
            tBinWidth = hdsigmadt_MC->GetBinWidth(i);
            t = hdsigmadt_MC->GetBinCenter(i); // GeV2
            delta =  sqrt(fabs(t)); // GeV
                //hint: 2\delta d\delta = dt
            
            dsigmadt_MC = hdsigmadt_MC->GetBinContent(i);  //nb/GeV^2
            dsigmadt_REC = hdsigmadt_REC->GetBinContent(i);  //nb/GeV^2
            
            dsigmadt_MC/=1e7; //convert to fm^2/GeV2
            dsigmadt_REC/=1e7; //convert to fm^2/GeV2
           
            //hint: your task is to figure out what the below 2 (3) quantities are. 
                // see later result1 and result2 what you need.
            double bessel=999.;//TMath::BesselJ0(); remember the quantity is dimensionless
            double amp_MC = 999.;//
            double amp_REC = 999.; //
            
            //this is the part 
            //we make modified amplitude in Eq.(22), 
            //by flipping the sign at each minimum point.
            if(t>t_cut)
                continue;
            if(t>0.0414)  { //1st minima
                amp_MC*=-1;
                amp_REC*=-1;
            }
            if(t>0.135)  {   //2nd minima
                amp_MC*=-1;
                amp_REC*=-1;
            }

            result1 =  amp_MC * bessel * tBinWidth/2 ;
            result2 =  amp_REC * bessel * tBinWidth/2 ;

            //perform the sum         
            F_b_MC += result1;
            F_b_REC += result2;
            
        }
        
        F_b_MC*=prefactor;F_b_MC/=hbarc;
        F_b_REC*=prefactor;F_b_REC/=hbarc;

        //for this exercise, we ignore the errors and just set to 1%
        hF_b_MC->SetBinContent(j, F_b_MC);hF_b_MC->SetBinError(j, F_b_MC*0.01);
        hF_b_REC->SetBinContent(j, F_b_REC);hF_b_REC->SetBinError(j, F_b_MC*0.01);
    }

    //plotting. you need to plot your results.
   
    //normalizing the F(b) to become profile.
    hF_b_MC->Scale(1.0 / hF_b_MC->Integral("width"));
    hF_b_REC->Scale(1.0 / hF_b_REC->Integral("width"));


    // define a canvas and then draw..(just in case you don't know ROOT;-))
    
    // hF_b_MC->Draw("same");
    // hF_b_REC->Draw("Psame");

    TLegend *leg = new TLegend(0.7, 0.75, 0.9, 0.85); // (x1, y1, x2, y2)
    leg->AddEntry(hF_b_MC, "MC t_{max} = 0.25", "P"); // change to your t_max for integration.
    leg->AddEntry(hF_b_REC, "REC t_{max} = 0.25", "P"); // change to your t_max for integration.
    leg->SetBorderSize(0);      // No border box
    leg->SetFillStyle(0);       // Transparent background
    leg->SetTextFont(42);       // Nice readable font
    leg->SetTextSize(0.03);     // Optional: adjust size
    // leg->Draw("same");

    cout << "True ion size = " <<  " [your answer] fm " << endl;
    cout << "Reco ion size = " <<  " [your answer] fm " << endl;
    cout << "All done. Bye." << endl;

}



