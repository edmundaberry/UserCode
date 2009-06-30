#include <iostream>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TFile.h>
#include <TLine.h>

//-------------------------------------------------
// User set information 
//-------------------------------------------------

const int  MAX_N_EVENTS = -1;
const int  GUN_ENERGY   = 100;
const bool SAVE_PLOTS   = true;  

//-------------------------------------------------
// Plot names
//-------------------------------------------------

const char* CCT_TOTAL_EHAD   = "cct_totalEHad";
const char* CCT_TOTAL_EEM    = "cct_totalEEm";
const char* CCT_TOTAL_E      = "cct_totalE";
			     
const char* DCT_TOTAL_EHAD   = "dct_totalEHad";
const char* DCT_TOTAL_EEM    = "dct_totalEEm";
const char* DCT_TOTAL_E      = "dct_totalE";
			     
const char* TPG_TOTAL_EHAD   = "tpg_totalEHad";
const char* TPG_TOTAL_EEM    = "tpg_totalEEm";
const char* TPG_TOTAL_E      = "tpg_totalE";

const char* RATIO_TOTAL_EHAD = "ratio_totalEHad";
const char* RATIO_TOTAL_EEM  = "ratio_totalEEm";
const char* RATIO_TOTAL_E    = "ratio_totalE";

//-------------------------------------------------
// All other Global Variables
//-------------------------------------------------

char PLOT_FILE_NAME[500];
char PARTICLE_NAME [100];
char POSITION_NAME [100];

TFile *PLOT_FILE;

//-------------------------------------------------
// Set the file names
//-------------------------------------------------

void setFileNames(){

  sprintf( PLOT_FILE_NAME, "Plots.root");
  
}


bool style_energy_comp( const char* plotName_1, const char* legend_entry_1,
			const char* plotName_2, const char* legend_entry_2,
			const char* plotName_3, const char* legend_entry_3,
			const char *save_name,
			const char* title, const bool draw_gun_line  ){
  
  TH1F *plot1 = (TH1F*) PLOT_FILE -> Get(plotName_1);
  TH1F *plot2 = (TH1F*) PLOT_FILE -> Get(plotName_2);
  TH1F *plot3 = (TH1F*) PLOT_FILE -> Get(plotName_3);

  if (!plot1) {
    std::cout << "Could not get a TH1F of the name " << plotName_1 << std::endl;
    std::cout << "from the file at " << PLOT_FILE -> GetName() << std::endl;
    return false;
  }

  if (!plot2) {
    std::cout << "Could not get a TH1F of the name " << plotName_2 << std::endl;
    std::cout << "from the file at " << PLOT_FILE -> GetName() << std::endl;
    return false;
  }
  
  if (!plot3) {
    std::cout << "Could not get a TH1F of the name " << plotName_3 << std::endl;
    std::cout << "from the file at " << PLOT_FILE -> GetName() << std::endl;
    return false;
  }

  TCanvas *canvas = new TCanvas();
  canvas -> cd();

  char yTitle[100];
  sprintf(yTitle, "Counts per %1.2f GeV", plot1 -> GetBinWidth(1));

  plot3 -> Draw();
  plot2 -> Draw("Same");
  plot1 -> Draw("Same");

  plot1 -> GetXaxis() -> SetTitle("Energy [GeV]");
  plot2 -> GetXaxis() -> SetTitle("Energy [GeV]");
  plot3 -> GetXaxis() -> SetTitle("Energy [GeV]");

  plot1 -> GetYaxis() -> SetTitle(yTitle);
  plot2 -> GetYaxis() -> SetTitle(yTitle);
  plot3 -> GetYaxis() -> SetTitle(yTitle);

  plot1 -> SetLineColor(kGreen);
  plot1 -> SetLineWidth(3.0);

  plot2 -> SetLineColor(kViolet);
  plot2 -> SetLineWidth(3.0);

  plot3 -> SetLineColor(kBlue);
  plot3 -> SetLineWidth(6.0);

  float plot1_max = plot1 -> GetMaximum();
  float plot2_max = plot2 -> GetMaximum();
  float plot3_max = plot3 -> GetMaximum();
  
  float temp = TMath::Max(plot1_max, plot2_max);
  float max  = TMath::Max(temp, plot3_max);

  plot1 -> SetMaximum(1.04 * max);
  plot2 -> SetMaximum(1.04 * max);
  plot3 -> SetMaximum(1.04 * max);
  
  plot3 -> Draw();
  plot2 -> Draw("Same");
  plot1 -> Draw("Same");

  canvas -> Update();
  double xmin = plot1 -> GetXaxis() -> GetXmin();
  double ymin = gPad  -> GetUymin();
  double ymax = gPad  -> GetUymax();

  TLine *gun_line = new TLine(2* GUN_ENERGY, ymin, 2* GUN_ENERGY, ymax);
  gun_line -> SetLineStyle(kDashed);
  gun_line -> SetLineWidth(3.0);
  gun_line -> SetLineColor(kRed);
  
  if (draw_gun_line) gun_line -> Draw("same");

  TLegend *leg = new TLegend(0.6,0.75,0.875,0.87);
  leg -> SetFillColor(0);
  leg -> AddEntry(plot1, legend_entry_1, "l");
  leg -> AddEntry(plot2, legend_entry_2, "l");
  leg -> AddEntry(plot3, legend_entry_3, "l");
  if (draw_gun_line) leg -> AddEntry(gun_line ,"Actual energy from particle gun", "l");

  TLatex *titleLatex = new TLatex();
  titleLatex -> SetTextFont (42);
  titleLatex -> SetTextAlign (12);
  titleLatex -> SetTextSize(0.055);

  titleLatex -> DrawLatex(xmin,ymin + (ymax - ymin)*1.08,title);
  
  leg -> Draw();

  if (SAVE_PLOTS) canvas -> SaveAs(save_name);

  return true;
}

bool setFiles(){

  PLOT_FILE = new TFile(PLOT_FILE_NAME);
  
  if (!PLOT_FILE){
    std::cout << "Could not get plot file." << std::endl;
    std::cout << "I looked here:" << std::endl;
    std::cout << PLOT_FILE_NAME << std::endl;
    return false;
  }

  return true;

}


void make_plots(){

  char title[500];
  char titleDetail[500]; 

  setFileNames();

  bool filesSet = setFiles(); if (!filesSet) return;
  
  sprintf(titleDetail,"- 2x%d GeV electron gun", GUN_ENERGY);


  sprintf(title,"Total event hadronic energy comparison %s", titleDetail);
  style_energy_comp ( CCT_TOTAL_EHAD, "Energy from created calotowers",
		      DCT_TOTAL_EHAD, "Energy from default calotowers",
		      TPG_TOTAL_EHAD, "Energy from trigger primitives", 
		      "hadEnergyComparison.gif",
		      title, false);

  sprintf(title,"Total event E&M energy comparison %s", titleDetail);
  style_energy_comp ( CCT_TOTAL_EEM, "Energy from created calotowers",
		      DCT_TOTAL_EEM, "Energy from default calotowers",
		      TPG_TOTAL_EEM, "Energy from trigger primitives", 
		      "emEnergyComparison.gif",
		      title, false);
  
  sprintf(title,"Total event energy comparison %s", titleDetail);
  style_energy_comp ( CCT_TOTAL_E, "Energy from created calotowers",
		      DCT_TOTAL_E, "Energy from default calotowers",
		      TPG_TOTAL_E, "Energy from trigger primitives", 
		      "totalEnergyComparison.gif",
		      title, true);

}
