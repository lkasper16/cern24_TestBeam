#define trdclass_cern24_cxx
#include "trdclass_cern24.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TSocket.h"
#include "TMarker.h"
#include "TMultiGraph.h"
#include "PlotLib.C"
#include <TError.h>
#include <TStopwatch.h>
#include "GNN/gnn_model.h"
#include "GNN/gnn_model.cpp"
#include "GNN/toGraph.cpp"
#include <iostream>
#include <fstream>

#define NPRT 1000
#define USE_TRK
#define  MAX_PRINT 10
//
#define USE_GNN  1
#define USE_FIT  1
#define USE_CLUST 1
#define USE_PULSE 0
//
#define USE_125_RAW
#define USE_250_PULSE
#define MAX_CLUST 500
//
#define SAVE_TRACK_HITS
#define SAVE_PDF
#define WRITE_CSV
#define DEBUG 0
//
//-- For single evt clustering display, uncomment BOTH:
//#define SHOW_EVTbyEVT
//#define SHOW_EVT_DISPLAY

void WriteToCSV(std::ofstream &csvFile, float v1, float v2, float v3, float v4, float v5, float v6, float v7, float v8, float v9, float v10) {
  csvFile<<v1<<","<<v2<<","<<v3<<","<<v4<<","<<v5<<","<<v6<<","<<v7<<","<<v8<<","<<v9<<","<<v10<<std::endl;
}

//===================================
//      TRD DAQ Channel Mapping
//===================================

// -- GEMTRD mapping --
int GetGEMChan(int ch, int slot) {
  int cardNumber = ch/24;
  int cardChannel = ch-cardNumber*24;
  float dchan = cardChannel+cardNumber*24+(slot-3)*72.;
  if (slot==5 || slot==6 || slot==7 || (slot==8 && ch<24)) {
    if ((dchan-144.)==224. || (dchan-144.)==207. || (dchan-144.)==208.) { return -1; } else {return dchan-144.;}
  }
  return -1;
}

// -- MMG-1 mapping --
int GetMMG1Chan(int ch, int slot, int runNum) {
  int cardNumber = ch/24;
  int cardChannel = ch-cardNumber*24;
    float dchan = cardChannel+cardNumber*24+(slot-3)*72.;
  if (runNum>4450) {
    if (slot==9 || slot==10 || (slot==8&&ch>23)) {
      if ((dchan-384.)==16. || (dchan-384.)==45.) { return -1; } else {return dchan - 384.;}
    }
    if (slot==13&&ch<48) {
      if ((dchan-528.)==224.) { return -1; } else {return dchan - 528.;}
    }
  }
  return -1;
}
//============ END DAQ Channel Mapping ============

void trdclass_cern24::Loop() {
  
  if (fChain == 0) return;

  //==================================================================================================
  //            B o o k    H i s t o g r a m s
  //==================================================================================================
  
  gErrorIgnoreLevel = kBreak; // Suppress warning messages from empty chi^2 fit data
  TList *HistList = new TList();
  #ifdef WRITE_CSV
    std::ofstream csvFile("EventByEvent.csv");
  #endif
  //-----------------  (canvas 0) Event Display ----------
  #ifdef SHOW_EVT_DISPLAY
    char c0Title[256]; sprintf(c0Title,"Event_Display_Run=%d",RunNum);
    TCanvas *c0 = new TCanvas("DISP",c0Title,1100,200,1500,1300);
    c0->Divide(4,3); c0->cd(1);
    //-----------------  canvas 2 FPGA Display ----------
    char c2Title[256]; sprintf(c2Title,"FPGA_Event_Display_Run=%d",RunNum);
    TCanvas *c2 = new TCanvas("FPGA",c2Title,1000,100,1500,1300);
    c2->Divide(5,2); c2->cd(1);
    char c3Title[256]; sprintf(c3Title,"DQM_Event_Display_Run=%d",RunNum);
    TCanvas *c3 = new TCanvas("PID",c3Title,1000,100,1500,1300);
    c3->Divide(5,1); //c3->cd(1);
  #endif
  
  // Track fit in time
  TF1 fx1("fx1","pol1",40,150);
  TF1 fx2("fx2","pol1",40,150);
  f125_fit = new TH2F("f125_fit","GEM-TRD Track Fit; Time Response (8ns) ; Channel ",250,0.5,250.5,256,-0.5,255.5);
  mmg1_f125_fit = new TH2F("mmg1_f125_fit","MMG1-TRD Track Fit; Time Response (8ns) ; Channel ",250,0.5,250.5,256,-0.5,255.5);
  //-- TRD - GEMTRKR alignment --------
  double xx1=-37., yy1=-55.,  xx2=53., yy2=44.;
  double aa=(yy2-yy1)/(xx2-xx1);
  double bb=yy1-aa*xx1;
  TF1 ftrk("ftrk","[0]*x+[1]",-55.,55.);
  ftrk.SetParameter(0,aa);
  ftrk.SetParameter(1,bb);
  TF1 ftrkr("ftrk","(x-[1])/[0]",0.,255.);
  ftrkr.SetParameter(0,aa);
  ftrkr.SetParameter(1,bb);
  //-----------------------
  //float x1 = 0;
  float z1 = 0;
  //float x2 = 0;
  float z2 = 1071;
  //float x3 = 0;
  float z3 = 1153.4;
  //float xgem = 0;
  float zgem = 571;
  //float xmmg1 = 0;
  float zmmg1 = 308;
  
  //----------------  GEM TRK fuducial area selection (box cut) ----------------------
  double xbc1=-50., xbc2=+50., ybc1=-50., ybc2=+50.;
  //double gemtrk_x2ch=-999.;
  //TLine peak_line[100];
  //----------------------------------------------------------------------------------
  
  hcount= new TH1D("hcount","Count",3,0,3);                                     HistList->Add(hcount);
  hcount->SetStats(0);   hcount->SetFillColor(38);   hcount->SetMinimum(1.);
  #if ROOT_VERSION_CODE > ROOT_VERSION(6,0,0)
    hcount->SetCanExtend(TH1::kXaxis);
  #else
    hcount->SetBit(TH1::kCanRebin);
  #endif
  
  h250_size = new TH1F("h250_size"," fa250 Raw data size",4096,0.5,4095.5);                                        HistList->Add(h250_size);
  int nx0=120;      int ny0=256; //528;
  double Ymin=0.;    double Ymax=ny0*0.4; // +528.;
  double Xmin=0.;    double Xmax=30.;
  mhevt  = new TH2F("mhevt","MMG1-TRD Event display; z pos,mm; y pos,mm ",nx0+100.,Xmin,Xmax,ny0,Ymin,Ymax); mhevt->SetStats(0); mhevt->SetMaximum(10.); mhevt->SetMinimum(-5.);
  mhevti = new TH2F("mhevti","MMG1: ML-FPGA response; z pos,mm; y pos,mm ",nx0+100.,Xmin,Xmax,ny0,Ymin,Ymax);  mhevti->SetStats(0); mhevti->SetMaximum(10.);
  mhevtf = new TH2F("mhevtf","MMG1: Clusters for FPGA ; z pos,mm; y pos,mm ",nx0+100.,Xmin,Xmax,ny0,Ymin,Ymax);  mhevtf->SetStats(0); mhevtf->SetMaximum(10.);

  hevt  = new TH2F("hevt","GEM-TRD Event display; z pos,mm; y pos,mm ",nx0,Xmin,Xmax,ny0,Ymin,Ymax); hevt->SetStats(0); hevt->SetMaximum(10.); hevt->SetMinimum(-5.);
  hevtc = new TH2F("hevtc"," Clustering ; FADC bins; GEM strips",nx0,-0.5,nx0-0.5,ny0,-0.5,ny0-0.5);
  hevtc->SetStats(0);   hevtc->SetMinimum(0.07); hevtc->SetMaximum(40.);
  hevti = new TH2F("hevti","GEM: ML-FPGA response; z pos,mm; y pos,mm ",nx0,Xmin,Xmax,ny0,Ymin,Ymax);  hevti->SetStats(0); hevti->SetMaximum(10.);
  hevtf = new TH2F("hevtf","GEM: Clusters for FPGA ; z pos,mm; y pos,mm ",nx0,Xmin,Xmax,ny0,Ymin,Ymax);  hevtf->SetStats(0); hevtf->SetMaximum(10.);
  //hevtL = new TH2F("hevtL"," test ; z pos,mm; y pos,mm ",nx0,0.,200,ny0,-0.,+530.);  hevtL->SetStats(0); hevtL->SetMaximum(10.);
  hevtk  = new TH2F("hevtk"," Event display; z pos,mm; y pos,mm ",nx0,Xmin,Xmax,ny0,Ymin,Ymax); /*hevtk->SetStats(0);*/ hevtk->SetMaximum(10.);
  hevtck = new TH2F("hevtck"," Clustering ; FADC bins; GEM strips",nx0,-0.5,nx0-0.5,ny0,-0.5,ny0-0.5);

  aver2d_e = new TH2F("aver2d_e","aver-rms ",100,0.,240., 100,0.,10.);   HistList->Add(aver2d_e);
  aver2d_p = new TH2F("aver2d_p","aver-rms",100,0.,240., 100,0.,10.);    HistList->Add(aver2d_p);

  //-- Calorimeter
  hCal_pulse  = new TH1F("hCal_pulse"," Calorimeter pulse ",300,0.,300.);                                         HistList->Add(hCal_pulse);
  hCal_sum    = new TH1F("hCal_sum"," Calorimeter Sum (GeV)",200,0.,6000.);                                       HistList->Add(hCal_sum);
  hCal_sum_el = new TH1F("hCal_sum_el"," Calorimeter Sum for electrons",200,0.,6000.);                            HistList->Add(hCal_sum_el);
  hCal_sum_pi = new TH1F("hCal_sum_pi"," Calorimeter Sum for pions",200,0.,6000.);                                HistList->Add(hCal_sum_pi);
  hCal_occ = new TH1F("hCal_occ"," Calorimeter Occupancy",9,-0.5,8.5);                                            HistList->Add(hCal_occ);
  hCal_Cher  = new TH2F("hCal_Cher"," Calorimeter vs Cherenkov ; Cherenkov ; Calorimeter ",100,0.,2000.,100,0.,6000.);                      HistList->Add(hCal_Cher);
  hCal_Presh  = new TH2F("hCal_Presh"," Calorimeter vs Preshower; Preshower ; Calorimeter ",100,0.,9000.,100,0.,6000.);                    HistList->Add(hCal_Presh);

 //-- Preshower
  hPresh_pulse  = new TH1F("hPresh_pulse"," Preshower pulse ",300,0.,300.);                                       HistList->Add(hPresh_pulse);
  hPresh_sum    = new TH1F("hPresh_sum"," Preshower Sum (GeV)",200,0.,9000.);                                     HistList->Add(hPresh_sum);
  hPresh_sum_el = new TH1F("hPresh_sum_el"," Preshower Sum for electrons",200,0.,9000.);                          HistList->Add(hPresh_sum_el);
  hPresh_sum_pi = new TH1F("hPresh_sum_pi"," Preshower Sum for pions",200,0.,9000.);                              HistList->Add(hPresh_sum_pi);

  //-- Cherenkov
  hCher_pulse  = new TH1F("hCher_pulse"," Cherenkov pulse ",300,0.,300.);                                        HistList->Add(hCher_pulse);
  hCher_sum    = new TH1F("hCher_sum"," Cherenkov Sum (GeV)",200,0.,2000.);                                      HistList->Add(hCher_sum);
  hCher_sum_el = new TH1F("hCher_sum_el"," Cherenkov Sum for electrons",200,0.,2000.);                           HistList->Add(hCher_sum_el);
  hCher_sum_pi = new TH1F("hCher_sum_pi"," Cherenkov Sum for pions",200,0.,2000.);                               HistList->Add(hCher_sum_pi);

  //-- Multiplicity
  hMult_pulse  = new TH1F("hMult_pulse"," Multiplicity pulse ",300,0.,300.);                                     HistList->Add(hMult_pulse);
  hMult_sum    = new TH1F("hMult_sum"," Multiplicity Sum (GeV)",200,0.,7000.);                                   HistList->Add(hMult_sum);
  hMult_sum_el = new TH1F("hMult_sum_el"," Multiplicity Sum for electrons",200,0.,7000.);                        HistList->Add(hMult_sum_el);
  hMult_sum_pi = new TH1F("hMult_sum_pi"," Multiplicity Sum for pions",200,0.,7000.);                            HistList->Add(hMult_sum_pi);
  
  //////
  f125_el = new TH1F("f125_el","GEM-TRD f125 Peak Amp for Electrons ; ADC Amplitude ; Counts ",100,0.,4096);                  HistList->Add(f125_el);
  f125_pi = new TH1F("f125_pi","GEM-TRD f125 Peak Amp for Pions ; ADC Amplitude ; Counts ",100,0.,4096);                  HistList->Add(f125_pi);
  f125_el_max_late = new TH1F("f125_el_max_late","GEM-TRD f125 Max Amp for Electrons ; ADC Amplitude ; Counts ",100,0.,4096);           HistList->Add(f125_el_max_late);
  f125_pi_max_late = new TH1F("f125_pi_max_late","GEM-TRD f125 Max Amp for Pions ; ADC Amplitude ; Counts ",100,0.,4096);           HistList->Add(f125_pi_max_late);
  f125_el_max_early = new TH1F("f125_el_max_early","GEM-TRD f125 Max Amp for Electrons ; ADC Amplitude ; Counts ",100,0.,4096);           HistList->Add(f125_el_max_early);
  f125_pi_max_early = new TH1F("f125_pi_max_early","GEM-TRD f125 Max Amp for Pions ; ADC Amplitude ; Counts ",100,0.,4096);           HistList->Add(f125_pi_max_early);
  f125_el_max = new TH1F("f125_el_max","GEM-TRD f125 Max Amp for Electrons ; ADC Amplitude ; Counts ",100,0.,4096);           HistList->Add(f125_el_max);
  f125_pi_max = new TH1F("f125_pi_max","GEM-TRD f125 Max Amp for Pions ; ADC Amplitude ; Counts ",100,0.,4096);           HistList->Add(f125_pi_max);
  mmg1_f125_el = new TH1F("mmg1_f125_el","MMG1-TRD f125 Peak Amp for Electrons ; ADC Amplitude ; Counts ",100,0.,4096);       HistList->Add(mmg1_f125_el);
  mmg1_f125_pi = new TH1F("mmg1_f125_pi","MMG1-TRD f125 Peak Amp for Pions ; ADC Amplitude ; Counts ",100,0.,4096);       HistList->Add(mmg1_f125_pi);
  mmg1_f125_el_max_late = new TH1F("mmg1_f125_el_max_late","MMG1-TRD f125 Max Amp for Electrons ; ADC Amplitude ; Counts ",100,0.,4096);           HistList->Add(mmg1_f125_el_max_late);
  mmg1_f125_pi_max_late = new TH1F("mmg1_f125_pi_max_late","MMG1-TRD f125 Max Amp for Pions ; ADC Amplitude ; Counts ",100,0.,4096);           HistList->Add(mmg1_f125_pi_max_late);
  mmg1_f125_el_max_early = new TH1F("mmg1_f125_el_max_early","MMG1-TRD f125 Max Amp for Electrons ; ADC Amplitude ; Counts ",100,0.,4096);           HistList->Add(mmg1_f125_el_max_early);
  mmg1_f125_pi_max_early = new TH1F("mmg1_f125_pi_max_early","MMG1-TRD f125 Max Amp for Pions ; ADC Amplitude ; Counts ",100,0.,4096);           HistList->Add(mmg1_f125_pi_max_early);
  mmg1_f125_el_max = new TH1F("mmg1_f125_el_max","MMG1-TRD f125 Max Amp for Electrons ; ADC Amplitude ; Counts ",100,0.,4096);           HistList->Add(mmg1_f125_el_max);
  mmg1_f125_pi_max = new TH1F("mmg1_f125_pi_max","MMG1-TRD f125 Max Amp for Pions ; ADC Amplitude ; Counts ",100,0.,4096);           HistList->Add(mmg1_f125_pi_max);
  gem_el_eff = new TH2F("gem_el_eff","GEM-TRD Electron Efficiency ; X (fADC) [mm] ; Y (SRS) [mm]",110,-55.,55.,110,-55.,55.);    HistList->Add(gem_el_eff);
  gem_pi_eff = new TH2F("gem_pi_eff","GEM-TRD Pion Efficiency ; X (fADC) [mm] ; Y (SRS) [mm]",110,-55.,55.,110,-55.,55.);    HistList->Add(gem_pi_eff);
  mmg1_el_eff = new TH2F("mmg1_el_eff","MMG1-TRD Electron Efficiency ; X (fADC) [mm] ; Y (SRS) [mm]",110,-55.,55.,110,-55.,55.);    HistList->Add(mmg1_el_eff);
  mmg1_pi_eff = new TH2F("mmg1_pi_eff","MMG1-TRD Pion Efficiency ; X (fADC) [mm] ; Y (SRS) [mm]",110,-55.,55.,110,-55.,55.);    HistList->Add(mmg1_pi_eff);
  gem_mmg1_doubleX = new TH2F("gem_mmg1_doubleX","TRD X Correlation (fADC) ; GEM-TRD X Chan ; MMG1-TRD X Chan ",256,-0.5,255.5,256,-0.5,255.5);    HistList->Add(gem_mmg1_doubleX);
  gem_mmg1_doubleY = new TH2F("gem_mmg1_doubleY","TRD Y Correlation (SRS) ; GEM-TRD Y Chan ; MMG1-TRD Y Chan ",256,-0.5,255.5,256,-0.5,255.5);    HistList->Add(gem_mmg1_doubleY);
  
  // --- SRS ---
  //--GEMTracker 1
  hgemtrkr_1_peak_xy = new TH2F("hgemtrkr_1_peak_xy","GEM-TRKR1 Peak X-Y Correlation (mm); Peak X [mm]; Peak Y [mm] ",110,-55.,55.,110,-55.,55.);    HistList->Add(hgemtrkr_1_peak_xy);
  hgemtrkr_1_atlas_xy = new TH2F("hgemtrkr_1_atlas_xy","GEM-TRKR1 X-Y Correlation (ATLAS TRIGGER); Peak X [mm]; Peak Y [mm] ",110,-55.,55.,110,-55.,55.);    HistList->Add(hgemtrkr_1_atlas_xy);
  hgemtrkr_1_peak_x = new TH1F("hgemtrkr_1_peak_x"," GEM-TRKR1 Peak X ; X [mm] ",110,-55.,55.);  HistList->Add(hgemtrkr_1_peak_x);
  hgemtrkr_1_peak_y = new TH1F("hgemtrkr_1_peak_y"," GEM-TRKR1 Peak Y ; Y [mm] ",110,-55.,55.);  HistList->Add(hgemtrkr_1_peak_y);
  hgemtrkr_1_peak_x_height = new TH1F("hgemtrkr_1_peak_x_height"," GEM-TRKR1 Peak Amplitudes in X; ADC Value ",100,0.,4096.);  HistList->Add(hgemtrkr_1_peak_x_height);
  hgemtrkr_1_peak_y_height = new TH1F("hgemtrkr_1_peak_y_height"," GEM-TRKR1 Peak Amplitudes in Y ; ADC Value ",100,0.,4096.);  HistList->Add(hgemtrkr_1_peak_y_height);
  //--GEMTracker 2
  hgemtrkr_2_peak_xy = new TH2F("hgemtrkr_2_peak_xy","GEM-TRKR2 Peak X-Y Correlation (mm); Peak X [mm]; Peak Y [mm] ",110,-55.,55.,110,-55.,55.);    HistList->Add(hgemtrkr_2_peak_xy);
  hgemtrkr_2_atlas_xy = new TH2F("hgemtrkr_2_atlas_xy","GEM-TRKR2 X-Y Correlation (ATLAS TRIGGER); Peak X [mm]; Peak Y [mm] ",110,-55.,55.,110,-55.,55.);    HistList->Add(hgemtrkr_2_atlas_xy);
  hgemtrkr_2_peak_x = new TH1F("hgemtrkr_2_peak_x"," GEM-TRKR2 Peak X ; X [mm] ",110,-55.,55.);  HistList->Add(hgemtrkr_2_peak_x);
  hgemtrkr_2_peak_y = new TH1F("hgemtrkr_2_peak_y"," GEM-TRKR2 Peak Y Pos ; Y [mm] ",110,-55.,55.);  HistList->Add(hgemtrkr_2_peak_y);
  hgemtrkr_2_peak_x_height = new TH1F("hgemtrkr_2_peak_x_height"," GEM-TRKR2 Peak Amplitudes in X; ADC Value ",100,0.,4096.);  HistList->Add(hgemtrkr_2_peak_x_height);
  hgemtrkr_2_peak_y_height = new TH1F("hgemtrkr_2_peak_y_height"," GEM-TRKR2 Peak Amplitudes in Y ; ADC Value ",100,0.,4096.);  HistList->Add(hgemtrkr_2_peak_y_height);
  //--GEMTracker 3
  hgemtrkr_3_peak_xy = new TH2F("hgemtrkr_3_peak_xy","GEM-TRKR3 Peak X-Y Correlation (mm); Peak X [mm]; Peak Y [mm] ",110,-55.,55.,110,-55.,55.);    HistList->Add(hgemtrkr_3_peak_xy);
  hgemtrkr_3_atlas_xy = new TH2F("hgemtrkr_3_atlas_xy","GEM-TRKR3 X-Y Correlation (ATLAS TRIGGER); Peak X [mm]; Peak Y [mm] ",110,-55.,55.,110,-55.,55.);    HistList->Add(hgemtrkr_3_atlas_xy);
  hgemtrkr_3_peak_x = new TH1F("hgemtrkr_3_peak_x"," GEM-TRKR3 Peak X Pos ; X [mm] ",110,-55.,55.);  HistList->Add(hgemtrkr_3_peak_x);
  hgemtrkr_3_peak_y = new TH1F("hgemtrkr_3_peak_y"," GEM-TRKR3 Peak Y Pos ; Y [mm] ",110,-55.,55.);  HistList->Add(hgemtrkr_3_peak_y);
  hgemtrkr_3_peak_x_height = new TH1F("hgemtrkr_3_peak_x_height"," GEM-TRKR3 Peak Amplitudes in X; ADC Value ",100,0.,4096.);  HistList->Add(hgemtrkr_3_peak_x_height);
  hgemtrkr_3_peak_y_height = new TH1F("hgemtrkr_3_peak_y_height"," GEM-TRKR3 Peak Amplitudes in Y ; ADC Value ",100,0.,4096.);  HistList->Add(hgemtrkr_3_peak_y_height);
  
  mmg1_peak_y = new TH1F("mmg1_peak_y"," MMG1-TRD Peak Y Pos (SRS) ; Y [mm] ",110,-55.,55.);  HistList->Add(mmg1_peak_y);
  hmmg1_peak_y_height = new TH1F("hmmg1_peak_y_height"," MMG1-TRD Peak Amplitudes in Y ; ADC Value ",100,0.,4096.);  HistList->Add(hmmg1_peak_y_height);
  gem_peak_y = new TH1F("gem_peak_y"," GEM-TRD Peak Y Pos (SRS) ; Y [mm] ",110,-55.,55.);  HistList->Add(gem_peak_y);
  hgem_peak_y_height = new TH1F("hgem_peak_y_height"," GEM-TRD Peak Amplitudes in Y ; ADC Value ",100,0.,4096.);  HistList->Add(hgem_peak_y_height);
  
  hgemtrkr_1_gem = new TH2F("hgemtrkr_1_gem","GEM-TRKR1 * GEM Y Correlation; GEM-TRD Y [mm]; GEM-TRKR1 Y [mm] ",110,-55.,55.,110,-55.,55.);    HistList->Add(hgemtrkr_1_gem);
  hgemtrkr_1_mmg1 = new TH2F("hgemtrkr_1_mmg1","GEM-TRKR1 & MMG1 Y Correlation; MMG1-TRD Y [mm]; GEM-TRKR1 Y [mm] ",110,-55.,55.,110,-55.,55.);    HistList->Add(hgemtrkr_1_mmg1);
  
  //--External Tracking
  TH1F *f125_tracker_hits = new TH1F("f125_tracker_hits","GEM-TRD Track Extr. Hits; X Chan [mm]",100.,-0.5,99.5);   HistList->Add(f125_tracker_hits);
  TH1F *f125_tracker_nohits = new TH1F("f125_tracker_nohits","",100.,-0.5,99.5);   HistList->Add(f125_tracker_nohits);
  TH1F *f125_tracker_eff = new TH1F("f125_tracker_eff","GEM-TRD Track Extr. Eff; X Chan [mm]",100.,-0.5,99.5);   HistList->Add(f125_tracker_eff);
  TH1F *f125_tracker_eff2d = new TH1F("f125_tracker_eff2d","",100.,-0.5,99.5);   HistList->Add(f125_tracker_eff2d);
  TH1F *gem_residuals = new TH1F("gem_residuals","GEM-TRD Residual Hits; X Chan [mm]",100.,-0.5,99.5);   HistList->Add(gem_residuals);
  TH1F *mmg1_f125_tracker_hits = new TH1F("mmg1_f125_tracker_hits","MMG1-TRD Track Extr. Hits; X Chan [mm]",100.,-0.5,99.5);   HistList->Add(mmg1_f125_tracker_hits);
  TH1F *mmg1_f125_tracker_nohits = new TH1F("mmg1_f125_tracker_nohits","",100.,-0.5,99.5);   HistList->Add(mmg1_f125_tracker_nohits);
  TH1F *mmg1_f125_tracker_eff = new TH1F("mmg1_f125_tracker_eff","MMG1-TRD Track Extr. Hits; X Chan [mm]",100.,-0.5,99.5);   HistList->Add(mmg1_f125_tracker_eff);
  TH1F *mmg1_f125_tracker_eff2d = new TH1F("mmg1_f125_tracker_eff2d","",100.,-0.5,99.5);   HistList->Add(mmg1_f125_tracker_eff2d);
  TH1F *mmg1_residuals = new TH1F("mmg1_residuals","MMG1-TRD Residual Hits; X Chan [mm]",100.,-0.5,99.5);   HistList->Add(mmg1_residuals);
  
  //---- GEM-TRD --
  float GEM_THRESH=140.; //175
  float MMG1_THRESH=125.; //150

  f125_el_evt_display = new TH2F("f125_el_evt_display","GEM-TRD track for Electrons ; Time Response (8ns) ; Channel ",100,100.5,200.5,200,-0.5,249.5);  HistList->Add(f125_el_evt_display);
  f125_el_evt_display->SetStats(0);  f125_el_evt_display->SetMinimum(GEM_THRESH);  f125_el_evt_display->SetMaximum(1000.);
  f125_el_raw = new TH2F("f125_el_raw","GEM-TRD raw for Electrons ; Time Response (8ns) ; Channel ",100,100.5,200.5,200,-0.5,249.5);    HistList->Add(f125_el_raw);
  f125_el_raw->SetStats(0);  f125_el_raw->SetMinimum(GEM_THRESH);   f125_el_raw->SetMaximum(1000.);
  f125_pi_raw = new TH2F("f125_pi_raw","GEM-TRD raw for Pions ; Time Response (8ns) ; Channel ",100,100.5,200.5,200,-0.5,249.5);    HistList->Add(f125_pi_raw);
  f125_pi_raw->SetStats(0);  f125_pi_raw->SetMinimum(GEM_THRESH);   f125_pi_raw->SetMaximum(1000.);
  
  //f125_el_fit = new TH2F("f125_el_fit","GEM-TRD track for Electrons ; Time Response (8ns) ; Channel ",250,0.5,250.5,240,0.5,240.5);     HistList->Add(f125_el_fit);
  f125_amp2ds = new TH2F("f125_amp2ds","GEM-TRD ADC Amp in Time (EXT. TRK) ; Time Response (8ns) ; X Channel ",200,0.5,200.5,240,0.5,240.5);           HistList->Add(f125_amp2ds);
  f125_el_amp2d = new TH2F("f125_el_amp2d","GEM-TRD ADC Amp in Time for Electrons ; Time Response (8ns) ; X Channel ",200,0.5,200.5,240,0.5,240.5);           HistList->Add(f125_el_amp2d);
  f125_pi_amp2d = new TH2F("f125_pi_amp2d","GEM-TRD ADC Amp in Time for Pions ; Time Response (8ns) ; X Channel ",200,0.5,200.5,240,0.5,240.5);           HistList->Add(f125_pi_amp2d);
  mmg1_f125_amp2ds = new TH2F("mmg1_f125_amp2ds","MMG1-TRD ADC Amp in Time (EXT. TRK) ; Time Response (8ns) ; X Channel ",200,0.5,200.5,240,0.5,240.5);           HistList->Add(mmg1_f125_amp2ds);
  mmg1_f125_el_amp2d = new TH2F("mmg1_f125_el_amp2d","MMG1-TRD ADC Amp in Time for Electrons ; Time Response (8ns) ; X Channel ",200,0.5,200.5,240,0.5,240.5);   HistList->Add(mmg1_f125_el_amp2d);
  mmg1_f125_pi_amp2d = new TH2F("mmg1_f125_pi_amp2d","MMG1-TRD ADC Amp in Time for Pions ; Time Response (8ns) ; X Channel ",200,0.5,200.5,240,0.5,240.5);   HistList->Add(mmg1_f125_pi_amp2d);
  
  f125_el_clu2d = new TH2F("f125_el_clu2d","GEM-TRD Amp for Electrons (Clusters)",200,0.5,200.5,240,0.5,240.5);               HistList->Add(f125_el_clu2d);
  f125_pi_clu2d = new TH2F("f125_pi_clu2d","GEM-TRD Amp for Pions (Clusters)",200,0.5,200.5,240,0.5,240.5);               HistList->Add(f125_pi_clu2d);
  mmg1_f125_el_clu2d = new TH2F("mmg1_f125_el_clu2d","MMG1-TRD Amp for Electrons (Clusters)",200,0.5,200.5,240,0.5,240.5);    HistList->Add(mmg1_f125_el_clu2d);
  mmg1_f125_pi_clu2d = new TH2F("mmg1_f125_pi_clu2d","MMG1-TRD Amp for Pions (Clusters)",200,0.5,200.5,240,0.5,240.5);    HistList->Add(mmg1_f125_pi_clu2d);
  
  //=========================================
  gem_zHist = new  TH1F("gem_zHist", "gem_zHist", 20, 80., 200.);
  mmg1_zHist = new  TH1F("mmg1_zHist", "mmg1_zHist", 20, 80., 200.);
  
  TFile* fHits;
  #ifdef SAVE_TRACK_HITS
    char hitsFileName[256]; sprintf(hitsFileName, "RootOutput/cern24/trd_singleTrackHits_Run_%06d.root", RunNum);
    fHits = new TFile(hitsFileName, "RECREATE");
    //-- GEM-TRD
    EVENT_VECT_GEM = new TTree("gem_hits","GEM TTree with single track hit info");
    EVENT_VECT_GEM->Branch("event_num",&event_num,"event_num/I");
    EVENT_VECT_GEM->Branch("nhit",&gem_nhit,"gem_nhit/I");
    EVENT_VECT_GEM->Branch("ecal_energy",&Ecal_Energy,"Ecal_Energy/float");
    EVENT_VECT_GEM->Branch("presh_energy",&Presh_Energy,"Presh_Energy/float");
    EVENT_VECT_GEM->Branch("mult_energy",&Mult_Energy,"Mult_Energy/float");
    EVENT_VECT_GEM->Branch("xpos",&gem_xpos);
    EVENT_VECT_GEM->Branch("zpos",&gem_zpos);
    EVENT_VECT_GEM->Branch("dedx",&gem_dedx);
    EVENT_VECT_GEM->Branch("parID",&gem_parID);
    EVENT_VECT_GEM->Branch("zHist",&gem_zHist_vect);
    EVENT_VECT_GEM->Branch("xposc",&clu_xpos);
    EVENT_VECT_GEM->Branch("zposc",&clu_zpos);
    EVENT_VECT_GEM->Branch("dedxc",&clu_dedx);
    EVENT_VECT_GEM->Branch("widthc",&clu_width);
    EVENT_VECT_GEM->Branch("xch_max",&gem_xch_max);
    //EVENT_VECT_GEM->Branch("xchtrkr",&gemtrkr_xch_max);
    //EVENT_VECT_GEM->Branch("ychtrkr",&gemtrkr_ych_max);
    EVENT_VECT_GEM->Branch("amp_max",&gem_amp_max);
    EVENT_VECT_GEM->Branch("time_max",&gem_time_max);
    //EVENT_VECT_GEM->Branch("xamptrkr",&gemtrkr_xamp_max);
    //EVENT_VECT_GEM->Branch("yamptrkr",&gemtrkr_yamp_max);
    EVENT_VECT_GEM->Branch("chi2",&gem_chi2cc);
    EVENT_VECT_GEM->Branch("Fint",&gem_integral);
    EVENT_VECT_GEM->Branch("a0",&a0);
    EVENT_VECT_GEM->Branch("a1",&a1);
    //-- MMG1-TRD
    EVENT_VECT_MMG1 = new TTree("mmg1_hits","MMG1 TTree with single track hit info");
    EVENT_VECT_MMG1->Branch("event_num",&event_num,"event_num/I");
    EVENT_VECT_MMG1->Branch("nhit",&mmg1_nhit,"mmg1_nhit/I");
    EVENT_VECT_MMG1->Branch("xpos",&mmg1_xpos);
    EVENT_VECT_MMG1->Branch("zpos",&mmg1_zpos);
    EVENT_VECT_MMG1->Branch("dedx",&mmg1_dedx);
    EVENT_VECT_MMG1->Branch("parID",&mmg1_parID);
    EVENT_VECT_MMG1->Branch("zHist",&mmg1_zHist_vect);
    EVENT_VECT_MMG1->Branch("xposc",&mmg1_clu_xpos);
    EVENT_VECT_MMG1->Branch("zposc",&mmg1_clu_zpos);
    EVENT_VECT_MMG1->Branch("dedxc",&mmg1_clu_dedx);
    EVENT_VECT_MMG1->Branch("widthc",&mmg1_clu_width);
    EVENT_VECT_MMG1->Branch("xch_max",&mmg1_xch_max);
    EVENT_VECT_MMG1->Branch("time_max",&mmg1_amp_max);
    EVENT_VECT_MMG1->Branch("amp_max",&mmg1_time_max);
    EVENT_VECT_MMG1->Branch("chi2",&mmg1_chi2cc);
    EVENT_VECT_MMG1->Branch("Fint",&mmg1_integral);
    EVENT_VECT_MMG1->Branch("a0",&mmg1_a0);
    EVENT_VECT_MMG1->Branch("a1",&mmg1_a1);
    
  #endif
  
  TStopwatch timer;
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  if (MaxEvt>0) nentries=MaxEvt;  //-- limit number of events for test
  
  //==================================================================================================
  //                      E v e n t    L o o p
  //==================================================================================================
  printf("***>>>  Begin Event Loop - 1st evt=%lld, Last evt=%lld \n",FirstEvt,MaxEvt);
  timer.Start();
  Long64_t jentry=0;
  int el_count=0, pi_count=0, atlas_trigger_count=0;
  int ShowEvent=0;
  
  for (jentry=FirstEvt; jentry<nentries; jentry++) { //-- Event Loop --
    Count("EVT");
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);
    nbytes += nb;
    if (!(jentry%NPRT))
      printf("------- evt=%llu  f125_raw_count=%llu f125_pulse_count=%llu f250_wraw_count=%llu, srs_peak_count=%llu \n",jentry,f125_wraw_count, f125_pulse_count, f250_wraw_count, gem_peak_count);
    event_num = jentry;
    gem_nhit=0;
    mmg1_nhit=0;
    
    //-- GEM-TRD
    gem_xpos.clear();
    gem_zpos.clear();
    gem_dedx.clear();
    gem_parID.clear();
    gem_zHist->Reset();
    gem_zHist_vect.clear();
    clu_xpos.clear();
    clu_zpos.clear();
    clu_dedx.clear();
    clu_width.clear();
    gem_amp_max.clear();
    gem_xch_max.clear();
    gem_time_max.clear();
    gem_chi2cc.clear();
    gem_integral.clear();
    
    //-- MMG1-TRD
    mmg1_xpos.clear();
    mmg1_zpos.clear();
    mmg1_dedx.clear();
    mmg1_parID.clear();
    mmg1_zHist->Reset();
    mmg1_zHist_vect.clear();
    mmg1_clu_xpos.clear();
    mmg1_clu_zpos.clear();
    mmg1_clu_dedx.clear();
    mmg1_clu_width.clear();
    mmg1_amp_max.clear();
    mmg1_xch_max.clear();
    mmg1_time_max.clear();
    mmg1_chi2cc.clear();
    mmg1_integral.clear();
    
    //==================================================================================================
    //                    Process Fa250  Pulse data
    //==================================================================================================
    
    #ifdef USE_250_PULSE
      
      double cal_energy=-1., cher_energy=-1., presh_energy=-1., mult_counter_energy=-1.;
      bool electron_tag=false, pion_tag=false, atlas_trigger=false;
      for (ULong64_t i=0;i<f250_pulse_count; i++) {
      #ifdef VERBOSE
        printf("F250:: i=%d  sl=%d, ch=%d,  npk=%d  time=%d amp=%d ped=%f \n",i,f250_pulse_slot->at(i),f250_pulse_channel->at(i),f250_pulse_pulse_number->at(i),f250_pulse_course_time->at(i),f250_pulse_pulse_peak->at(i),f250_pulse_pedestal->at(i)/4.);
      #endif
        if (f250_pulse_channel->at(i)==0) {cal_energy=f250_pulse_pulse_peak->at(i)-f250_pulse_pedestal->at(i)/4.;}
        if (f250_pulse_channel->at(i)==1) {presh_energy=f250_pulse_pulse_peak->at(i)-f250_pulse_pedestal->at(i)/4.;}
        if (f250_pulse_channel->at(i)==2) {cher_energy=f250_pulse_pulse_peak->at(i)-f250_pulse_pedestal->at(i)/4.;}
        if (f250_pulse_channel->at(i)==3) {if ((f250_pulse_pulse_peak->at(i)-f250_pulse_pedestal->at(i)/4.)>=1000.) {electron_tag=true;  Count("ATLel");  }}
        if (f250_pulse_channel->at(i)==4) {if ((f250_pulse_pulse_peak->at(i)-f250_pulse_pedestal->at(i)/4.)>=1000.) {pion_tag=true;      Count("ATLpi");  }}
        if (f250_pulse_channel->at(i)==5) {if ((f250_pulse_pulse_peak->at(i)-f250_pulse_pedestal->at(i)/4.)>=1000.) {atlas_trigger=true; Count("ATLtrig");}}
        if (f250_pulse_channel->at(i)==7) {mult_counter_energy=f250_pulse_pulse_peak->at(i)-f250_pulse_pedestal->at(i)/4.;}
      }
      #ifdef VERBOSE
        printf("F250: calorimeter_en=%f, preshower_en=%f, cherenkov_en=%f,  mult_counter_en=%f  el_tag=%d pi_tag=%d atlas_trig=%d \n", cal_energy, presh_energy, cher_energy, mult_counter_energy, electron_tag, pion_tag, atlas_trigger);
      #endif
      if (electron_tag) el_count++;
      if (pion_tag) pi_count++;
      if (atlas_trigger) atlas_trigger_count++;
      
    #endif
    
    //==================================================================================================
    //                    Process Fa250  RAW data
    //==================================================================================================
    
    h250_size->Fill(f250_wraw_count);
    double CalSum=0;
    bool electron=false;
    bool pion=false;
    int ch_cal     = 0;
    int ch_presh   = 1;
    int ch_cher    = 2;
    int ch_el_tag  = 3;
    int ch_pi_tag  = 4;
    int ch_atl_trig= 5;
    int ch_mult    = 7;
    double calIntegral=0;
    double preshIntegral=0;
    double cherIntegral=0;
    double multIntegral=0;
    
    for (ULong64_t i=0; i<f250_wraw_count; i++) { // --- fadc250 channels loop
      #ifdef VERBOSE
        if (jentry<MAX_PRINT) printf("F250:: i=%lld  sl=%d, ch=%d, idx=%d, cnt=%d \n",i,f250_wraw_slot->at(i),f250_wraw_channel->at(i),f250_wraw_samples_index->at(i),f250_wraw_samples_count->at(i));
      #endif
      int fadc_chan = f250_wraw_channel->at(i);
      int fadc_window = f250_wraw_samples_count->at(i);
      hCal_occ->Fill(fadc_chan+0.);
      int amax250=0;
      int tmax250=9;
      if (fadc_chan==ch_cal) hCal_pulse->Reset();
      if (fadc_chan==ch_presh) hPresh_pulse->Reset();
      if (fadc_chan==ch_cher) hCher_pulse->Reset();
      if (fadc_chan==ch_mult) hMult_pulse->Reset();
      
      int ped=100; //CHANGE LATER...
      for (int si=0; si<fadc_window; si++) {
        int adc = f250_wraw_samples->at(f250_wraw_samples_index->at(i)+si);
        if (fadc_chan==ch_cal)   hCal_pulse->Fill(si+0.5,adc+0.);
        if (fadc_chan==ch_cher)  hCher_pulse->Fill(si+0.5,adc+0.);
        if (fadc_chan==ch_presh) hPresh_pulse->Fill(si+0.5,adc+0.);
        if (fadc_chan==ch_mult)  hMult_pulse->Fill(si+0.5,adc+0.);
        if (37 < si && si < 47) calIntegral+=adc-ped;  // -- Cal window
	        if (fadc_chan==ch_presh && 42 < si && si < 48) preshIntegral+=adc-ped;  // -- Presh window
	        if (fadc_chan==ch_cher  && 50 < si && si < 53) cherIntegral+=adc-ped;  // -- Cherenkov window
	        if (fadc_chan==ch_mult  && 44 < si && si < 48) multIntegral+=adc-ped;  // -- Multiplicity window
        if (adc>amax250) {
          amax250=adc;
          tmax250=si;
        }
      } // --  end of fadc window samples loop --
      if (fadc_chan==ch_cal) { CalSum=calIntegral; }
    } // -- end of fadc250 channels loop --
    
    Ecal_Energy=CalSum;
    Presh_Energy=preshIntegral;
    Mult_Energy=multIntegral;
    hCal_Presh->Fill(preshIntegral,CalSum);
    hCal_Cher->Fill(cherIntegral,CalSum);
    hCal_sum->Fill(CalSum); hPresh_sum->Fill(preshIntegral);  hCher_sum->Fill(cherIntegral);  hMult_sum->Fill(multIntegral);
    if(electron_tag)  {  hCal_sum_el->Fill(CalSum); hPresh_sum_el->Fill(preshIntegral);  hCher_sum_el->Fill(cherIntegral);  hMult_sum_el->Fill(multIntegral); }
    if(pion_tag)      {  hCal_sum_pi->Fill(CalSum); hPresh_sum_pi->Fill(preshIntegral);  hCher_sum_pi->Fill(cherIntegral);  hMult_sum_pi->Fill(multIntegral); }
    
    if (!atlas_trigger) continue;
    if (!electron_tag && !pion_tag) continue;
    
    //==================================================================================================
    //                    Process SRS data
    //==================================================================================================
    
    //===========================================================
    //  GEMTracker (SRS) Correlation with TRD Prototypes
    //===========================================================
    
    ULong64_t gt_1_idx_x = 0, gt_1_idx_y=0, gt_2_idx_x = 0, gt_2_idx_y = 0, gt_3_idx_x = 0, gt_3_idx_y = 0, mmg1_idx_y=0, gem_idx_y=0;
    
    double gemtrkr1_xamp_max=-1., gemtrkr1_xch_max=-1.;
    double gemtrkr1_yamp_max=-1., gemtrkr1_ych_max=-1.;
    double gemtrkr2_xamp_max=-1., gemtrkr2_xch_max=-1.;
    double gemtrkr2_yamp_max=-1., gemtrkr2_ych_max=-1.;
    double gemtrkr3_xamp_max=-1., gemtrkr3_xch_max=-1.;
    double gemtrkr3_yamp_max=-1., gemtrkr3_ych_max=-1.;
    
    double gemtrkr_1_peak_pos_y[gem_peak_count];
    double gemtrkr_1_peak_pos_x[gem_peak_count];
    double gemtrkr_1_peak_x_height[gem_peak_count];
    double gemtrkr_1_peak_y_height[gem_peak_count];
    double gemtrkr_1_peak_ych[gem_peak_count];
    double gemtrkr_1_peak_xch[gem_peak_count];
    double gemtrkr_2_peak_pos_y[gem_peak_count];
    double gemtrkr_2_peak_pos_x[gem_peak_count];
    double gemtrkr_2_peak_x_height[gem_peak_count];
    double gemtrkr_2_peak_y_height[gem_peak_count];
    double gemtrkr_2_peak_ych[gem_peak_count];
    double gemtrkr_2_peak_xch[gem_peak_count];
    double gemtrkr_3_peak_pos_y[gem_peak_count];
    double gemtrkr_3_peak_pos_x[gem_peak_count];
    double gemtrkr_3_peak_x_height[gem_peak_count];
    double gemtrkr_3_peak_y_height[gem_peak_count];
    double gemtrkr_3_peak_ych[gem_peak_count];
    double gemtrkr_3_peak_xch[gem_peak_count];
    double mmg1_peak_pos_y[gem_peak_count];
    double mmg1_peak_y_height[gem_peak_count];
    double gem_peak_pos_y[gem_peak_count];
    double gem_peak_y_height[gem_peak_count];
    
    for (ULong64_t i=0; i<gem_peak_count; i++) {
      gemtrkr_1_peak_pos_y[gem_peak_count] = -1000;
      gemtrkr_1_peak_pos_x[gem_peak_count] = -1000;
      gemtrkr_1_peak_x_height[gem_peak_count] = -1000;
      gemtrkr_1_peak_y_height[gem_peak_count] = -1000;
      gemtrkr_1_peak_ych[gem_peak_count] = -1000;
      gemtrkr_1_peak_xch[gem_peak_count] = -1000;
      gemtrkr_2_peak_pos_y[gem_peak_count] = -1000;
      gemtrkr_2_peak_pos_x[gem_peak_count] = -1000;
      gemtrkr_2_peak_x_height[gem_peak_count] = -1000;
      gemtrkr_2_peak_y_height[gem_peak_count] = -1000;
      gemtrkr_2_peak_ych[gem_peak_count] = -1000;
      gemtrkr_2_peak_xch[gem_peak_count] = -1000;
      gemtrkr_3_peak_pos_y[gem_peak_count] = -1000;
      gemtrkr_3_peak_pos_x[gem_peak_count] = -1000;
      gemtrkr_3_peak_x_height[gem_peak_count] = -1000;
      gemtrkr_3_peak_y_height[gem_peak_count] = -1000;
      gemtrkr_3_peak_ych[gem_peak_count] = -1000;
      gemtrkr_3_peak_xch[gem_peak_count] = -1000;
      mmg1_peak_pos_y[gem_peak_count] = -1000;
      mmg1_peak_y_height[gem_peak_count] = -1000;
      gem_peak_pos_y[gem_peak_count] = -1000;
      gem_peak_y_height[gem_peak_count] = -1000;
    }
   
    for (ULong64_t i=0; i<gem_peak_count; i++) { //-- SRS Peaks Loop
      
      if (gem_peak_plane_name->at(i) == "GEMTR1X") {
        gemtrkr_1_peak_pos_x[gt_1_idx_x] = gem_peak_real_pos->at(i);
        if (gemtrkr_1_peak_pos_x[gt_1_idx_x]<0) gemtrkr_1_peak_pos_x[gt_1_idx_x]+=50.; else if (gemtrkr_1_peak_pos_x[gt_1_idx_x]>0) gemtrkr_1_peak_pos_x[gt_1_idx_x]-=50.;  gemtrkr_1_peak_pos_x[gt_1_idx_x]*=-1.;
        gemtrkr_1_peak_x_height[gt_1_idx_x] = gem_peak_height->at(i);
        gemtrkr_1_peak_xch[gt_1_idx_x] = gem_peak_index->at(i);
        gt_1_idx_x++; Count("gt1_x");
      } if (gem_peak_plane_name->at(i) == "GEMTR1Y") {
          gemtrkr_1_peak_pos_y[gt_1_idx_y] = gem_peak_real_pos->at(i);
          if (gemtrkr_1_peak_pos_y[gt_1_idx_y]<0) gemtrkr_1_peak_pos_y[gt_1_idx_y]+=50.; else if (gemtrkr_1_peak_pos_y[gt_1_idx_y]>0) gemtrkr_1_peak_pos_y[gt_1_idx_y]-=50.;  gemtrkr_1_peak_pos_y[gt_1_idx_y]*=-1.;
          gemtrkr_1_peak_y_height[gt_1_idx_y] = gem_peak_height->at(i);
          gemtrkr_1_peak_ych[gt_1_idx_y] = gem_peak_index->at(i);
          gt_1_idx_y++; Count("gt1_y");
      } if (gem_peak_plane_name->at(i) == "GEMTR2X") {
        gemtrkr_2_peak_pos_x[gt_2_idx_x] = gem_peak_real_pos->at(i);
        if (gemtrkr_2_peak_pos_x[gt_2_idx_x]<0) gemtrkr_2_peak_pos_x[gt_2_idx_x]+=50.; else if (gemtrkr_2_peak_pos_x[gt_2_idx_x]>0) gemtrkr_2_peak_pos_x[gt_2_idx_x]-=50.;  gemtrkr_2_peak_pos_x[gt_2_idx_x]*=-1.;
        gemtrkr_2_peak_x_height[gt_2_idx_x] = gem_peak_height->at(i);
        gemtrkr_2_peak_xch[gt_2_idx_x] = gem_peak_index->at(i);
        gt_2_idx_x++; Count("gt2_x");
      } if (gem_peak_plane_name->at(i) == "GEMTR2Y") {
          gemtrkr_2_peak_pos_y[gt_2_idx_y] = gem_peak_real_pos->at(i);
          if (gemtrkr_2_peak_pos_y[gt_2_idx_y]<0) gemtrkr_2_peak_pos_y[gt_2_idx_y]+=50.; else if (gemtrkr_2_peak_pos_y[gt_2_idx_y]>0) gemtrkr_2_peak_pos_y[gt_2_idx_y]-=50.;  gemtrkr_2_peak_pos_y[gt_2_idx_y]*=-1.;
          gemtrkr_2_peak_y_height[gt_2_idx_y] = gem_peak_height->at(i);
          gemtrkr_2_peak_ych[gt_2_idx_y] = gem_peak_index->at(i);
          gt_2_idx_y++; Count("gt2_y");
      } if (gem_peak_plane_name->at(i) == "GEMTR3X") {
        gemtrkr_3_peak_pos_x[gt_3_idx_x] = gem_peak_real_pos->at(i);
        if (gemtrkr_3_peak_pos_x[gt_3_idx_x]<0) gemtrkr_3_peak_pos_x[gt_3_idx_x]+=50.; else if (gemtrkr_3_peak_pos_x[gt_3_idx_x]>0) gemtrkr_3_peak_pos_x[gt_3_idx_x]-=50.;  gemtrkr_3_peak_pos_x[gt_3_idx_x]*=-1.;
        gemtrkr_3_peak_x_height[gt_3_idx_x] = gem_peak_height->at(i);
        gemtrkr_3_peak_xch[gt_3_idx_x] = gem_peak_index->at(i);
        gt_3_idx_x++; Count("gt3_x");
      } if (gem_peak_plane_name->at(i) == "GEMTR3Y") {
          gemtrkr_3_peak_pos_y[gt_3_idx_y] = gem_peak_real_pos->at(i);
          if (gemtrkr_3_peak_pos_y[gt_3_idx_y]<0) gemtrkr_3_peak_pos_y[gt_3_idx_y]+=50.; else if (gemtrkr_3_peak_pos_y[gt_3_idx_y]>0) gemtrkr_3_peak_pos_y[gt_3_idx_y]-=50.;  gemtrkr_3_peak_pos_y[gt_3_idx_y]*=-1.;
          gemtrkr_3_peak_y_height[gt_3_idx_y] = gem_peak_height->at(i);
          gemtrkr_3_peak_ych[gt_3_idx_y] = gem_peak_index->at(i);
          gt_3_idx_y++; Count("gt3_y");
      } if (gem_peak_plane_name->at(i) == "MMG1TRDY") {
        mmg1_peak_pos_y[mmg1_idx_y] = gem_peak_real_pos->at(i);
        if (mmg1_peak_pos_y[mmg1_idx_y]<0) mmg1_peak_pos_y[mmg1_idx_y]+=50.; else if (mmg1_peak_pos_y[mmg1_idx_y]>0) mmg1_peak_pos_y[mmg1_idx_y]-=50.;  mmg1_peak_pos_y[mmg1_idx_y]*=-1.;
        mmg1_peak_y_height[mmg1_idx_y] = gem_peak_height->at(i);
        mmg1_idx_y++; Count("mmg1_y");
      } if (gem_peak_plane_name->at(i) == "VU_GEMTRDY") {
        gem_peak_pos_y[gem_idx_y] = gem_peak_real_pos->at(i);
        if (gem_peak_pos_y[gem_idx_y]<0) gem_peak_pos_y[gem_idx_y]+=50.; else if (gem_peak_pos_y[gem_idx_y]>0) gem_peak_pos_y[gem_idx_y]-=50.;  gem_peak_pos_y[gem_idx_y]*=-1.;
        gem_peak_y_height[gem_idx_y] = gem_peak_height->at(i);
        gem_idx_y++; Count("gem_y");
      }
    }
    
    for (ULong64_t j=0; j<gt_1_idx_y; j++) {
      hgemtrkr_1_peak_y->Fill(gemtrkr_1_peak_pos_y[j]);
      hgemtrkr_1_peak_y_height->Fill(gemtrkr_1_peak_y_height[j]);
      if (gemtrkr_1_peak_y_height[j]>gemtrkr1_yamp_max) {
        gemtrkr1_yamp_max=gemtrkr_1_peak_y_height[j];
        gemtrkr1_ych_max=gemtrkr_1_peak_ych[j];
      }
      for (ULong64_t k=0; k<gt_1_idx_x; k++) {
          hgemtrkr_1_peak_xy->Fill(gemtrkr_1_peak_pos_x[k], gemtrkr_1_peak_pos_y[j]);
          if (atlas_trigger) hgemtrkr_1_atlas_xy->Fill(gemtrkr_1_peak_pos_x[k], gemtrkr_1_peak_pos_y[j]);
      }
    }
    for (ULong64_t j=0; j<gt_2_idx_y; j++) {
      hgemtrkr_2_peak_y->Fill(gemtrkr_2_peak_pos_y[j]);
      hgemtrkr_2_peak_y_height->Fill(gemtrkr_2_peak_y_height[j]);
      if (gemtrkr_2_peak_y_height[j]>gemtrkr2_yamp_max) {
        gemtrkr2_yamp_max=gemtrkr_2_peak_y_height[j];
        gemtrkr2_ych_max=gemtrkr_2_peak_ych[j];
      }
        for (ULong64_t k=0; k<gt_2_idx_x; k++) {
          hgemtrkr_2_peak_xy->Fill(gemtrkr_2_peak_pos_x[k], gemtrkr_2_peak_pos_y[j]);
          if (atlas_trigger) hgemtrkr_2_atlas_xy->Fill(gemtrkr_2_peak_pos_x[k], gemtrkr_2_peak_pos_y[j]);
      }
    }
    for (ULong64_t j=0; j<gt_3_idx_y; j++) {
      hgemtrkr_3_peak_y->Fill(gemtrkr_3_peak_pos_y[j]);
      hgemtrkr_3_peak_y_height->Fill(gemtrkr_3_peak_y_height[j]);
      if (gemtrkr_3_peak_y_height[j]>gemtrkr3_yamp_max) {
        gemtrkr3_yamp_max=gemtrkr_3_peak_y_height[j];
        gemtrkr3_ych_max=gemtrkr_3_peak_ych[j];
      }
        for (ULong64_t k=0; k<gt_3_idx_x; k++) {
          hgemtrkr_3_peak_xy->Fill(gemtrkr_3_peak_pos_x[k], gemtrkr_3_peak_pos_y[j]);
          if (atlas_trigger) hgemtrkr_3_atlas_xy->Fill(gemtrkr_3_peak_pos_x[k], gemtrkr_3_peak_pos_y[j]);
      }
    }
    for (ULong64_t k=0; k<gt_1_idx_x; k++) {
      hgemtrkr_1_peak_x->Fill(gemtrkr_1_peak_pos_x[k]);
      hgemtrkr_1_peak_x_height->Fill(gemtrkr_1_peak_x_height[k]);
      if (gemtrkr_1_peak_x_height[k]>gemtrkr1_xamp_max) {
        gemtrkr1_xamp_max=gemtrkr_1_peak_x_height[k];
        gemtrkr1_xch_max=gemtrkr_1_peak_xch[k];
      }
    }
    for (ULong64_t k=0; k<gt_2_idx_x; k++) {
      hgemtrkr_2_peak_x->Fill(gemtrkr_2_peak_pos_x[k]);
      hgemtrkr_2_peak_x_height->Fill(gemtrkr_2_peak_x_height[k]);
      if (gemtrkr_2_peak_x_height[k]>gemtrkr2_xamp_max) {
        gemtrkr2_xamp_max=gemtrkr_2_peak_x_height[k];
        gemtrkr2_xch_max=gemtrkr_2_peak_xch[k];
      }
    }
    for (ULong64_t k=0; k<gt_3_idx_x; k++) {
      hgemtrkr_3_peak_x->Fill(gemtrkr_3_peak_pos_x[k]);
      hgemtrkr_3_peak_x_height->Fill(gemtrkr_3_peak_x_height[k]);
      if (gemtrkr_3_peak_x_height[k]>gemtrkr3_xamp_max) {
        gemtrkr3_xamp_max=gemtrkr_3_peak_x_height[k];
        gemtrkr3_xch_max=gemtrkr_3_peak_xch[k];
      }
    }
    for (ULong64_t k=0; k<mmg1_idx_y; k++) {
      mmg1_peak_y->Fill(mmg1_peak_pos_y[k]);
      hmmg1_peak_y_height->Fill(mmg1_peak_y_height[k]);
      for (ULong64_t j=0; j<gt_1_idx_y; j++) {
        hgemtrkr_1_mmg1->Fill(mmg1_peak_pos_y[k], gemtrkr_1_peak_pos_y[j]);
      }
    }
    for (ULong64_t k=0; k<gem_idx_y; k++) {
      gem_peak_y->Fill(gem_peak_pos_y[k]);
      hgem_peak_y_height->Fill(gem_peak_y_height[k]);
      for (ULong64_t j; j<mmg1_idx_y; j++) {
        gem_mmg1_doubleY->Fill(gem_peak_pos_y[k], mmg1_peak_pos_y[j]);
      }
      for (ULong64_t j=0; j<gt_1_idx_y; j++) {
        hgemtrkr_1_gem->Fill(gem_peak_pos_y[k], gemtrkr_1_peak_pos_y[j]);
      }
    }
    //=============== END GEMTracker (SRS) Correlations ==============
    
    //==================================================================================================
    //                    Process Fa125  Pulse  data
    //==================================================================================================
    
    if (gt_1_idx_x>0 && gt_2_idx_x>0) { //--External tracking condition
      Count("trk_hit");
      float x1=gemtrkr1_xch_max*0.4;
      float x2=gemtrkr2_xch_max*0.4;
      //float x3=gemtrkr3_xch_max*0.4;
      float a=(x2-x1)/(z2-z1);
      float b=((x1)*z2-(x2)*z1)/(z2-z1);
      float gem_extr = a*zgem+b;
      float mmg1_extr = a*zmmg1+b;
      f125_tracker_hits->Fill(gem_extr);
      mmg1_f125_tracker_hits->Fill(mmg1_extr);
      bool match = false, match_mmg1 = false;
      
      for (ULong64_t i=0; i<f125_pulse_count; i++) {
        
        float peak_amp = f125_pulse_peak_amp->at(i);
        float ped = f125_pulse_pedestal->at(i);
        if (0 > ped || ped > 200 ) ped = 100;
        float amp = peak_amp-ped;
        if (amp<0) amp=0;
        float time = f125_pulse_peak_time->at(i);
        int fADCSlot = f125_pulse_slot->at(i);
        int fADCChan = f125_pulse_channel->at(i);
        int gemChan = GetGEMChan(fADCChan, fADCSlot);
        int mmg1Chan = GetMMG1Chan(fADCChan, fADCSlot, RunNum);
        
        if (gemChan>-1 && amp>GEM_THRESH) {
          float xgem = gemChan*0.4;
          gem_residuals->Fill(xgem-gem_extr);
          if (abs(xgem-gem_extr)<30) { //within 30 mm
            if (!match) {
              f125_tracker_eff->Fill(gem_extr);
              match = true;
            }
            f125_amp2ds->Fill(time,gemChan,amp);
          }
        }
        if (mmg1Chan>-1 && amp>MMG1_THRESH) {
          float xmmg1 = mmg1Chan*0.4;
          mmg1_residuals->Fill(xmmg1-mmg1_extr);
          if (abs(xmmg1-mmg1_extr)<30) { //within 30 mm
            if (!match_mmg1) {
              mmg1_f125_tracker_eff->Fill(mmg1_extr);
              match_mmg1 = true;
            }
            mmg1_f125_amp2ds->Fill(time,mmg1Chan,amp);
          }
        }
      } //---- End Fadc125 Pulse Loop ----
    } //-- End external tracker condition
    
    #if 1
      //==================================================================================================
      //                    Again Process Fa125  Pulse  data
      //==================================================================================================
      
      //--------- Event Histogram Filling ---------
      f125_fit->Reset();
      mmg1_f125_fit->Reset();
      double chi2cc_gem=-999., chi2cc_mmg1=-999.;
      double integral_gem=0., integral_mmg1=0.;
      hevtk->Reset();
      hevtck->Reset();
      double gem_ampmax=-1., mmg1_ampmax=-1., gem_xchmax=-1, mmg1_xchmax=-1;
      int gem_timemax=0, mmg1_timemax=0;
      double xaver=0, xaver2=0;
      int naver=0;
      float ped_max=0;
      
      for (ULong64_t i=0; i<f125_pulse_count; i++) {
        
      	float peak_amp = f125_pulse_peak_amp->at(i);
      	float ped = f125_pulse_pedestal->at(i);
      	if (0 > ped || ped > 200 ) ped = 100;
      	float amp=peak_amp-ped;
      	float time=f125_pulse_peak_time->at(i);
      	int fADCSlot = f125_pulse_slot->at(i);
      	int fADCChan = f125_pulse_channel->at(i);
      	
      	int gemChan = GetGEMChan(fADCChan, fADCSlot);
      	int mmg1Chan = GetMMG1Chan(fADCChan, fADCSlot, RunNum);
        
      	if (amp<0) amp=0;
      	
      	if (gemChan>-1 && amp>GEM_THRESH) {
          //for (ULong64_t j=0; j<f125_pulse_count; j++) {
          if (mmg1Chan>-1 && amp>MMG1_THRESH) {
            gem_mmg1_doubleX->Fill(gemChan,mmg1Chan,1.);
          //}
          }
          
          f125_fit->Fill(time,gemChan,amp);
          //------------ for Clusters -------------
          int TimeWindowStart = 45;
          int time0=int(time)-TimeWindowStart;
          if ( 0 < time0 && time0 < 130 ) { // --- drop early and late hits ---
            hevtck->SetBinContent(time0,gemChan,amp/10);
            hevtk->SetBinContent(time0,gemChan,amp/10);
          }
    	    if (gem_ampmax<amp) {
            gem_ampmax=amp;
            gem_timemax=time;
            gem_xchmax=gemChan;
          }
          if (electron_tag) {
      	    f125_el->Fill(amp);
            f125_el_amp2d->Fill(time,gemChan,amp);
      	    f125_el_clu2d->Fill(time,gemChan,1.);
      	    gem_xpos.push_back(gemChan);
      	    gem_dedx.push_back(amp);
      	    gem_zpos.push_back(time);
      	    gem_parID.push_back(1);
      	    gem_nhit++;
      	    gem_zHist->Fill(time, amp);
          } else if (pion_tag) {
            f125_pi->Fill(amp);
            f125_pi_amp2d->Fill(time,gemChan,amp);
            f125_pi_clu2d->Fill(time,gemChan,1.);
            gem_xpos.push_back(gemChan);
            gem_dedx.push_back(amp);
            gem_zpos.push_back(time);
            gem_parID.push_back(0);
            gem_nhit++;
            //gem_zHist->Fill(time, amp);
          }
    	  }
    	  if (amp>MMG1_THRESH && mmg1Chan>-1) {
          mmg1_f125_fit->Fill(time,mmg1Chan,amp);
    	    if (mmg1_ampmax<amp) {
            mmg1_ampmax=amp;
            mmg1_timemax=time;
            mmg1_xchmax=mmg1Chan;
          }
          if (electron_tag) {
      	    mmg1_f125_el_amp2d->Fill(time,mmg1Chan,amp);
      	    mmg1_f125_el->Fill(amp);
      	    mmg1_f125_el_clu2d->Fill(time,mmg1Chan,1.);
      	    mmg1_xpos.push_back(mmg1Chan);
      	    mmg1_dedx.push_back(amp);
      	    mmg1_zpos.push_back(time);
      	    mmg1_parID.push_back(1);
      	    mmg1_nhit++;
      	    mmg1_zHist->Fill(time, amp);
          } else if (pion_tag) {
            mmg1_f125_pi_amp2d->Fill(time,mmg1Chan,amp);
            mmg1_f125_pi->Fill(amp);
            mmg1_f125_pi_clu2d->Fill(time,mmg1Chan,1.);
            mmg1_xpos.push_back(mmg1Chan);
            mmg1_dedx.push_back(amp);
            mmg1_zpos.push_back(time);
            mmg1_parID.push_back(0);
            mmg1_nhit++;
            //mmg1_zHist->Fill(time, amp);
          }
    	  }
    	} //--- end Fa125 Pulse Loop ---
      
      //===================================================================
      //                    Chi^2 Fit Calculation
      //===================================================================
      
      char f125Title[80]; sprintf(f125Title,"GEM Display Event: %lld   Run: %d; z pos,mm; y pos,mm ",jentry,RunNum);
      f125_fit->SetTitle(f125Title);
      if (f125_fit->GetEntries()!=0) {
        std::pair<Double_t, Double_t> fitResult  = TrkFit(f125_fit,fx1,"fx1",1);
        chi2cc_gem = fitResult.first;
        integral_gem = fitResult.second;
      }
      double a0 = fx1.GetParameter(0);
      double a1 = fx1.GetParameter(1);
  
      char mmg1f125Title[80]; sprintf(mmg1f125Title,"MMG1 Display Event: %lld   Run: %d; z pos,mm; y pos,mm ",jentry,RunNum);
      mmg1_f125_fit->SetTitle(mmg1f125Title);
      if (mmg1_f125_fit->GetEntries()!=0) {
        std::pair<Double_t, Double_t> fitResult  = TrkFit(mmg1_f125_fit,fx2,"fx2",1);
        chi2cc_mmg1 = fitResult.first;
        integral_mmg1 = fitResult.second;
      }
      double a0_mmg1 = fx2.GetParameter(0);
      double a1_mmg1 = fx2.GetParameter(1);
      
      //=================== Loop to  calculate average and RMS  ================
      if (naver>0) xaver=xaver/naver;
      #ifdef VERBOSE
        if (fabs(ped_max-100.) > 50. && electron_tag ) {
          printf("Ev=%lld, CalSum=%f xaver2=%f amp_max=%f, ped_max=%f,  el=%d pi=%d \n",jentry,CalSum,xaver2,gem_amp_max,ped_max,electron_tag,pion_tag);
        }
      #endif
      for (ULong64_t i=0;i<gem_xpos.size();i++) {
        if ( 110 < gem_xpos.at(i) && gem_xpos.at(i) < 130 && 65 < gem_zpos.at(i) && gem_zpos.at(i) <135 )  xaver2+=((gem_xpos.at(i)-xaver)*(gem_xpos.at(i)-xaver));
      }
      if (naver>0) xaver2=sqrt(xaver2/naver);
      if (electron_tag) {
        aver2d_e->Fill(xaver,xaver2);
      } else if (pion_tag) {
        aver2d_p->Fill(xaver,xaver2);
      }
      #ifdef VERBOSE
        if ( xaver2<2. && electron_tag ) {
          printf("xaver2=%f el=%d pi=%d \n",xaver2,electron_tag,pion_tag);
        }
      #endif
      
      //==================== Max Amplitude histos ============================
      if (electron_tag==1) {
        if (gem_ampmax>0) {
          f125_el_max->Fill(gem_ampmax);
          if (gem_timemax>130) f125_el_max_late->Fill(gem_ampmax);
          if (gem_timemax<85) f125_el_max_early->Fill(gem_ampmax);
        }
        if (mmg1_ampmax>0) {
          mmg1_f125_el_max->Fill(mmg1_ampmax);
          if (mmg1_timemax>140) mmg1_f125_el_max_late->Fill(mmg1_ampmax);
          if (mmg1_timemax<90) mmg1_f125_el_max_early->Fill(mmg1_ampmax);
        }
      }
      if (pion_tag==1) {
        if (gem_ampmax>0) {
          f125_pi_max->Fill(gem_ampmax);
          if (gem_timemax>130) f125_pi_max_late->Fill(gem_ampmax);
          if (gem_timemax<85) f125_pi_max_early->Fill(gem_ampmax);
        }
        if (mmg1_ampmax>0) {
          mmg1_f125_pi_max->Fill(mmg1_ampmax);
          if (mmg1_timemax>140) mmg1_f125_pi_max_late->Fill(mmg1_ampmax);
          if (mmg1_timemax<90) mmg1_f125_pi_max_early->Fill(mmg1_ampmax);
        }
      }
      gem_amp_max.push_back(gem_ampmax);
      gem_time_max.push_back(gem_timemax);
      gem_xch_max.push_back(gem_xchmax);
      gem_chi2cc.push_back(chi2cc_gem);
      gem_integral.push_back(integral_gem);
      mmg1_amp_max.push_back(mmg1_ampmax);
      mmg1_time_max.push_back(mmg1_timemax);
      mmg1_xch_max.push_back(mmg1_xchmax);
      mmg1_chi2cc.push_back(chi2cc_mmg1);
      mmg1_integral.push_back(integral_mmg1);
      
      for (int i=1; i<21; i++) {
        gem_zHist_vect.push_back(gem_zHist->GetBinContent(i));
        mmg1_zHist_vect.push_back(mmg1_zHist->GetBinContent(i));
      }
    #endif
    //======================= End Process Fa125 Pulse data ================================
    
    //==================================================================================================
    //                    Process Fa125  RAW data
    //==================================================================================================
    
    #ifdef USE_125_RAW
      #ifdef VERBOSE
        if (jentry<MAX_PRINT) printf("------------------ Fadc125  wraw_count = %llu ---------\n", f125_wraw_count);
      #endif
      mhevt->Reset();
      mhevtf->Reset();
      //mhevti->Reset();
      hevt->Reset();
      hevtc->Reset();
      //hevti->Reset();
      hevtf->Reset();
      
      for (ULong64_t i=0;i<f125_wraw_count; i++) { // --- fadc125 channels loop
        
        int fadc_window = f125_wraw_samples_count->at(i);
        int fADCSlot = f125_wraw_slot->at(i);
        int fADCChan = f125_wraw_channel->at(i);
        int gemChan = GetGEMChan(fADCChan, fADCSlot);
        int mmg1Chan = GetMMG1Chan(fADCChan, fADCSlot, RunNum);
        int amax=0;
        int tmax=0;
        //if (gemChan<0) continue;
        double DEDX_THR = GEM_THRESH;
        int TimeWindowStart = 45;
        int TimeMin = 0;
        int TimeMax = 130;
        
        for (int si=0; si<fadc_window; si++) {
          int time=si;
          int adc = f125_wraw_samples->at(f125_wraw_samples_index->at(i)+si);
          if (adc>4090) printf("!!!!!!!!!!!!!!!!!!!!!! ADC 125 overflow: %d \n",adc);
          if (adc>amax) {
            amax=adc;
            tmax=si;
          }
          if (adc>DEDX_THR) {
            if (gemChan>-1) {
              if (electron_tag)  {
                f125_el_raw->Fill(time,gemChan,adc);
              } else if (pion_tag) {
                f125_pi_raw->Fill(time,gemChan,adc);
              }
              time-=TimeWindowStart;
              if ( TimeMin > time || time > TimeMax ) continue; // --- drop early and late hits ---
              
              hevtc->SetBinContent(time,gemChan,adc/100.);
              hevt->SetBinContent(time,gemChan,adc/100.);
            }
            if (mmg1Chan>-1) {
              mhevt->SetBinContent(time-35,mmg1Chan,adc/100.);
            }
          }
        } // --  end of samples loop
      } // -- end of fadc125 raw channels loop
      #ifdef SHOW_EVT_DISPLAY
          #if (USE_PULSE>0)
            c2->cd(1);   hevtk->Draw("box");
          #else
            c2->cd(1); hevt->Draw("colz");
            c2->cd(6); mhevt->Draw("colz");
          #endif
          c2->cd(2);   hevtf->Draw("text");
          c2->cd(7);   mhevtf->Draw("text");
          c2->Modified(); c2->Update();
      #endif
      
      //==================================================================================================
      //            Begin NN Clustering & Track Fitting
      //==================================================================================================
      #if (USE_CLUST>0)
        // -------------------------------   hist dist clustering         ------------------------
        float clust_Xpos[MAX_CLUST];
        float clust_Ypos[MAX_CLUST];
        float clust_Zpos[MAX_CLUST];
        float clust_dEdx[MAX_CLUST];
        float clust_Size[MAX_CLUST];
        float clust_Width[MAX_CLUST][3];  // y1, y2, dy ; strips
        float clust_Length[MAX_CLUST][3]; // x1, x2, dx ; time
        float hits_Xpos[500];
        float hits_Ypos[500];
        float hits_Zpos[500];
        float hits_dEdx[500];
        float hits_Size[MAX_CLUST];
        float hits_Width[MAX_CLUST];  // y1, y2, dy ; strips
        float hits_Length[MAX_CLUST]; // x1, x2, dx ; time
        
        for (int k=0; k<MAX_CLUST; k++) {
          clust_Xpos[k]=0; clust_Ypos[k]=0; clust_Zpos[k]=0; clust_dEdx[k]=0;  clust_Size[k]=0;
          clust_Width[k][0]=999999;   	clust_Width[k][1]=-999999;   	clust_Width[k][2]=0;
          clust_Length[k][0]=999999;  	clust_Length[k][1]=-999999;  	clust_Length[k][2]=0;
        }
        int nclust=0;
        #if (USE_PULSE>0)
          TH2F* hp = hevtk; // -- hevt and hevtc should be same bin size
          TH2F* hpc = hevtck;
        #else
          TH2F* hp = hevt; // -- hevt and hevtc should be same bin size
          TH2F* hpc = hevtc;
        #endif
        int nx=hp->GetNbinsX();    int ny=hp->GetNbinsY();
        double xmi=hp->GetXaxis()->GetBinLowEdge(1);     double xma=hp->GetXaxis()->GetBinUpEdge(nx);
        double ymi=hp->GetYaxis()->GetBinLowEdge(1);     double yma=hp->GetYaxis()->GetBinUpEdge(ny);
        double binx = (xma-xmi)/nx;      double biny = (yma-ymi)/ny;
        #ifdef VERBOSE
          printf("nx=%d,ny=%d,xmi=%f,xma=%f,ymi=%f,yma=%f\n",nx,ny,xmi,xma,ymi,yma);
        #endif
        #if (USE_PULSE>0)
          float CL_DIST=3.3; // mm
          double THR2 = 0.01;
        #else
          float CL_DIST=2.9; // mm
          double THR2 = 0.2;
        #endif
        
        for (int ix=0; ix<nx; ix++) {  //-------------------- clustering loop ------------------------------------
          for (int iy=0; iy<ny; iy++) {
            double c1 = hpc->GetBinContent(ix,iy);                    // energy
            double x1=double(ix)/double(nx)*(xma-xmi)+xmi-binx/2.;    // drift time
            double y1=double(iy)/double(ny)*(yma-ymi)+ymi-biny/2.;    // X strip
            if (c1<THR2) continue;
            if (nclust==0) {
              clust_Xpos[nclust]=y1; clust_Ypos[nclust]=0; clust_Zpos[nclust]=x1;  clust_dEdx[nclust]=c1;  clust_Size[nclust]=1;
              clust_Width[nclust][0]=y1;   	clust_Width[nclust][1]=y1;   	clust_Width[nclust][2]=0;
              clust_Length[nclust][0]=x1;  	clust_Length[nclust][1]=x1;  	clust_Length[nclust][2]=0;
              nclust++; continue;
            }
            int added=0;
            for (int k=0; k<nclust; k++) {
              double dist=sqrt(pow((y1-clust_Xpos[k]),2.)+pow((x1-clust_Zpos[k]),2.)); //--- dist hit to clusters
              if (dist<CL_DIST) {
                clust_Xpos[k]=(y1*c1+clust_Xpos[k]*clust_dEdx[k])/(c1+clust_dEdx[k]);  //--  new X pos
                clust_Zpos[k]=(x1*c1+clust_Zpos[k]*clust_dEdx[k])/(c1+clust_dEdx[k]);  //--  new Z pos
                clust_dEdx[k]=c1+clust_dEdx[k];  // new dEdx
                clust_Size[k]=1+clust_Size[k];  // clust size in pixels
                if (y1<clust_Width[k][0]) clust_Width[k][0]=y1; if (y1>clust_Width[k][1]) clust_Width[k][1]=y1; clust_Width[k][2]=clust_Width[k][1]-clust_Width[k][0];
                if (x1<clust_Length[k][0]) clust_Length[k][0]=x1;if (x1>clust_Length[k][1]) clust_Length[k][1]=x1;clust_Length[k][2]=clust_Length[k][1]-clust_Length[k][0];
                hpc->SetBinContent(ix,iy,k+1.);
                added=1; break;
              }
            }
            if (added==0) {
              if (nclust+1>=MAX_CLUST) continue;
              clust_Xpos[nclust]=y1; clust_Ypos[nclust]=0; clust_Zpos[nclust]=x1;  clust_dEdx[nclust]=c1;  clust_Size[nclust]=1;
              clust_Width[nclust][0]=y1;   	clust_Width[nclust][1]=y1;   	clust_Width[nclust][2]=0;
              clust_Length[nclust][0]=x1;  	clust_Length[nclust][1]=x1;  	clust_Length[nclust][2]=0;
              nclust++;
            }
          }
        } //----------------------------------- end  clustering loop -----------------------------------------------
        
        #if (USE_PULSE>0)
          int MinClustSize=1;
          double MinClustWidth=0.00;
          double MinClustLength=0.0;
          double MaxClustLength=5.;
          double zStart =  0.; // mm
          double zEnd   = 29.; // mm
        #else
          int MinClustSize=1;//5;
          double MinClustWidth=0.001;
          double MinClustLength=0.01;
          double MaxClustLength=5.;
          double zStart =  0.; // mm
          double zEnd   = 30.; // mm
          double dEmin  = 2.; //
        #endif
        
        int ii=0;
        #ifdef VERBOSE
          printf("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
          printf("                Xpos   Ypos   Zpos       E    Width  Length   Size \n");
        #endif
        for (int k=0; k<nclust; k++) {
          #ifdef VERBOSE
            if (k<30) printf("%2d Clust(%2d): %6.1f %6.1f %6.1f %8.1f %6.2f %6.2f %8.1f  ",k,k+1,clust_Xpos[k],clust_Ypos[k],clust_Zpos[k],clust_dEdx[k],clust_Width[k][2],clust_Length[k][2],clust_Size[k]);
          #endif
          //-------------  Cluster Filter -----------------
          if ((clust_Size[k] >= MinClustSize && zStart < clust_Zpos[k] && clust_Zpos[k] < zEnd && clust_Width[k][2]>MinClustWidth) || clust_Length[k][2]<MaxClustLength ) {
        	  hits_Xpos[ii]=clust_Xpos[k];
        	  hits_Ypos[ii]=clust_Ypos[k];
        	  hits_Zpos[ii]=clust_Zpos[k];
        	  hits_dEdx[ii]=clust_dEdx[k];
            hits_Width[ii]=clust_Width[k][2];
            hits_Length[ii]=clust_Length[k][2];
        	  ii++;
            #ifdef VERBOSE
              if (k<30) printf("\n");
            #endif
        	  } else {
            #ifdef VERBOSE
              if (k<30) printf(" <--- skip \n");
            #endif
            }
        }
        int nhits=ii;
        // ----------------------- end hist dist clustering ---------------------------------
        
        //=================================== Draw HITS and CLUST  ============================================
        #ifdef SHOW_EVTbyEVT
          char hevtTitle[80]; sprintf(hevtTitle,"GEM Display Event: %lld   Run: %d e=%d #pi=%d; z pos,mm; y pos,mm ",jentry,RunNum,electron_tag,pion_tag);
          hevt->SetTitle(hevtTitle);
          hevtk->SetTitle(hevtTitle);
          char mhevtTitle[80]; sprintf(mhevtTitle,"MMG1 Display Event: %lld   Run: %d e=%d #pi=%d; z pos,mm; y pos,mm ",jentry,RunNum,electron_tag,pion_tag);
          mhevt->SetTitle(mhevtTitle);
          #ifdef VERBOSE
            printf("hits_SIZE=%d  Clust size = %d \n",nhits,nclust);
          #endif
          if (jentry<NPRT || ShowEvent>0 ) {
          	c2->cd(1); gPad->Modified(); gPad->Update();
          	int COLMAP[]={1,2,3,4,6,5};
          	int pmt=22 ,pmt0 = 20; // PM type
            int max2draw = nclust; //std::min(nclust, 33); // draw max 33 clust
            for(int i = 0; i < max2draw; i++){
          	  TMarker m = TMarker(clust_Zpos[i],clust_Xpos[i],pmt);
          	  int tcol=2; //min(tracks[i],6);
          	  if (clust_Size[i]<MinClustSize) pmt=22; else pmt=pmt0;
          	  int mcol = COLMAP[tcol-1];   m.SetMarkerColor(mcol);   m.SetMarkerStyle(pmt);
          	  m.SetMarkerSize(0.7+clust_dEdx[i]/300);
          	  m.DrawClone();  //gPad->Modified(); gPad->Update();
          	}
            gPad->Modified(); gPad->Update();
          }
        #endif
        
        #if (USE_GNN==1)   // GNN MC
          //----------------------------------------------------------
          //--   Send to Model simulation
          //----------------------------------------------------------
          
          #ifdef VERBOSE
            printf("**> Start Model simulation nclust=%d nhits=%d \n",nclust,nhits);
          #endif
          std::vector<int> tracks(nhits, 0);
          std::vector<float> Xcl;
          std::vector<float> Zcl;
          Xcl.clear();
          Zcl.clear();
          for (int n=0; n<nhits; n++) {
          	Xcl.push_back(hits_Xpos[n]);
          	Zcl.push_back(hits_Zpos[n]);
          }
          doPattern(Xcl, Zcl, tracks);  //---- call GNN ---
          #ifdef VERBOSE
            printf("**> End Model simulation \n"); //===================================================
          #endif
          #ifdef SHOW_EVTbyEVT
            if (jentry<NPRT || ShowEvent>0 ) {
              c2->cd(2); gPad->Modified(); gPad->Update();
              int COLMAP[]={1,2,3,4,6,5};
              for(ULong64_t i = 0; i < tracks.size(); i++) {
                #ifdef VERBOSE
                  if (i<30) printf("i=%d trk=%d |  %8.2f,%8.2f\n",i, tracks[i], Xcl[i], Zcl[i]);
              	#endif
                TMarker m = TMarker(hits_Zpos[i],hits_Xpos[i],24);
              	int tcol=min(tracks[i],6);
              	int mcol = COLMAP[tcol-1];   m.SetMarkerColor(mcol);   m.SetMarkerStyle(41);     m.SetMarkerSize(1.5);
            	  m.DrawClone();  gPad->Modified(); gPad->Update();
              }
              #ifdef VERBOSE
                printf("\n\n");
                printf("**> End Cluster Plot \n");
              #endif
            }
          #endif
          //--------------------------------------------------
          //----           Track fitting                 -----
          //--------------------------------------------------
          
          #ifdef VERBOSE
            printf("==> GNN: tracks sort  : trk_siz=%ld \r\n", tracks.size());
          #endif
          //-----------------   tracks sorting -------------
          std::vector<std::vector<float>> TRACKS;  // -- [Nnodes][N_params]
          TRACKS.resize(nhits);
          std::vector<float> hit_coord(2,0);
          std::vector<int>  TRACKS_N(nhits, 0); // [Nnodes]
          for (int i = 0; i < nhits; i++)  { TRACKS_N[i] = 0;  }
          std::vector<float> xz(2,0);
          for (int i2 = 0; i2 < nhits; i2++) {
          	int num =  tracks[i2];
          	int num2 = std::max(0, std::min(num, nhits - 1));
            #ifdef VERBOSE
              if (i2<20) printf("==> lstm3:track sort i=%d  : num=%d(%d) x=%f z=%f \n", i2, num, num2,  Xcl[i2],Zcl[i2]);
          	#endif
            xz[0]=Xcl[i2]; 	xz[1]=Zcl[i2];
          	TRACKS[num2].push_back(Xcl[i2]); TRACKS[num2].push_back(Zcl[i2]);
          	TRACKS_N[num2]++;
          }
          
          #if (DEBUG > 1)
            for (int i2 = 0; i2 < nhits; i2++) {
              printf(" trdID=%d n_hits=%d v_size=%d \n",i2,TRACKS_N[i2],TRACKS[i2].size());
              for (ULong64_t i3 = 0; i3 < TRACKS[i2].size(); i3+=2) {
                printf(" trkID=%d  hit=%d x=%f z=%f \n",i2,i3/2,TRACKS[i2].at(i3),TRACKS[i2].at(i3+1));
              }
              if ( TRACKS_N[i2]>0) printf("\n");
            }
          #endif
          //--- end tracks sorting ---
          
          #if (USE_FIT==1)
            //-----------------------------------
            //---       linear fitting        ---
            //-----------------------------------
            static TMultiGraph *mg;
            if (mg != NULL ) delete mg;
            mg = new TMultiGraph();
            mg->SetTitle(" ML-FPGA response; z pos,mm; y pos,mm ");
            int NTRACKS=0;
            int MIN_HITS=2;
            Double_t p0,p1;
            
            for (int i2 = 1; i2 < nhits; i2++) {  // tracks loop; zero track -> noise
              
             if (TRACKS_N[i2]<MIN_HITS) continue;   //---- select 2 (x,z) and more hits on track ----
            	#ifdef VERBOSE
                printf("==> fit: start trk: %d \r\n", i2);
              #endif
            	std::vector<Double_t> x;
            	std::vector<Double_t> y;
            	for (int i3 = 0; i3 < (int)TRACKS[i2].size(); i3+=2) {
            	  #ifdef VERBOSE
                  printf(" trkID=%d  hit=%d x=%f z=%f \n",i2,i3/2,TRACKS[i2].at(i3),TRACKS[i2].at(i3+1));
            	  #endif
                x.push_back(TRACKS[i2].at(i3+1));
            	  y.push_back(TRACKS[i2].at(i3));
            	}
              #ifdef SHOW_EVTbyEVT
              	gErrorIgnoreLevel = kBreak; // Suppress warning messages from empty chi^2 fit data
              	TGraph *g = new TGraph(TRACKS_N[i2], &x[0], &y[0]);  g->SetMarkerStyle(21); g->SetMarkerColor(i2);
              	TF1 *f = new TF1("f", "[1] * x + [0]");
              	g->Fit(f,"Q");
                //  --- get fit parameters ---
                TF1 *ffunc=g->GetFunction("f");
                p0=ffunc->GetParameter(0);
                p1=ffunc->GetParameter(1);
                Double_t chi2x_nn = ffunc->GetChisquare();
                Double_t Ndfx_nn = ffunc->GetNDF();
                double chi2nn=chi2x_nn/Ndfx_nn;
                #ifdef VERBOSE
                  printf("+++++>>  Track = %d fit: p0=%f p1=%f (%f deg) ff(15)=%f chi2nn=%f \n",i2,p0,p1,p1/3.1415*180.,ffunc->Eval(15.),chi2nn);
                #endif
            	  mg->Add(g,"p");
              #endif
              NTRACKS++;
            }  //  end tracks loop
            #ifdef SHOW_EVTbyEVT
              if (NTRACKS!=1) continue;  // --- skip event ----
              if (nhits<3)    continue;  // --- skip event ----
              if (jentry<NPRT || ShowEvent>0 ) {
                c2->cd(3); mg->Draw("AP");
                mg->GetXaxis()->SetLimits(Xmin,Xmax);
                mg->SetMinimum(Ymin);
                mg->SetMaximum(Ymax);
                gPad->Modified(); gPad->Update();
              }
            #endif
          #endif // USE_FIT
        #endif // USE_GNN MC
        
        //******************************************************************************
        #ifdef SHOW_EVTbyEVT
            cout<<"Event#="<<event_num<<" Electron="<<electron_tag<<"  Pion="<<pion_tag<<"  AtlasTrig="<<atlas_trigger<<" CherenkovEn="<<cher_energy<<" CalorimeterEn="<<cal_energy<<" PreshowerEn="<<presh_energy<<" CounterEn="<<mult_counter_energy<<" #ofTracks="<<NTRACKS<<endl;
            #ifdef WRITE_CSV
              WriteToCSV(csvFile,event_num,electron_tag,pion_tag,atlas_trigger,cher_energy,cal_energy,presh_energy,mult_counter_energy,NTRACKS,chi2cc_gem);
            #endif
            c3->cd(1);  hCal_sum->Draw();         gPad->Modified();   gPad->Update();
            c3->cd(2);  hCal_pulse->Draw("hist"); gPad->Modified();   gPad->Update();
            c3->cd(3);  hCher_pulse->Draw("hist"); gPad->Modified();  gPad->Update();
            c3->cd(4);  hPresh_pulse->Draw("hist"); gPad->Modified(); gPad->Update();
            c3->cd(5);  hMult_pulse->Draw("hist"); gPad->Modified(); gPad->Update();
            
            c2->cd(4);  f125_fit->Draw("box");    gPad->Modified(); gPad->Update();
            c2->cd(5);  hevt->Draw("colz");       gPad->Modified();   gPad->Update();
            c2->cd(9);  mmg1_f125_fit->Draw("box"); gPad->Modified(); gPad->Update();
            c2->cd(10);  mhevt->Draw("colz");       gPad->Modified(); gPad->Update();
            printf(" all done, click middle pad ... a0=%f a1=%f (%f deg)  fx1(150)=%f chi2cc_gem=%f  \n",a0,a1,a1/3.1415*180.,fx1.Eval(150.),chi2cc_gem);
            if (electron_tag || pion_tag) c2->cd(1); gPad->WaitPrimitive();
          ShowEvent=0;
        #endif
      #endif   // (USE_CLUST)
    #endif   //=======================  End Fa125 RAW process Loop  =====================================
    
    if (NTRACKS==1) Count("singleTRK");
    if (NTRACKS==1 && electron_tag) Count("snTRKel");
    if (NTRACKS==1 && pion_tag) Count("snTRKpi");
    
    //=====================================================================================
    //===                Fill Root TTree Hits                                            ===
    //=====================================================================================
    
    #ifdef SAVE_TRACK_HITS
      if (gem_nhit>0) EVENT_VECT_GEM->Fill();
      if (mmg1_nhit>0)EVENT_VECT_MMG1->Fill();
      for (int n=0; n<nhits; n++) {
        clu_xpos.push_back(hits_Xpos[n]);
        clu_zpos.push_back(hits_Zpos[n]);
        clu_dedx.push_back(hits_dEdx[n]);
        clu_width.push_back(hits_Width[n]);
        //clu_length.push_back(hits_Length[n]);
        /*
        if (electron_tag) {
          f125_el_clu2d->Fill(hits_Zpos[n],hits_Xpos[n],hits_dEdx[n]);
        }
        if (pion_tag) {
          f125_pi_clu2d->Fill(hits_Zpos[n],hits_Xpos[n],hits_dEdx[n]);
        }
        */
      }
    #endif
  } // ------------------------ END of event loop  ------------------------------
  
  timer.Stop();
  cout<<"***>>> End Event Loop, Elapsed Time:"<<endl; timer.Print();
  #ifdef WRITE_CSV
    csvFile.close();
  #endif
  cout<<" Total events = "<<jentry<<endl;
  cout<<" # of electrons="<<el_count<<" # of pions="<<pi_count<<" # of atlas triggers="<<atlas_trigger_count<<endl;
  
  //---Drift time distribution plot ---
  TH1D *f125_drift = f125_el_amp2d->ProjectionX("f125_drift",110,140);
  TH1D *f125_drift_c = new TH1D(*f125_drift);
  double gemDriftScale = 1./f125_drift_c->GetEntries();
  f125_drift_c->Scale(gemDriftScale);
  TH1D *mmg1_f125_drift = mmg1_f125_el_amp2d->ProjectionX("mmg1_f125_drift",110,140);
  TH1D *mmg1_f125_drift_c = new TH1D(*mmg1_f125_drift);
  double mmg1DriftScale = 1./mmg1_f125_drift_c->GetEntries();
  mmg1_f125_drift_c->Scale(mmg1DriftScale);
  f125_drift_c->SetLineColor(4);  HistList->Add(f125_drift_c);
  mmg1_f125_drift_c->SetLineColor(2); HistList->Add(mmg1_f125_drift_c);
  f125_drift_c->GetXaxis()->SetTitle("Time Response (8ns)");
  f125_drift_c->SetTitle("Drift Time Distribution");
  TLegend *l1 = new TLegend(0.75,0.65,0.9,0.9);
  l1->SetNColumns(2);
  l1->AddEntry(f125_drift_c, "GEM", "l");
  l1->AddEntry(mmg1_f125_drift_c, "MMG1", "l");
  
  TCanvas *c4 = new TCanvas("c4","Drift Time Distribution", 1200, 800);
  c4->cd();
  f125_drift_c->Draw();
  mmg1_f125_drift_c->Draw("same");
  l1->Draw();
  
  //=====================================================================================
  //===                 S A V E   H I S T O G R A M S                                ====
  //=====================================================================================
  TFile* fOut;
  char rootFileName[256]; sprintf(rootFileName, "RootOutput/cern24/Run_%06d_Output.root", RunNum);
  fOut = new TFile(rootFileName, "RECREATE");
  fOut->cd();
  cout<<"Writing Output File: "<<rootFileName<<endl;
  HistList->Write("HistDQM", TObject::kSingleKey);
  c4->Write();
  fOut->Close();
  delete fOut;
  
  //=====================================================================================
  //===                 S A V E   T R A C K   H I T   T T R E E                      ====
  //=====================================================================================
  #ifdef SAVE_TRACK_HITS
    printf("Writing TTree Hit Info File... \n");
    fHits->cd();
    EVENT_VECT_GEM->Write();
    EVENT_VECT_MMG1->Write();
    fHits->Close();
    printf("TTree File Written & Closed OK \n");
  #endif

  //=====================================================================================
  //===                 P L O T  H I S T O G R A M S                                  ===
  //=====================================================================================
  const char *OutputDir="RootOutput/cern24";
  #ifdef SAVE_PDF
    int NN_MODE = 8;
    char ctit[120];
    sprintf(G_DIR,"%s/Run_%06d",OutputDir,RunNum);
    sprintf(ctit,"File=%s",G_DIR);
    bool COMPACT=false;
    TCanvas *cc;
    int nxd=3;
    int nyd=5;
    char pdfname[120];  sprintf(pdfname,"%s_evdisp.pdf",G_DIR);  //c0->Print(pdfname);
    
    //---------------------  page 1 --------------------
    htitle(" Count ");   // if (!COMPACT) cc=NextPlot(0,0);
    cc=NextPlot(nxd,nyd);  gPad->SetLogy(); hcount->Draw();
    
    //---------------------  page 2a --------------------
    htitle(" SRS GEM-TRK  ");   if (!COMPACT) cc=NextPlot(0,0);
    cc=NextPlot(nxd,nyd); hgemtrkr_1_peak_xy->Draw("colz");
    cc=NextPlot(nxd,nyd); hgemtrkr_1_atlas_xy->Draw("colz");
    cc=NextPlot(nxd,nyd); hgemtrkr_1_peak_x->Draw();
    cc=NextPlot(nxd,nyd); hgemtrkr_1_peak_y->Draw();
    cc=NextPlot(nxd,nyd); hgemtrkr_1_peak_x_height->Draw();
    cc=NextPlot(nxd,nyd); hgemtrkr_1_peak_y_height->Draw();
    
    htitle(" SRS GEM-TRK  ");   if (!COMPACT) cc=NextPlot(0,0);
    cc=NextPlot(nxd,nyd); hgemtrkr_2_peak_xy->Draw("colz");
    cc=NextPlot(nxd,nyd); hgemtrkr_2_atlas_xy->Draw("colz");
    cc=NextPlot(nxd,nyd); hgemtrkr_2_peak_x->Draw();
    cc=NextPlot(nxd,nyd); hgemtrkr_2_peak_y->Draw();
    cc=NextPlot(nxd,nyd); hgemtrkr_2_peak_x_height->Draw();
    cc=NextPlot(nxd,nyd); hgemtrkr_2_peak_y_height->Draw();
    
    htitle(" SRS GEM-TRK  ");   if (!COMPACT) cc=NextPlot(0,0);
    cc=NextPlot(nxd,nyd); hgemtrkr_3_peak_xy->Draw("colz");
    cc=NextPlot(nxd,nyd); hgemtrkr_3_atlas_xy->Draw("colz");
    cc=NextPlot(nxd,nyd); hgemtrkr_3_peak_x->Draw();
    cc=NextPlot(nxd,nyd); hgemtrkr_3_peak_y->Draw();
    cc=NextPlot(nxd,nyd); hgemtrkr_3_peak_x_height->Draw();
    cc=NextPlot(nxd,nyd); hgemtrkr_3_peak_y_height->Draw();
    
    htitle(" SRS GEM-TRK  ");   if (!COMPACT) cc=NextPlot(0,0);
    nxd=2; nyd=4;
    cc=NextPlot(nxd,nyd); mmg1_peak_y->Draw();
    cc=NextPlot(nxd,nyd); hmmg1_peak_y_height->Draw();
    cc=NextPlot(nxd,nyd); gem_peak_y->Draw();
    cc=NextPlot(nxd,nyd); hgem_peak_y_height->Draw();
    //TBox fbox(xbc1,ybc1,xbc2,ybc2);  //---- draw box cut ---
    //fbox.Draw("same");
    //fbox.SetLineColor(kRed);
    //fbox.SetFillStyle(0);
    //fbox.SetLineWidth(1);
    
    htitle(" SRS GEM-TRK  ");   if (!COMPACT) cc=NextPlot(0,0);
    nxd=2; nyd=4;
    cc=NextPlot(nxd,nyd); gem_mmg1_doubleY->Draw("colz");
    cc=NextPlot(nxd,nyd); gem_mmg1_doubleX->Draw("colz");
    cc=NextPlot(nxd,nyd); hgemtrkr_1_gem->Draw("colz");
    cc=NextPlot(nxd,nyd); hgemtrkr_1_mmg1->Draw("colz");
    
   //---------------------  page 3 --------------------
    htitle("  TRD (fa125) Amp Distributions ");    if (!COMPACT) cc=NextPlot(0,0);
    cc=NextPlot(nxd,nyd);   gPad->SetLogy();  f125_el->Draw();
    cc=NextPlot(nxd,nyd);   gPad->SetLogy();  mmg1_f125_el->Draw();
    cc=NextPlot(nxd,nyd);   gPad->SetLogy();  f125_el_max->Draw();
    cc=NextPlot(nxd,nyd);   gPad->SetLogy();  f125_pi_max->Draw();
    cc=NextPlot(nxd,nyd);   gPad->SetLogy();  mmg1_f125_el_max->Draw();
    cc=NextPlot(nxd,nyd);   gPad->SetLogy();  mmg1_f125_pi_max->Draw();
    
   //---------------------  page 3a --------------------
    htitle("  GEM-TRD (fa125) Amp 2D");    if (!COMPACT) cc=NextPlot(0,0);
    nxd=2; nyd=4;
    cc=NextPlot(nxd,nyd);   f125_el_amp2d->Draw("colz");
    cc=NextPlot(nxd,nyd);   f125_pi_amp2d->Draw("colz");
    cc=NextPlot(nxd,nyd);   f125_el_clu2d->Draw("colz");
    cc=NextPlot(nxd,nyd);   f125_pi_clu2d->Draw("colz");
    /*
    f125_el_amp2d->Scale(1./el_count);  f125_el_amp2d->SetMinimum(0.); f125_el_amp2d->SetMaximum(50.); // f125_el_amp2d->SetStats(0);
    f125_pi_amp2d->Scale(1./pi_count);  f125_pi_amp2d->SetMinimum(0.); f125_pi_amp2d->SetMaximum(50.); // f125_pi_amp2d->SetStats(0);
    f125_el_clu2d->Scale(1./el_count);  f125_el_clu2d->SetMinimum(0.); f125_el_clu2d->SetMaximum(50.); // f125_el_clu2d->SetStats(0);
    f125_pi_clu2d->Scale(1./pi_count);  f125_pi_clu2d->SetMinimum(0.); f125_pi_clu2d->SetMaximum(50.); // f125_pi_clu2d->SetStats(0);
  
  */
    
    htitle("  MMG1-TRD (fa125) Amp 2D");    if (!COMPACT) cc=NextPlot(0,0);
    nxd=2; nyd=4;
    cc=NextPlot(nxd,nyd);   mmg1_f125_el_amp2d->Draw("colz");
    cc=NextPlot(nxd,nyd);   mmg1_f125_pi_amp2d->Draw("colz");
    cc=NextPlot(nxd,nyd);   mmg1_f125_el_clu2d->Draw("colz");
    cc=NextPlot(nxd,nyd);   mmg1_f125_pi_clu2d->Draw("colz");
   
   //---------------------  page 3a --------------------
    htitle("  External Tracking");    if (!COMPACT) cc=NextPlot(0,0);
    //nxd=2; nyd=5;
    cc=NextPlot(nxd,nyd);   f125_tracker_hits->Draw();
    cc=NextPlot(nxd,nyd);   mmg1_f125_tracker_hits->Draw();
    cc=NextPlot(nxd,nyd);   f125_tracker_eff->Draw();
    cc=NextPlot(nxd,nyd);   mmg1_f125_tracker_eff->Draw();
    cc=NextPlot(nxd,nyd);   f125_amp2ds->Draw("colz");
    cc=NextPlot(nxd,nyd);   mmg1_f125_amp2ds->Draw("colz");
    
    htitle("  External Tracking");    if (!COMPACT) cc=NextPlot(0,0);
    nxd=2; nyd=4;
    cc=NextPlot(nxd,nyd);   if (f125_tracker_hits->GetEntries()!=0) { f125_tracker_hits->Sumw2(); f125_tracker_eff->Sumw2();  f125_tracker_eff->Divide(f125_tracker_hits);  f125_tracker_eff->Draw(); }
    cc=NextPlot(nxd,nyd);   if (mmg1_f125_tracker_hits->GetEntries()!=0) { mmg1_f125_tracker_hits->Sumw2(); mmg1_f125_tracker_eff->Sumw2(); mmg1_f125_tracker_eff->Divide(mmg1_f125_tracker_hits);  mmg1_f125_tracker_eff->Draw(); }
    cc=NextPlot(nxd,nyd);   gem_residuals->Draw();
    cc=NextPlot(nxd,nyd);   mmg1_residuals->Draw();
    
   //------------- MAX COMPARISONS ---------------
    htitle("  TRD (fa125) Max Amp Comparisons ");    if (!COMPACT) cc=NextPlot(0,0);
    f125_el_max->Scale(1./(f125_el_max->GetEntries()));  f125_el_max->SetLineColor(2);
    f125_pi_max->Scale(1./(f125_pi_max->GetEntries()));  f125_pi_max->GetYaxis()->SetTitle("Counts / numEntries"); f125_pi_max->SetStats(0); f125_el_max->SetStats(0);
    f125_pi_max->SetTitle("GEM-TRD Max Amplitude (All Time)"); TLegend *l2 = new TLegend(0.75,0.75,0.9,0.9);  l2->AddEntry(f125_el_max, "e", "l");  l2->AddEntry(f125_pi_max, "pi", "l");
    cc=NextPlot(nxd,nyd);   f125_pi_max->Draw();  f125_el_max->Draw("same");  l2->Draw();
    
    mmg1_f125_el_max->Scale(1./(mmg1_f125_el_max->GetEntries()));  mmg1_f125_el_max->SetLineColor(2);
    mmg1_f125_pi_max->Scale(1./(mmg1_f125_pi_max->GetEntries()));  mmg1_f125_pi_max->GetYaxis()->SetTitle("Counts / numEntries");  mmg1_f125_pi_max->SetStats(0);  mmg1_f125_el_max->SetStats(0);
    mmg1_f125_pi_max->SetTitle("MMG1-TRD Max Amplitude (All Time)"); TLegend *l3 = new TLegend(0.75,0.75,0.9,0.9);  l3->AddEntry(mmg1_f125_el_max, "e", "l");  l3->AddEntry(mmg1_f125_pi_max, "pi", "l");
    cc=NextPlot(nxd,nyd);   mmg1_f125_pi_max->Draw();  mmg1_f125_el_max->Draw("same");  l3->Draw();
    
    f125_el_max_late->Scale(1./(f125_el_max_late->GetEntries()));  f125_el_max_late->SetLineColor(2);
    f125_pi_max_late->Scale(1./(f125_pi_max_late->GetEntries()));  f125_pi_max_late->GetYaxis()->SetTitle("Counts / numEntries"); f125_pi_max_late->SetStats(0); f125_el_max_late->SetStats(0);
    f125_pi_max_late->SetTitle("GEM-TRD Max Amplitude (Late Time)"); TLegend *l4 = new TLegend(0.75,0.75,0.9,0.9);  l4->AddEntry(f125_el_max_late, "e", "l");  l4->AddEntry(f125_pi_max_late, "pi", "l");
    cc=NextPlot(nxd,nyd);   f125_pi_max_late->Draw();  f125_el_max_late->Draw("same");  l4->Draw();
    
    mmg1_f125_el_max_late->Scale(1./(mmg1_f125_el_max_late->GetEntries()));  mmg1_f125_el_max_late->SetLineColor(2);
    mmg1_f125_pi_max_late->Scale(1./(mmg1_f125_pi_max_late->GetEntries()));  mmg1_f125_pi_max_late->GetYaxis()->SetTitle("Counts / numEntries");  mmg1_f125_pi_max_late->SetStats(0);  mmg1_f125_el_max_late->SetStats(0);
    mmg1_f125_pi_max_late->SetTitle("MMG1-TRD Max Amplitude (Late Time)"); TLegend *l6 = new TLegend(0.75,0.75,0.9,0.9);  l6->AddEntry(mmg1_f125_el_max_late, "e", "l");  l6->AddEntry(mmg1_f125_pi_max_late, "pi", "l");
    cc=NextPlot(nxd,nyd);   mmg1_f125_pi_max_late->Draw();  mmg1_f125_el_max_late->Draw("same");  l6->Draw();
    
    f125_el_max_early->Scale(1./(f125_el_max_early->GetEntries()));  f125_el_max_early->SetLineColor(2);
    f125_pi_max_early->Scale(1./(f125_pi_max_early->GetEntries()));  f125_pi_max_early->GetYaxis()->SetTitle("Counts / numEntries"); f125_pi_max_early->SetStats(0); f125_el_max_early->SetStats(0);
    f125_pi_max_early->SetTitle("GEM-TRD Max Amplitude (Early Time)"); TLegend *l5 = new TLegend(0.75,0.75,0.9,0.9);  l5->AddEntry(f125_el_max_early, "e", "l");  l5->AddEntry(f125_pi_max_early, "pi", "l");
    cc=NextPlot(nxd,nyd);   f125_pi_max_early->Draw();  f125_el_max_early->Draw("same");  l5->Draw();
    
    mmg1_f125_el_max_early->Scale(1./(mmg1_f125_el_max_early->GetEntries()));  mmg1_f125_el_max_early->SetLineColor(2);
    mmg1_f125_pi_max_early->Scale(1./(mmg1_f125_pi_max_early->GetEntries()));  mmg1_f125_pi_max_early->GetYaxis()->SetTitle("Counts / numEntries");  mmg1_f125_pi_max_early->SetStats(0);  mmg1_f125_el_max_early->SetStats(0);
    mmg1_f125_pi_max_early->SetTitle("MMG1-TRD Max Amplitude (Early Time)"); TLegend *l7 = new TLegend(0.75,0.75,0.9,0.9);  l7->AddEntry(mmg1_f125_el_max_early, "e", "l");  l7->AddEntry(mmg1_f125_pi_max_early, "pi", "l");
    cc=NextPlot(nxd,nyd);   mmg1_f125_pi_max_early->Draw();  mmg1_f125_el_max_early->Draw("same");  l7->Draw();
    
    //--- close PDF file ----
    cc=NextPlot(-1,-1);
  #endif
  cout<<"========== END OF RUN "<<RunNum<<" =========="<<endl;
}
