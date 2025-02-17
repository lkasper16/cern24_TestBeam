#define trdclass_cern24_cxx
#include "trdclass_cern24.h"
#include "PlotLib.C"
#include "GNN/gnn_model.h"
#include "GNN/gnn_model.cpp"
#include "GNN/toGraph.cpp"

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
#define MAX_NODES 100
#define USE_MAXPOS 1
//
#define SAVE_TRACK_HITS
#define SAVE_PDF
//#define WRITE_CSV
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
    char c0Title[256];
    sprintf(c0Title,"Event_Display_Run=%d",RunNum);
    TCanvas *c0 = new TCanvas("DISP",c0Title,1100,200,1500,1300);
    c0->Divide(4,3); c0->cd(1);
    //-----------------  canvas 2 FPGA Display ----------
    char c2Title[256];
    sprintf(c2Title,"FPGA_Event_Display_Run=%d",RunNum);
    TCanvas *c2 = new TCanvas("FPGA",c2Title,1000,100,1500,1300);
    c2->Divide(5,2); c2->cd(1);
    char c3Title[256];
    sprintf(c3Title,"DQM_Event_Display_Run=%d",RunNum);
    TCanvas *c3 = new TCanvas("PID",c3Title,1000,100,1100,900);
    c3->Divide(5,1); //c3->cd(1);
  #endif
  
  // Track fit in time
  TF1 fx1("fx1","pol1",40,150);
  TF1 fx2("fx2","pol1",40,150);
  f125_fit = new TH2F("f125_fit","GEM-TRD Track Fit; Time Response (8ns) ; X Channel",200,0.5,200.5,240,-0.5,239.5);
  mmg1_f125_fit = new TH2F("mmg1_f125_fit","MMG1-TRD Track Fit; Time Response (8ns) ; X Channel",200,0.5,200.5,240,-0.5,239.5);
  //-- TRD - GEMTRKR alignment --------
  double xx1=-37., yy1=-55.,  xx2=53., yy2=44.;
  double aa=(yy2-yy1)/(xx2-xx1);
  double bb=yy1-aa*xx1;
  TF1 ftrk("ftrk","[0]*x+[1]",-55.,55.);
  ftrk.SetParameter(0,aa);
  ftrk.SetParameter(1,bb);
  TF1 ftrkr("ftrkr","(x-[1])/[0]",0.,255.);
  ftrkr.SetParameter(0,aa);
  ftrkr.SetParameter(1,bb);
  //-----------------------
  float z1 = 0;
  float z2 = 1071;
  float z3 = 1153.4;
  float zgem = 571;
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
  
  hgem_nhits = new TH1F("hgem_nhits","GEM-TRD Y Hits (SRS)",100,0,100);  HistList->Add(hgem_nhits);
  hmmg1_nhits = new TH1F("hmmg1_nhits","MMG1-TRD Hits (SRS)",100,0,100);  HistList->Add(hmmg1_nhits);
  hgt1_nhits = new TH1F("hgt1_nhits","GEM-TRKR1 Hits",12,0,12);  HistList->Add(hgt1_nhits);
  hgt2_nhits = new TH1F("hgt2_nhits","GEM-TRKR2 Hits",12,0,12);  HistList->Add(hgt2_nhits);
  hgt3_nhits = new TH1F("hgt3_nhits","GEM-TRKR3 Hits",12,0,12);  HistList->Add(hgt3_nhits);
  
  cout<<"**************************RunNum="<<RunNum<<endl;
  int nx0=100;
  int mfac=110;
  if (RunNum>5284.) {nx0=130; mfac=70;} //-- Second Xe bottle
  cout<<"**************************nx0="<<nx0<<endl;
  int ny0=256;
  double Ymin=0.;    double Ymax=ny0*0.4;
  double Xmin=0.;    double Xmax=30.; //double Xmax=26.;
  mhevt  = new TH2F("mhevt","MMG1-TRD Event Display; z pos [mm]; y pos [mm]",nx0+mfac,Xmin,Xmax,ny0,Ymin,Ymax); mhevt->SetStats(0); mhevt->SetMaximum(10.); mhevt->SetMinimum(0.);
  mhevtc = new TH2F("mhevtc","Clustering; FADC bins; MMG1 strips",nx0+mfac,-0.5,(nx0+mfac)-0.5,ny0,-0.5,ny0-0.5);  mhevtc->SetStats(0);   mhevtc->SetMinimum(0.07); mhevtc->SetMaximum(40.);
  mhevtf = new TH2F("mhevtf","MMG1: Clusters for FPGA; z pos [mm]; y pos [mm]",nx0+mfac,Xmin,Xmax,ny0,Ymin,Ymax);  mhevtf->SetStats(0); mhevtf->SetMaximum(10.);
  
  hevt  = new TH2F("hevt","GEM-TRD Event display; z pos [mm]; y pos [mm]",nx0,Xmin,Xmax,ny0,Ymin,Ymax); hevt->SetStats(0); hevt->SetMaximum(10.); hevt->SetMinimum(0.);
  hevtc = new TH2F("hevtc","Clustering; FADC bins; GEM strips",nx0,-0.5,nx0-0.5,ny0,-0.5,ny0-0.5); hevtc->SetStats(0);   hevtc->SetMinimum(0.07); hevtc->SetMaximum(40.);
  hevtf = new TH2F("hevtf","GEM: Clusters for FPGA; z pos [mm]; y pos [mm]",nx0,Xmin,Xmax,ny0,Ymin,Ymax);  hevtf->SetStats(0); hevtf->SetMaximum(10.);
  #if (USE_PULSE>0)
    hevtk  = new TH2F("hevtk","Event display; z pos [mm]; y pos [mm]",nx0,Xmin,Xmax,ny0,Ymin,Ymax); /*hevtk->SetStats(0);*/ hevtk->SetMaximum(10.);
    hevtck = new TH2F("hevtck","Clustering; FADC bins; GEM strips",nx0,-0.5,nx0-0.5,ny0,-0.5,ny0-0.5);
  #endif
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
  gem_mmg1_doubleX = new TH2F("gem_mmg1_doubleX","TRD X Correlation (fADC) ; GEM-TRD X Chan ; MMG1-TRD X Chan ",240,-0.5,239.5,240,-0.5,239.5);    gem_mmg1_doubleX->SetStats(0); HistList->Add(gem_mmg1_doubleX);
  gem_mmg1_doubleY = new TH2F("gem_mmg1_doubleY","TRD Y Correlation (SRS) ; GEM-TRD Y Chan ; MMG1-TRD Y Chan ",256,-0.5,255.5,256,-0.5,255.5);    gem_mmg1_doubleY->SetStats(0);  HistList->Add(gem_mmg1_doubleY);
  
  // --- SRS ---
  hgemtrkr_1_max_xch = new TH1F("hgemtrkr_1_max_xch"," GEM-TRKR1 Max X Position ; X Chan [mm]",128,-0.2,102.2);  HistList->Add(hgemtrkr_1_max_xch);
  hgemtrkr_1_max_xamp = new TH1F("hgemtrkr_1_max_xamp"," GEM-TRKR1 Max X Amp ; ADC Amp ",100,0.,4096);  HistList->Add(hgemtrkr_1_max_xamp);
  hgemtrkr_2_max_xch = new TH1F("hgemtrkr_2_max_xch"," GEM-TRKR2 Max X Position ; X Chan [mm]",128,-0.2,102.2);  HistList->Add(hgemtrkr_2_max_xch);
  hgemtrkr_2_max_xamp = new TH1F("hgemtrkr_2_max_xamp"," GEM-TRKR2 Max X Amp ; ADC Amp ",100,0.,4096);  HistList->Add(hgemtrkr_2_max_xamp);
  hgemtrkr_3_max_xch = new TH1F("hgemtrkr_3_max_xch"," GEM-TRKR3 Max X Position ; X Chan [mm]",128,-0.2,102.2);  HistList->Add(hgemtrkr_3_max_xch);
  hgemtrkr_3_max_xamp = new TH1F("hgemtrkr_3_max_xamp"," GEM-TRKR3 Max X Amp ; ADC Amp ",100,0.,4096);  HistList->Add(hgemtrkr_3_max_xamp);
  //--GEMTracker 1
  hgemtrkr_1_peak_xy = new TH2F("hgemtrkr_1_peak_xy","GEM-TRKR1 Peak X-Y Correlation; Peak X [mm]; Peak Y [mm] ",128,-0.2,102.2,128,-0.2,102.2);    hgemtrkr_1_peak_xy->SetStats(0); HistList->Add(hgemtrkr_1_peak_xy);
  hgemtrkr_1_max_xy = new TH2F("hgemtrkr_1_max_xy","GEM-TRKR1 X-Y Correlation for Max Hits; Peak X [mm]; Peak Y [mm] ",128,-0.2,102.2,128,-0.2,102.2);    hgemtrkr_1_max_xy->SetStats(0); HistList->Add(hgemtrkr_1_max_xy);
  hgemtrkr_1_peak_x = new TH1F("hgemtrkr_1_peak_x"," GEM-TRKR1 Peak X Position; X [mm] ",128,-0.2,102.2);  HistList->Add(hgemtrkr_1_peak_x);
  hgemtrkr_1_peak_y = new TH1F("hgemtrkr_1_peak_y"," GEM-TRKR1 Peak Y Position; Y [mm] ",128,-0.2,102.2);  HistList->Add(hgemtrkr_1_peak_y);
  hgemtrkr_1_peak_x_height = new TH1F("hgemtrkr_1_peak_x_height"," GEM-TRKR1 Peak Amplitudes in X; ADC Value ",100,0.,4096.);  HistList->Add(hgemtrkr_1_peak_x_height);
  hgemtrkr_1_peak_y_height = new TH1F("hgemtrkr_1_peak_y_height"," GEM-TRKR1 Peak Amplitudes in Y ; ADC Value ",100,0.,4096.);  HistList->Add(hgemtrkr_1_peak_y_height);
  //--GEMTracker 2
  hgemtrkr_2_peak_xy = new TH2F("hgemtrkr_2_peak_xy","GEM-TRKR2 Peak X-Y Correlation; Peak X [mm]; Peak Y [mm] ",128,-0.2,102.2,128,-0.2,102.2); hgemtrkr_2_peak_xy->SetStats(0); HistList->Add(hgemtrkr_2_peak_xy);
  hgemtrkr_2_max_xy = new TH2F("hgemtrkr_2_max_xy","GEM-TRKR2 X-Y Correlation for Max Hits; Peak X [mm]; Peak Y [mm] ",128,-0.2,102.2,128,-0.2,102.2); hgemtrkr_2_max_xy->SetStats(0);   HistList->Add(hgemtrkr_2_max_xy);
  hgemtrkr_2_peak_x = new TH1F("hgemtrkr_2_peak_x"," GEM-TRKR2 Peak X Position; X [mm] ",128,-0.2,102.2);  HistList->Add(hgemtrkr_2_peak_x);
  hgemtrkr_2_peak_y = new TH1F("hgemtrkr_2_peak_y"," GEM-TRKR2 Peak Y Position ; Y [mm] ",128,-0.2,102.2);  HistList->Add(hgemtrkr_2_peak_y);
  hgemtrkr_2_peak_x_height = new TH1F("hgemtrkr_2_peak_x_height"," GEM-TRKR2 Peak Amplitudes in X; ADC Value ",100,0.,4096.);  HistList->Add(hgemtrkr_2_peak_x_height);
  hgemtrkr_2_peak_y_height = new TH1F("hgemtrkr_2_peak_y_height"," GEM-TRKR2 Peak Amplitudes in Y ; ADC Value ",100,0.,4096.);  HistList->Add(hgemtrkr_2_peak_y_height);
  //--GEMTracker 3
  hgemtrkr_3_peak_xy = new TH2F("hgemtrkr_3_peak_xy","GEM-TRKR3 Peak X-Y Correlation; Peak X [mm]; Peak Y [mm] ",128,-0.2,102.2,128,-0.2,102.2); hgemtrkr_3_peak_xy->SetStats(0); HistList->Add(hgemtrkr_3_peak_xy);
  hgemtrkr_3_max_xy = new TH2F("hgemtrkr_3_max_xy","GEM-TRKR3 X-Y Correlation for Max Hits; Peak X [mm]; Peak Y [mm] ",128,-0.2,102.2,128,-0.2,102.2); hgemtrkr_3_max_xy->SetStats(0);   HistList->Add(hgemtrkr_3_max_xy);
  hgemtrkr_3_peak_x = new TH1F("hgemtrkr_3_peak_x"," GEM-TRKR3 Peak X Position ; X [mm] ",128,-0.2,102.2);  HistList->Add(hgemtrkr_3_peak_x);
  hgemtrkr_3_peak_y = new TH1F("hgemtrkr_3_peak_y"," GEM-TRKR3 Peak Y Position ; Y [mm] ",128,-0.2,102.2);  HistList->Add(hgemtrkr_3_peak_y);
  hgemtrkr_3_peak_x_height = new TH1F("hgemtrkr_3_peak_x_height"," GEM-TRKR3 Peak Amplitudes in X; ADC Value ",100,0.,4096.);  HistList->Add(hgemtrkr_3_peak_x_height);
  hgemtrkr_3_peak_y_height = new TH1F("hgemtrkr_3_peak_y_height"," GEM-TRKR3 Peak Amplitudes in Y ; ADC Value ",100,0.,4096.);  HistList->Add(hgemtrkr_3_peak_y_height);
  
  mmg1_peak_y = new TH1F("mmg1_peak_y"," MMG1-TRD Peak Y Position (SRS) ; Y [mm] ",128,-0.2,102.2);  HistList->Add(mmg1_peak_y);
  hmmg1_peak_y_height = new TH1F("hmmg1_peak_y_height"," MMG1-TRD Peak Amplitudes in Y ; ADC Value ",100,0.,4096.);  HistList->Add(hmmg1_peak_y_height);
  hmmg1_peak_y_height_el = new TH1F("hmmg1_peak_y_height_el"," MMG1-TRD Peak Amplitudes in Y (Electrons); ADC Value ",100,0.,4096.);  HistList->Add(hmmg1_peak_y_height_el);
  hmmg1_peak_y_height_pi = new TH1F("hmmg1_peak_y_height_pi"," MMG1-TRD Peak Amplitudes in Y (Pions); ADC Value ",100,0.,4096.);  HistList->Add(hmmg1_peak_y_height_pi);
  gem_peak_y = new TH1F("gem_peak_y"," GEM-TRD Peak Y Position (SRS) ; Y [mm] ",128,-0.2,102.2);  HistList->Add(gem_peak_y);
  hgem_peak_y_height = new TH1F("hgem_peak_y_height"," GEM-TRD Peak Amplitudes in Y ; ADC Value ",100,0.,4096.);  HistList->Add(hgem_peak_y_height);
  hgem_peak_y_height_el = new TH1F("hgem_peak_y_height_el"," GEM-TRD Peak Amplitudes in Y (Electrons); ADC Value ",100,0.,4096.);  HistList->Add(hgem_peak_y_height_el);
  hgem_peak_y_height_pi = new TH1F("hgem_peak_y_height_pi"," GEM-TRD Peak Amplitudes in Y (Pions); ADC Value ",100,0.,4096.);  HistList->Add(hgem_peak_y_height_pi);
  
  hgemtrkr_1_gem = new TH2F("hgemtrkr_1_gem","GEM-TRKR1 & GEM-TRD Y Correlation; GEM-TRD Y [mm]; GEM-TRKR1 Y [mm] ",128,-0.2,102.2,128,-0.2,102.2); hgemtrkr_1_gem->SetStats(0); HistList->Add(hgemtrkr_1_gem);
  hgemtrkr_1_mmg1 = new TH2F("hgemtrkr_1_mmg1","GEM-TRKR1 & MMG1 Y Correlation; MMG1-TRD Y [mm]; GEM-TRKR1 Y [mm] ",128,-0.2,102.2,128,-0.2,102.2); hgemtrkr_1_mmg1->SetStats(0); HistList->Add(hgemtrkr_1_mmg1);
  
  //--External Tracking
  TH1F *f125_el_tracker_hits = new TH1F("f125_el_tracker_hits","GEM-TRD Track Extr. Hits for Electrons; X Chan [mm]",128,-0.2,102.2);   HistList->Add(f125_el_tracker_hits);
  TH1F *f125_pi_tracker_hits = new TH1F("f125_pi_tracker_hits","GEM-TRD Track Extr. Hits for Pions; X Chan [mm]",128,-0.2,102.2);   HistList->Add(f125_pi_tracker_hits);
  TH1F *f125_el_tracker_eff = new TH1F("f125_el_tracker_eff","GEM-TRD Track Eff. for Electrons; X Chan [mm]",128,-0.2,102.2);   HistList->Add(f125_el_tracker_eff);
  TH1F *f125_pi_tracker_eff = new TH1F("f125_pi_tracker_eff","GEM-TRD Track Eff. for Pions; X Chan [mm]",128,-0.2,102.2);   HistList->Add(f125_pi_tracker_eff);
  TH1F *gem_residuals = new TH1F("gem_residuals","GEM-TRD Residual Hits; X Chan [mm] (actual - expected)",250.,-25,25);   HistList->Add(gem_residuals);
  TH1F *gem_residualscorr = new TH1F("gem_residualscorr","GEM-TRD Residual Hits WITH CORR.; X Chan [mm] (actual - expected)",250.,-25,25);   HistList->Add(gem_residualscorr);
  TH2F *gem_residual_ch = new TH2F("gem_residual_ch","GEM-TRD Residual Hits vs Chan; X Chan [mm] (actual); X Chan [mm] (actual - expected)",65, 34.5, 62.5,35,-7.5,6.5); gem_residual_ch->SetStats(0); HistList->Add(gem_residual_ch);
  TH1F *mmg1_f125_el_tracker_hits = new TH1F("mmg1_f125_el_tracker_hits","MMG1-TRD Track Extr. Hits for Electrons; X Chan [mm]",128,-0.2,102.2);   HistList->Add(mmg1_f125_el_tracker_hits);
  TH1F *mmg1_f125_pi_tracker_hits = new TH1F("mmg1_f125_pi_tracker_hits","MMG1-TRD Track Extr. Hits for Pions; X Chan [mm]",128,-0.2,102.2);   HistList->Add(mmg1_f125_pi_tracker_hits);
  TH1F *mmg1_f125_el_tracker_eff = new TH1F("mmg1_f125_el_tracker_eff","MMG1-TRD Track Eff. for Electrons; X Chan [mm]",128,-0.2,102.2);   HistList->Add(mmg1_f125_el_tracker_eff);
  TH1F *mmg1_f125_pi_tracker_eff = new TH1F("mmg1_f125_pi_tracker_eff","MMG1-TRD Track Eff. for Pions; X Chan [mm]",128,-0.2,102.2);   HistList->Add(mmg1_f125_pi_tracker_eff);
  TH1F *mmg1_residuals = new TH1F("mmg1_residuals","MMG1-TRD Residual Hits; X Chan [mm] (actual - expected)",250.,-25,25);   HistList->Add(mmg1_residuals);
  TH1F *mmg1_residualscorr = new TH1F("mmg1_residualscorr","MMG1-TRD Residual Hits WITH CORR.; X Chan [mm] (actual - expected)",250.,-25,25);   HistList->Add(mmg1_residualscorr);
  TH2F *mmg1_residual_ch = new TH2F("mmg1_residual_ch","MMG1-TRD Residual Hits vs Chan; X Chan [mm] (actual); X Chan [mm] (actual - expected)",65, 36.5, 62.5,35,-7.5,6.5); mmg1_residual_ch->SetStats(0); HistList->Add(mmg1_residual_ch);
  
  //---- GEM-TRD --
  float GEM_THRESH=140.; //175
  float MMG1_THRESH=125.; //150
  
  //f125_amp2ds = new TH2F("f125_amp2ds","GEM-TRD ADC Amp in Time (EXT. TRK) ; Time Response (8ns) ; X Channel ",200,0.5,200.5,240,0.5,240.5); f125_amp2ds->SetStats(0); HistList->Add(f125_amp2ds);
  f125_el_amp2d = new TH2F("f125_el_amp2d","GEM-TRD ADC Amp in Time for Electrons ; Time Response (8ns) ; X Channel ",200,0.5,200.5,240,-0.5,239.5); f125_el_amp2d->SetStats(0); HistList->Add(f125_el_amp2d);
  f125_pi_amp2d = new TH2F("f125_pi_amp2d","GEM-TRD ADC Amp in Time for Pions ; Time Response (8ns) ; X Channel ",200,0.5,200.5,240,-0.5,239.5); f125_pi_amp2d->SetStats(0); HistList->Add(f125_pi_amp2d);
  //mmg1_f125_amp2ds = new TH2F("mmg1_f125_amp2ds","MMG1-TRD ADC Amp in Time (EXT. TRK) ; Time Response (8ns) ; X Channel ",200,0.5,200.5,240,0.5,240.5); mmg1_f125_amp2ds->SetStats(0); HistList->Add(mmg1_f125_amp2ds);
  mmg1_f125_el_amp2d = new TH2F("mmg1_f125_el_amp2d","MMG1-TRD ADC Amp in Time for Electrons ; Time Response (8ns) ; X Channel ",200,0.5,200.5,240,-0.5,239.5); mmg1_f125_el_amp2d->SetStats(0); HistList->Add(mmg1_f125_el_amp2d);
  mmg1_f125_pi_amp2d = new TH2F("mmg1_f125_pi_amp2d","MMG1-TRD ADC Amp in Time for Pions ; Time Response (8ns) ; X Channel ",200,0.5,200.5,240,-0.5,239.5); mmg1_f125_pi_amp2d->SetStats(0); HistList->Add(mmg1_f125_pi_amp2d);
  
  f125_el_clu2d = new TH2F("f125_el_clu2d","GEM-TRD Amp for Electrons (Clusters)",200,0.5,200.5,240,-0.5,239.5); f125_el_clu2d->SetStats(0); HistList->Add(f125_el_clu2d);
  f125_pi_clu2d = new TH2F("f125_pi_clu2d","GEM-TRD Amp for Pions (Clusters)",200,0.5,200.5,240,-0.5,239.5); f125_pi_clu2d->SetStats(0); HistList->Add(f125_pi_clu2d);
  mmg1_f125_el_clu2d = new TH2F("mmg1_f125_el_clu2d","MMG1-TRD Amp for Electrons (Clusters)",200,0.5,200.5,240,-0.5,239.5); mmg1_f125_el_clu2d->SetStats(0); HistList->Add(mmg1_f125_el_clu2d);
  mmg1_f125_pi_clu2d = new TH2F("mmg1_f125_pi_clu2d","MMG1-TRD Amp for Pions (Clusters)",200,0.5,200.5,240,-0.5,239.5); mmg1_f125_pi_clu2d->SetStats(0); HistList->Add(mmg1_f125_pi_clu2d);
  
  hgem_xy = new TH2F("hgem_xy","GEM-TRD X-Y Hit Display; X (fADC) [mm]; Y (SRS) [mm] ",128,-0.2,102.2,128,-0.2,102.2); HistList->Add(hgem_xy);
  hmmg1_xy = new TH2F("hmmg1_xy","MMG1-TRD X-Y Hit Display; X (fADC) [mm]; Y (SRS) [mm] ",128,-0.2,102.2,128,-0.2,102.2); HistList->Add(hmmg1_xy);
  gem_mmg1_xcorr = new TH2F("gem_mmg1_xcorr","GEM-TRD X Correlation With MMG1-TRD; GEM-TRD X [mm]; MMG1-TRD X [mm] ",128,-0.4,102.,128,-0.4,102.); HistList->Add(gem_mmg1_xcorr);
  gem_mmg1_ycorr = new TH2F("gem_mmg1_ycorr","GEM-TRD Y Correlation With MMG1-TRD; GEM-TRD Y [mm]; MMG1-TRD Y [mm] ",128,-0.2,102.2,128,-0.2,102.2); HistList->Add(gem_mmg1_ycorr);
  gem_mmg1_max_xcorr = new TH2F("gem_mmg1_max_xcorr","GEM-TRD Max X Correlation With MMG1-TRD; GEM-TRD X [mm]; MMG1-TRD X [mm] ",128,-0.4,102.,128,-0.4,102.); HistList->Add(gem_mmg1_max_xcorr);
  gem_gt1_xcorr = new TH2F("gem_gt1_xcorr","GEM-TRD X Correlation With GEMTRKR-1; GEM-TRD X [mm]; GEMTRKR-1 X [mm] ",128,-0.4,102.,128,-0.2,102.2); HistList->Add(gem_gt1_xcorr);
  gem_gt2_xcorr = new TH2F("gem_gt2_xcorr","GEM-TRD X Correlation With GEMTRKR-2; GEM-TRD X [mm]; GEMTRKR-2 X [mm] ",128,-0.4,102.,128,-0.2,102.2); HistList->Add(gem_gt2_xcorr);
  gem_gt3_xcorr = new TH2F("gem_gt3_xcorr","GEM-TRD X Correlation With GEMTRKR-3; GEM-TRD X [mm]; GEMTRKR-3 X [mm] ",128,-0.4,102.,128,-0.2,102.2); HistList->Add(gem_gt3_xcorr);
  mmg1_gt1_xcorr = new TH2F("mmg1_gt1_xcorr","MMG1-TRD X Correlation With GEMTRKR-1; MMG1-TRD X [mm]; GEMTRKR-1 X [mm] ",128,-0.4,102.,128,-0.2,102.2); HistList->Add(mmg1_gt1_xcorr);
  mmg1_gt2_xcorr = new TH2F("mmg1_gt2_xcorr","MMG1-TRD X Correlation With GEMTRKR-2; MMG1-TRD X [mm]; GEMTRKR-2 X [mm] ",128,-0.4,102.,128,-0.2,102.2); HistList->Add(mmg1_gt2_xcorr);
  mmg1_gt3_xcorr = new TH2F("mmg1_gt3_xcorr","MMG1-TRD X Correlation With GEMTRKR-3; MMG1-TRD X [mm]; GEMTRKR-3 X [mm] ",128,-0.4,102.,128,-0.2,102.2); HistList->Add(mmg1_gt3_xcorr);
  
  hgemClusterDiff_el = new TH1F("hgemClusterDiff_el","GEM Cluster Distance from MMG1 Track (Electrons); Distance [mm]",160,-20.,20.); HistList->Add(hgemClusterDiff_el);
  hgemClusterDiff_pi = new TH1F("hgemClusterDiff_pi","GEM Cluster Distance from MMG1 Track (Pions); Distance [mm]",160,-20.,20.); HistList->Add(hgemClusterDiff_pi);
  hmmg1ClusterDiff_el = new TH1F("hmmg1ClusterDiff_el","MMG1 Cluster Distance from GEM Track (Electrons); Distance [mm]",160,-20.,20.); HistList->Add(hmmg1ClusterDiff_el);
  hmmg1ClusterDiff_pi = new TH1F("hmmg1ClusterDiff_pi","MMG1 Cluster Distance from GEM Track (Pions); Distance [mm]",160,-20.,20.); HistList->Add(hmmg1ClusterDiff_pi);
  hgemPulseDiff_el = new TH1F("hgemPulseDiff_el","GEM Pulse Distance from MMG1 Track (Electrons); Distance [mm]",160,-20.,20.); HistList->Add(hgemPulseDiff_el);
  hgemPulseDiff_pi = new TH1F("hgemPulseDiff_pi","GEM Pulse Distance from MMG1 Track (Pions); Distance [mm]",160,-20.,20.); HistList->Add(hgemPulseDiff_pi);
  hmmg1PulseDiff_el = new TH1F("hmmg1PulseDiff_el","MMG1 Pulse Distance from GEM Track (Electrons); Distance [mm]",160,-20.,20.); HistList->Add(hmmg1PulseDiff_el);
  hmmg1PulseDiff_pi = new TH1F("hmmg1PulseDiff_pi","MMG1 Pulse Distance from GEM Track (Pions); Distance [mm]",160,-20.,20.); HistList->Add(hmmg1PulseDiff_pi);
  
  hClusterMaxdEdx_el = new TH1F("hClusterMaxdEdx_el","GEM Max Cluster Energy (Electrons); Cluster Energy ; Counts ",100,0.,4960); HistList->Add(hClusterMaxdEdx_el);
  hClusterMaxdEdx_pi = new TH1F("hClusterMaxdEdx_pi","GEM Max Cluster Energy (Pions); Cluster Energy ; Counts ",100,0.,4960); HistList->Add(hClusterMaxdEdx_pi);
  hClusterTotaldEdx_el = new TH1F("hClusterTotaldEdx_el","GEM Total Cluster Energy (Electrons); Total Cluster Energy ; Counts ",100,0.,4960); HistList->Add(hClusterTotaldEdx_el);
  hClusterTotaldEdx_pi = new TH1F("hClusterTotaldEdx_pi","GEM Total Cluster Energy (Pions); Total Cluster Energy ; Counts ",100,0.,4960); HistList->Add(hClusterTotaldEdx_pi);
  hmmg1ClusterMaxdEdx_el = new TH1F("hmmg1ClusterMaxdEdx_el","MMG1 Max Cluster Energy (Electrons); Cluster Energy ; Counts ",100,0.,4960); HistList->Add(hmmg1ClusterMaxdEdx_el);
  hmmg1ClusterMaxdEdx_pi = new TH1F("hmmg1ClusterMaxdEdx_pi","MMG1 Max Cluster Energy (Pions); Cluster Energy ; Counts ",100,0.,4960); HistList->Add(hmmg1ClusterMaxdEdx_pi);
  hmmg1ClusterTotaldEdx_el = new TH1F("hmmg1ClusterTotaldEdx_el","MMG1 Total Cluster Energy (Electrons); Total Cluster Energy ; Counts ",100,0.,4960); HistList->Add(hmmg1ClusterTotaldEdx_el);
  hmmg1ClusterTotaldEdx_pi = new TH1F("hmmg1ClusterTotaldEdx_pi","MMG1 Total Cluster Energy (Pions); Total Cluster Energy ; Counts ",100,0.,4960); HistList->Add(hmmg1ClusterTotaldEdx_pi);
  
  hchan_g_el = new TH1F("hchan_g_el","GEM Max X Position for Electrons; GEM X Max [mm]",128,-0.2,102.2);  HistList->Add(hchan_g_el);
  hchan_g_pi = new TH1F("hchan_g_pi","GEM Max X Position for Pions; GEM X Max [mm]",128,-0.2,102.2);  HistList->Add(hchan_g_pi);
  hchan_m_el = new TH1F("hchan_m_el","MMG1 Max X Position for Electrons; MMG1 X Max [mm]",128,-0.2,102.2);  HistList->Add(hchan_m_el);
  hchan_m_pi = new TH1F("hchan_m_pi","MMG1 Max X Position for Pions; MMG1 X Max [mm]",128,-0.2,102.2);  HistList->Add(hchan_m_pi);
  
  //=========================================
  
  TFile* fHits;
  #ifdef SAVE_TRACK_HITS
    #if ANALYZE_MERGED
    char hitsFileName[256]; sprintf(hitsFileName, "RootOutput/cern24/merged/trd_singleTrackHits_Run_%06d_%06dEntries.root",RunNum,nEntries);
    #else
    char hitsFileName[256]; sprintf(hitsFileName, "RootOutput/cern24/trd_singleTrackHits_Run_%06d.root", RunNum);
    #endif
    fHits = new TFile(hitsFileName, "RECREATE");
    gem_zHist = new  TH1F("gem_zHist", "gem_zHist", 20, 80., 200.);
    mmg1_zHist = new  TH1F("mmg1_zHist", "mmg1_zHist", 20, 80., 200.);
    //-- GEM-TRD
    EVENT_VECT_GEM = new TTree("gem_hits","GEM TTree with single track hit info");
    EVENT_VECT_GEM->Branch("event_num",&event_num,"event_num/I");
    EVENT_VECT_GEM->Branch("nhit",&gem_nhit,"gem_nhit/I");
    EVENT_VECT_GEM->Branch("nclu",&gem_nclu,"gem_nclu/I");
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
    EVENT_VECT_GEM->Branch("xposc_max",&clu_xpos_max,"clu_xpos_max/f");
    EVENT_VECT_GEM->Branch("zposc_max",&clu_zpos_max,"clu_zpos_max/f");
    EVENT_VECT_GEM->Branch("dedxc_max",&clu_dedx_max,"clu_dedx_max/f");
    EVENT_VECT_GEM->Branch("widthc_max",&clu_width_max,"clu_width_max/f");
    EVENT_VECT_GEM->Branch("dedxc_tot",&clu_dedx_tot,"clu_dedx_tot/f");
    EVENT_VECT_GEM->Branch("xch_max",&gem_xch_max,"gem_xch_max/I");
    //EVENT_VECT_GEM->Branch("xchtrkr",&gemtrkr_xch_max);
    //EVENT_VECT_GEM->Branch("ychtrkr",&gemtrkr_ych_max);
    EVENT_VECT_GEM->Branch("amp_max",&gem_amp_max,"gem_amp_max/f");
    EVENT_VECT_GEM->Branch("time_max",&gem_time_max,"gem_time_max/f");
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
    EVENT_VECT_MMG1->Branch("nclu",&mmg1_nclu,"mmg1_nclu/I");
    EVENT_VECT_MMG1->Branch("ecal_energy",&Ecal_Energy,"Ecal_Energy/float");
    EVENT_VECT_MMG1->Branch("presh_energy",&Presh_Energy,"Presh_Energy/float");
    EVENT_VECT_MMG1->Branch("mult_energy",&Mult_Energy,"Mult_Energy/float");
    EVENT_VECT_MMG1->Branch("xpos",&mmg1_xpos);
    EVENT_VECT_MMG1->Branch("zpos",&mmg1_zpos);
    EVENT_VECT_MMG1->Branch("dedx",&mmg1_dedx);
    EVENT_VECT_MMG1->Branch("parID",&mmg1_parID);
    EVENT_VECT_MMG1->Branch("zHist",&mmg1_zHist_vect);
    EVENT_VECT_MMG1->Branch("xposc",&mmg1_clu_xpos);
    EVENT_VECT_MMG1->Branch("zposc",&mmg1_clu_zpos);
    EVENT_VECT_MMG1->Branch("dedxc",&mmg1_clu_dedx);
    EVENT_VECT_MMG1->Branch("widthc",&mmg1_clu_width);
    EVENT_VECT_MMG1->Branch("xposc_max",&mmg1_clu_xpos_max,"mmg1_clu_xpos_max/f");
    EVENT_VECT_MMG1->Branch("zposc_max",&mmg1_clu_zpos_max,"mmg1_clu_zpos_max/f");
    EVENT_VECT_MMG1->Branch("dedxc_max",&mmg1_clu_dedx_max,"mmg1_clu_dedx_max/f");
    EVENT_VECT_MMG1->Branch("widthc_max",&mmg1_clu_width_max,"mmg1_clu_width_max/f");
    EVENT_VECT_MMG1->Branch("dedxc_tot",&mmg1_clu_dedx_tot,"mmg1_clu_dedx_tot/f");
    EVENT_VECT_MMG1->Branch("xch_max",&mmg1_xch_max,"mmg1_xch_max/I");
    EVENT_VECT_MMG1->Branch("amp_max",&mmg1_amp_max,"mmg1_amp_max/f");
    EVENT_VECT_MMG1->Branch("time_max",&mmg1_time_max,"mmg1_time_max/f");
    EVENT_VECT_MMG1->Branch("chi2",&mmg1_chi2cc);
    EVENT_VECT_MMG1->Branch("Fint",&mmg1_integral);
    EVENT_VECT_MMG1->Branch("a0",&mmg1_a0);
    EVENT_VECT_MMG1->Branch("a1",&mmg1_a1);
    
  #endif
  
  TStopwatch timer;
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes=0, nb=0;
  if (MaxEvt>0) nentries=MaxEvt;  //-- limit number of events for test
  
  //==================================================================================================
  //                      E v e n t    L o o p
  //==================================================================================================
  printf("===============  Begin Event Loop - 1st evt=%lld, Last evt=%lld =============== \n",FirstEvt,MaxEvt);
  timer.Start();
  Long64_t jentry=0;
  int el_count=0, pi_count=0, atlas_trigger_count=0;
  
  for (jentry=FirstEvt; jentry<nentries; jentry++) { //-- Event Loop --
    Count("EVT");
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);
    nbytes += nb;
    if (!(jentry%NPRT))
      printf("------- evt=%llu  f125_raw_count=%llu f125_pulse_count=%llu f250_wraw_count=%llu, srs_peak_count=%llu \n",jentry,f125_wraw_count, f125_pulse_count, f250_wraw_count, gem_peak_count);
    event_num=jentry;
    gem_nhit=0;
    mmg1_nhit=0;
    gem_nclu=0;
    mmg1_nclu=0;
    
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
    clu_xpos_max=0;
    clu_zpos_max=0;
    clu_dedx_max=0;
    clu_width_max=0;
    clu_dedx_tot=0;
    gem_amp_max=0;
    gem_xch_max=0;
    gem_time_max=0;
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
    mmg1_clu_xpos_max=0;
    mmg1_clu_zpos_max=0;
    mmg1_clu_dedx_max=0;
    mmg1_clu_width_max=0;
    mmg1_clu_dedx_tot=0;
    mmg1_amp_max=0;
    mmg1_xch_max=0;
    mmg1_time_max=0;
    mmg1_chi2cc.clear();
    mmg1_integral.clear();
    
    //==================================================================================================
    //                    Process Fa250 Pulse data (PID)
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
    
    double CalSum=0;
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
      if (fadc_chan==ch_cal) hCal_pulse->Reset();
      if (fadc_chan==ch_presh) hPresh_pulse->Reset();
      if (fadc_chan==ch_cher) hCher_pulse->Reset();
      if (fadc_chan==ch_mult) hMult_pulse->Reset();
      
      //--ped calculation for f250 raw data
      int nped = 0, ped=100;
      double ped_sum = 0.;
      if (fadc_window > 15) {
        for (int si=20; si<35; si++) {
          int ped_samp = f250_wraw_samples->at(f250_wraw_samples_index->at(i)+si);
          ped_sum += ped_samp;
          nped++;
        }
      }
      if (nped>0.) ped = ped_sum / nped;
      if (0. > ped || ped > 200 ) ped = 100;
      
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
    //  GEMTracker (SRS) Correlations with TRD Prototypes
    //===========================================================
    
    ULong64_t gt_1_idx_x = 0, gt_1_idx_y=0, gt_2_idx_x = 0, gt_2_idx_y = 0, gt_3_idx_x = 0, gt_3_idx_y = 0, mmg1_idx_y=0, gem_idx_y=0;
    int gt1_nhit=0, gt2_nhit=0, gt3_nhit=0;
    
    double gemtrkr1_xamp_max=-1., gemtrkr1_xch_max=-1000.;
    double gemtrkr1_yamp_max=-1., gemtrkr1_ych_max=-1000.;
    double gemtrkr2_xamp_max=-1., gemtrkr2_xch_max=-1000.;
    double gemtrkr2_yamp_max=-1., gemtrkr2_ych_max=-1000.;
    double gemtrkr3_xamp_max=-1., gemtrkr3_xch_max=-1000.;
    double gemtrkr3_yamp_max=-1., gemtrkr3_ych_max=-1000.;
    double mmg1_el_yamp_max=-1., mmg1_pi_yamp_max=-1., mmg1_el_ych_max=-1000., mmg1_pi_ych_max=-1000.;
    double gem_pi_yamp_max=-1., gem_el_yamp_max=-1., gem_el_ych_max=-1000., gem_pi_ych_max=-1000.;
    
    double gemtrkr_1_peak_pos_y[gem_peak_count];
    double gemtrkr_1_peak_pos_x[gem_peak_count];
    double gemtrkr_1_peak_x_height[gem_peak_count];
    double gemtrkr_1_peak_y_height[gem_peak_count];
    double gemtrkr_2_peak_pos_y[gem_peak_count];
    double gemtrkr_2_peak_pos_x[gem_peak_count];
    double gemtrkr_2_peak_x_height[gem_peak_count];
    double gemtrkr_2_peak_y_height[gem_peak_count];
    double gemtrkr_3_peak_pos_y[gem_peak_count];
    double gemtrkr_3_peak_pos_x[gem_peak_count];
    double gemtrkr_3_peak_x_height[gem_peak_count];
    double gemtrkr_3_peak_y_height[gem_peak_count];
    double mmg1_peak_pos_y[gem_peak_count];
    double mmg1_peak_y_height[gem_peak_count];
    double gem_peak_pos_y[gem_peak_count];
    double gem_peak_y_height[gem_peak_count];
    
    for (ULong64_t i=0; i<gem_peak_count; i++) {
      gemtrkr_1_peak_pos_y[gem_peak_count] = -1000;
      gemtrkr_1_peak_pos_x[gem_peak_count] = -1000;
      gemtrkr_1_peak_x_height[gem_peak_count] = -1000;
      gemtrkr_1_peak_y_height[gem_peak_count] = -1000;
      gemtrkr_2_peak_pos_y[gem_peak_count] = -1000;
      gemtrkr_2_peak_pos_x[gem_peak_count] = -1000;
      gemtrkr_2_peak_x_height[gem_peak_count] = -1000;
      gemtrkr_2_peak_y_height[gem_peak_count] = -1000;
      gemtrkr_3_peak_pos_y[gem_peak_count] = -1000;
      gemtrkr_3_peak_pos_x[gem_peak_count] = -1000;
      gemtrkr_3_peak_x_height[gem_peak_count] = -1000;
      gemtrkr_3_peak_y_height[gem_peak_count] = -1000;
      mmg1_peak_pos_y[gem_peak_count] = -1000;
      mmg1_peak_y_height[gem_peak_count] = -1000;
      gem_peak_pos_y[gem_peak_count] = -1000;
      gem_peak_y_height[gem_peak_count] = -1000;
    }
   
    for (ULong64_t i=0; i<gem_peak_count; i++) { //-- SRS Peaks Loop
      
      if (gem_peak_plane_name->at(i) == "GEMTR1X") {
        gemtrkr_1_peak_pos_x[gt_1_idx_x] = gem_peak_real_pos->at(i);
        if (gemtrkr_1_peak_pos_x[gt_1_idx_x]<0) gemtrkr_1_peak_pos_x[gt_1_idx_x]+=50.; else if (gemtrkr_1_peak_pos_x[gt_1_idx_x]>0) gemtrkr_1_peak_pos_x[gt_1_idx_x]-=50.;  gemtrkr_1_peak_pos_x[gt_1_idx_x]*=-1.;  gemtrkr_1_peak_pos_x[gt_1_idx_x]+=50.;
        gemtrkr_1_peak_x_height[gt_1_idx_x] = gem_peak_height->at(i);
        if (gemtrkr_1_peak_x_height[gt_1_idx_x]>1000.) gt_1_idx_x++; Count("gt1_x");
        //gt_1_idx_x++; Count("gt1_x");
      } if (gem_peak_plane_name->at(i) == "GEMTR1Y") {
          gemtrkr_1_peak_pos_y[gt_1_idx_y] = gem_peak_real_pos->at(i);
          if (gemtrkr_1_peak_pos_y[gt_1_idx_y]<0) gemtrkr_1_peak_pos_y[gt_1_idx_y]+=50.; else if (gemtrkr_1_peak_pos_y[gt_1_idx_y]>0) gemtrkr_1_peak_pos_y[gt_1_idx_y]-=50.;  gemtrkr_1_peak_pos_y[gt_1_idx_y]*=-1.;  gemtrkr_1_peak_pos_y[gt_1_idx_y]+=50.;
          gemtrkr_1_peak_y_height[gt_1_idx_y] = gem_peak_height->at(i);
          if (gemtrkr_1_peak_y_height[gt_1_idx_y]>1000.) gt_1_idx_y++; Count("gt1_y");
          //gt_1_idx_y++; Count("gt1_y");
      } if (gem_peak_plane_name->at(i) == "GEMTR2X") {
        gemtrkr_2_peak_pos_x[gt_2_idx_x] = gem_peak_real_pos->at(i);
        if (gemtrkr_2_peak_pos_x[gt_2_idx_x]<0) gemtrkr_2_peak_pos_x[gt_2_idx_x]+=50.; else if (gemtrkr_2_peak_pos_x[gt_2_idx_x]>0) gemtrkr_2_peak_pos_x[gt_2_idx_x]-=50.;  gemtrkr_2_peak_pos_x[gt_2_idx_x]*=-1.;  gemtrkr_2_peak_pos_x[gt_2_idx_x]+=50.;
        gemtrkr_2_peak_x_height[gt_2_idx_x] = gem_peak_height->at(i);
        if (gemtrkr_2_peak_x_height[gt_2_idx_x]>1000.) gt_2_idx_x++; Count("gt2_x");
        //gt_2_idx_x++; Count("gt2_x");
      } if (gem_peak_plane_name->at(i) == "GEMTR2Y") {
          gemtrkr_2_peak_pos_y[gt_2_idx_y] = gem_peak_real_pos->at(i);
          if (gemtrkr_2_peak_pos_y[gt_2_idx_y]<0) gemtrkr_2_peak_pos_y[gt_2_idx_y]+=50.; else if (gemtrkr_2_peak_pos_y[gt_2_idx_y]>0) gemtrkr_2_peak_pos_y[gt_2_idx_y]-=50.;  gemtrkr_2_peak_pos_y[gt_2_idx_y]*=-1.;  gemtrkr_2_peak_pos_y[gt_2_idx_y]+=50.;
          gemtrkr_2_peak_y_height[gt_2_idx_y] = gem_peak_height->at(i);
          if (gemtrkr_2_peak_y_height[gt_2_idx_y]>1000.) gt_2_idx_y++; Count("gt2_y");
          //gt_2_idx_y++; Count("gt2_y");
      } if (gem_peak_plane_name->at(i) == "GEMTR3X") {
        gemtrkr_3_peak_pos_x[gt_3_idx_x] = gem_peak_real_pos->at(i);
        if (gemtrkr_3_peak_pos_x[gt_3_idx_x]<0) gemtrkr_3_peak_pos_x[gt_3_idx_x]+=50.; else if (gemtrkr_3_peak_pos_x[gt_3_idx_x]>0) gemtrkr_3_peak_pos_x[gt_3_idx_x]-=50.;  gemtrkr_3_peak_pos_x[gt_3_idx_x]*=-1.;  gemtrkr_3_peak_pos_x[gt_3_idx_x]+=50.;
        gemtrkr_3_peak_x_height[gt_3_idx_x] = gem_peak_height->at(i);
        if (gemtrkr_3_peak_x_height[gt_3_idx_x]>1000.) gt_3_idx_x++; Count("gt3_x");
        //gt_3_idx_x++; Count("gt3_x");
      } if (gem_peak_plane_name->at(i) == "GEMTR3Y") {
          gemtrkr_3_peak_pos_y[gt_3_idx_y] = gem_peak_real_pos->at(i);
          if (gemtrkr_3_peak_pos_y[gt_3_idx_y]<0) gemtrkr_3_peak_pos_y[gt_3_idx_y]+=50.; else if (gemtrkr_3_peak_pos_y[gt_3_idx_y]>0) gemtrkr_3_peak_pos_y[gt_3_idx_y]-=50.;  gemtrkr_3_peak_pos_y[gt_3_idx_y]*=-1.;  gemtrkr_3_peak_pos_y[gt_3_idx_y]+=50.;
          gemtrkr_3_peak_y_height[gt_3_idx_y] = gem_peak_height->at(i);
          if (gemtrkr_3_peak_y_height[gt_3_idx_y]>1000.) gt_3_idx_y++; Count("gt3_y");
          //gt_3_idx_y++; Count("gt3_y");
      } if (gem_peak_plane_name->at(i) == "MMG1TRDY") {
        mmg1_peak_pos_y[mmg1_idx_y] = gem_peak_real_pos->at(i);
        if (mmg1_peak_pos_y[mmg1_idx_y]<0) mmg1_peak_pos_y[mmg1_idx_y]+=50.; else if (mmg1_peak_pos_y[mmg1_idx_y]>0) mmg1_peak_pos_y[mmg1_idx_y]-=50.;  mmg1_peak_pos_y[mmg1_idx_y]*=-1.;  mmg1_peak_pos_y[mmg1_idx_y]+=50.;
        mmg1_peak_y_height[mmg1_idx_y] = gem_peak_height->at(i);
        if (mmg1_peak_y_height[mmg1_idx_y]>1000.) mmg1_idx_y++; Count("mmg1_y");
        //mmg1_idx_y++; Count("mmg1_y");
      } if (gem_peak_plane_name->at(i) == "VU_GEMTRDY") {
        gem_peak_pos_y[gem_idx_y] = gem_peak_real_pos->at(i);
        if (gem_peak_pos_y[gem_idx_y]<0) gem_peak_pos_y[gem_idx_y]+=50.; else if (gem_peak_pos_y[gem_idx_y]>0) gem_peak_pos_y[gem_idx_y]-=50.;  gem_peak_pos_y[gem_idx_y]*=-1.;  gem_peak_pos_y[gem_idx_y]+=50.;
        gem_peak_y_height[gem_idx_y] = gem_peak_height->at(i);
        if (gem_peak_y_height[gem_idx_y]>1000.) gem_idx_y++; Count("gem_y");
        //gem_idx_y++; Count("gem_y");
      }
    }
    
    for (ULong64_t j=0; j<gt_1_idx_y; j++) {
      hgemtrkr_1_peak_y->Fill(gemtrkr_1_peak_pos_y[j]);
      hgemtrkr_1_peak_y_height->Fill(gemtrkr_1_peak_y_height[j]);
      if (gemtrkr_1_peak_y_height[j]>gemtrkr1_yamp_max) {
        gemtrkr1_yamp_max=gemtrkr_1_peak_y_height[j];
        gemtrkr1_ych_max=gemtrkr_1_peak_pos_y[j];
      }
      for (ULong64_t k=0; k<gt_1_idx_x; k++) {
        hgemtrkr_1_peak_xy->Fill(gemtrkr_1_peak_pos_x[k], gemtrkr_1_peak_pos_y[j]);
        gt1_nhit++;
      }
    }
    for (ULong64_t j=0; j<gt_2_idx_y; j++) {
      hgemtrkr_2_peak_y->Fill(gemtrkr_2_peak_pos_y[j]);
      hgemtrkr_2_peak_y_height->Fill(gemtrkr_2_peak_y_height[j]);
      if (gemtrkr_2_peak_y_height[j]>gemtrkr2_yamp_max) {
        gemtrkr2_yamp_max=gemtrkr_2_peak_y_height[j];
        gemtrkr2_ych_max=gemtrkr_2_peak_pos_y[j];
      }
      for (ULong64_t k=0; k<gt_2_idx_x; k++) {
        hgemtrkr_2_peak_xy->Fill(gemtrkr_2_peak_pos_x[k], gemtrkr_2_peak_pos_y[j]);
        gt2_nhit++;
      }
    }
    for (ULong64_t j=0; j<gt_3_idx_y; j++) {
      hgemtrkr_3_peak_y->Fill(gemtrkr_3_peak_pos_y[j]);
      hgemtrkr_3_peak_y_height->Fill(gemtrkr_3_peak_y_height[j]);
      if (gemtrkr_3_peak_y_height[j]>gemtrkr3_yamp_max) {
        gemtrkr3_yamp_max=gemtrkr_3_peak_y_height[j];
        gemtrkr3_ych_max=gemtrkr_3_peak_pos_y[j];
      }
      for (ULong64_t k=0; k<gt_3_idx_x; k++) {
        hgemtrkr_3_peak_xy->Fill(gemtrkr_3_peak_pos_x[k], gemtrkr_3_peak_pos_y[j]);
        gt3_nhit++;
      }
    }
    for (ULong64_t k=0; k<gt_1_idx_x; k++) {
      hgemtrkr_1_peak_x->Fill(gemtrkr_1_peak_pos_x[k]);
      hgemtrkr_1_peak_x_height->Fill(gemtrkr_1_peak_x_height[k]);
      if (gemtrkr_1_peak_x_height[k]>gemtrkr1_xamp_max) {
        gemtrkr1_xamp_max=gemtrkr_1_peak_x_height[k];
        gemtrkr1_xch_max=gemtrkr_1_peak_pos_x[k];
      }
    }
    for (ULong64_t k=0; k<gt_2_idx_x; k++) {
      hgemtrkr_2_peak_x->Fill(gemtrkr_2_peak_pos_x[k]);
      hgemtrkr_2_peak_x_height->Fill(gemtrkr_2_peak_x_height[k]);
      if (gemtrkr_2_peak_x_height[k]>gemtrkr2_xamp_max) {
        gemtrkr2_xamp_max=gemtrkr_2_peak_x_height[k];
        gemtrkr2_xch_max=gemtrkr_2_peak_pos_x[k];
      }
    }
    for (ULong64_t k=0; k<gt_3_idx_x; k++) {
      hgemtrkr_3_peak_x->Fill(gemtrkr_3_peak_pos_x[k]);
      hgemtrkr_3_peak_x_height->Fill(gemtrkr_3_peak_x_height[k]);
      if (gemtrkr_3_peak_x_height[k]>gemtrkr3_xamp_max) {
        gemtrkr3_xamp_max=gemtrkr_3_peak_x_height[k];
        gemtrkr3_xch_max=gemtrkr_3_peak_pos_x[k];
      }
    }
    for (ULong64_t k=0; k<mmg1_idx_y; k++) {
      mmg1_peak_y->Fill(mmg1_peak_pos_y[k]);
      hmmg1_peak_y_height->Fill(mmg1_peak_y_height[k]);
      if (electron_tag) hmmg1_peak_y_height_el->Fill(mmg1_peak_y_height[k]);
      else if (pion_tag) hmmg1_peak_y_height_pi->Fill(mmg1_peak_y_height[k]);
      if (electron_tag && mmg1_peak_y_height[k]>mmg1_el_yamp_max) {
        mmg1_el_yamp_max=mmg1_peak_y_height[k];
        mmg1_el_ych_max=mmg1_peak_pos_y[k];
      } else if (pion_tag && mmg1_peak_y_height[k]>mmg1_pi_yamp_max) {
        mmg1_pi_yamp_max=mmg1_peak_y_height[k];
        mmg1_pi_ych_max=mmg1_peak_pos_y[k];
      }
      for (ULong64_t j=0; j<gt_1_idx_y; j++) {
        hgemtrkr_1_mmg1->Fill(mmg1_peak_pos_y[k], gemtrkr_1_peak_pos_y[j]);
      }
    }
    for (ULong64_t k=0; k<gem_idx_y; k++) {
      gem_peak_y->Fill(gem_peak_pos_y[k]);
      hgem_peak_y_height->Fill(gem_peak_y_height[k]);
      if (electron_tag) hgem_peak_y_height_el->Fill(gem_peak_y_height[k]);
      else if (pion_tag) hgem_peak_y_height_pi->Fill(gem_peak_y_height[k]);
      if (electron_tag && gem_peak_y_height[k]>gem_el_yamp_max) {
        gem_el_yamp_max=gem_peak_y_height[k];
        gem_el_ych_max=gem_peak_pos_y[k];
      } else if (pion_tag && gem_peak_y_height[k]>gem_pi_yamp_max) {
        gem_pi_yamp_max=gem_peak_y_height[k];
        gem_pi_ych_max=gem_peak_pos_y[k];
      }
      for (ULong64_t j=0; j<mmg1_idx_y; j++) {
        gem_mmg1_doubleY->Fill(gem_peak_pos_y[k], mmg1_peak_pos_y[j]);
      }
      for (ULong64_t j=0; j<gt_1_idx_y; j++) {
        hgemtrkr_1_gem->Fill(gem_peak_pos_y[k], gemtrkr_1_peak_pos_y[j]);
      }
    }
    
    if (gemtrkr1_xch_max>0.) hgemtrkr_1_max_xch->Fill(gemtrkr1_xch_max);
    if (gemtrkr1_xamp_max>0.) hgemtrkr_1_max_xamp->Fill(gemtrkr1_xamp_max);
    if (gemtrkr2_xch_max>0.) hgemtrkr_2_max_xch->Fill(gemtrkr2_xch_max);
    if (gemtrkr2_xamp_max>0.) hgemtrkr_2_max_xamp->Fill(gemtrkr2_xamp_max);
    if (gemtrkr3_xch_max>0.) hgemtrkr_3_max_xch->Fill(gemtrkr3_xch_max);
    if (gemtrkr3_xamp_max>0.) hgemtrkr_3_max_xamp->Fill(gemtrkr3_xamp_max);
    if (gemtrkr1_xch_max>0. && gemtrkr1_ych_max>0.) hgemtrkr_1_max_xy->Fill(gemtrkr1_xch_max, gemtrkr1_ych_max);
    if (gemtrkr2_xch_max>0. && gemtrkr2_ych_max>0.) hgemtrkr_2_max_xy->Fill(gemtrkr2_xch_max, gemtrkr2_ych_max);
    if (gemtrkr3_xch_max>0. && gemtrkr3_ych_max>0.) hgemtrkr_3_max_xy->Fill(gemtrkr3_xch_max, gemtrkr3_ych_max);
    
    //=============== END SRS Data Processing & Correlations ==============
    
    //==================================================================================================
    //                    Process Fa125  Pulse  data
    //==================================================================================================
    
    #if 1
      
      f125_fit->Reset();
      mmg1_f125_fit->Reset();
      double chi2cc_gem=-999., chi2cc_mmg1=-999.;
      double integral_gem=0., integral_mmg1=0.;
      #if (USE_PULSE>0)
        hevtk->Reset();
        hevtck->Reset();
      #endif
      double gem_ampmax=-1., mmg1_ampmax=-1., gem_xchmax=-1, mmg1_xchmax=-1;
      double gem_el_amp_max=0., gem_el_chan_max=-1.;
      double gem_pi_amp_max=0., gem_pi_chan_max=-1.;
      double mmg1_el_amp_max=0., mmg1_el_chan_max=-1.;
      double mmg1_pi_amp_max=0., mmg1_pi_chan_max=-1.;
      int gem_timemax=0, mmg1_timemax=0;
      int gem_trk_hit=0, mmg1_trk_hit=0;
      ULong64_t gem_idx_x=0, mmg1_idx_x=0;
      double gem_pos_x[f125_pulse_count], mmg1_pos_x[f125_pulse_count], gem_amp_x[f125_pulse_count], mmg1_amp_x[f125_pulse_count];
      
      if (gt_3_idx_x>0 && gt_2_idx_x>0) { //--External tracking condition
        
        Count("trk_hit");
        float x1=gemtrkr3_xch_max;
        float x2=gemtrkr2_xch_max;
        float a=(x2-x1)/(z2-z1);
        float b=((x1)*z2-(x2)*z1)/(z2-z1);
        float gem_extr = a*zgem+b;
        float mmg1_extr = a*zmmg1+b;
        if (electron_tag) {
          f125_el_tracker_hits->Fill(gem_extr);
          mmg1_f125_el_tracker_hits->Fill(mmg1_extr);
        } else if (pion_tag) {
          f125_pi_tracker_hits->Fill(gem_extr);
          mmg1_f125_pi_tracker_hits->Fill(mmg1_extr);
        }
        bool match = false, match_mmg1 = false;
        
        for (ULong64_t i=0; i<f125_pulse_count; i++) { //--- Fadc125 Pulse Loop
          
        	float peak_amp = f125_pulse_peak_amp->at(i);
        	float ped = f125_pulse_pedestal->at(i);
        	if (0 > ped || ped > 200 ) ped = 100;
        	float amp=peak_amp-ped;
        	if (amp<0) amp=0;
        	float time=f125_pulse_peak_time->at(i);
        	int fADCSlot = f125_pulse_slot->at(i);
        	int fADCChan = f125_pulse_channel->at(i);
        	
        	int gemChan = GetGEMChan(fADCChan, fADCSlot);
        	int mmg1Chan = GetMMG1Chan(fADCChan, fADCSlot, RunNum);
        	
        	if (gemChan>-1 && amp>GEM_THRESH) {
            float gemChan_x = gemChan*0.4+3.2; // to [mm]
            float gem_correction = -2.7931 + (gemChan_x)*0.02845;
            gem_residuals->Fill((gemChan_x-gem_extr));
            gem_residualscorr->Fill((gemChan_x-gem_extr)-gem_correction);
            gem_residual_ch->Fill(gemChan_x, (gemChan_x-gem_extr-gem_correction));
            gem_trk_hit=0;
            if (abs(gemChan_x-gem_extr-gem_correction)<10) { //within 10 mm
              Count ("gem_trk_hit");
              if (!match) {
                if (electron_tag) {
                  f125_el_tracker_eff->Fill(gem_extr);
                } else if (pion_tag) {
                 f125_pi_tracker_eff->Fill(gem_extr);
                }
                match = true;
              }
              gem_trk_hit++;
            } //-- END within 10mm
      	  }
      	  if (amp>MMG1_THRESH && mmg1Chan>-1) {
            float mmg1Chan_x = mmg1Chan*0.4+3.2; //-- to [mm]
            float mmg1_correction = -2.8703 + (mmg1Chan_x)*-0.000274;
            mmg1_residuals->Fill((mmg1Chan_x-mmg1_extr));
            mmg1_residualscorr->Fill((mmg1Chan_x-mmg1_extr)-mmg1_correction);
            mmg1_residual_ch->Fill(mmg1Chan_x, (mmg1Chan_x-mmg1_extr-mmg1_correction));
            mmg1_trk_hit=0;
            if (abs(mmg1Chan_x-mmg1_extr-mmg1_correction)<10) { //within 10 mm
              Count ("mmg1_trk_hit");
              if (!match_mmg1) {
                if (electron_tag) {
                  mmg1_f125_el_tracker_eff->Fill(mmg1_extr);
                } else if (pion_tag) {
                  mmg1_f125_pi_tracker_eff->Fill(mmg1_extr);
                }
                match_mmg1 = true;
              }
              mmg1_trk_hit++;
            } //-- END within 10mm
      	  }
      	} //--- end Fa125 Pulse Loop ---
      } //-- End external tracker condition !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      //===================================================================
      //                    Chi^2 Fit Calculation
      //===================================================================
      
      char f125Title[80]; sprintf(f125Title,"GEM-TRD:  Event=%lld Run=%d; z pos [time *8ns]; y pos [ch #]",jentry,RunNum);
      f125_fit->SetTitle(f125Title);
      if (f125_fit->GetEntries()!=0) {
        std::pair<Double_t, Double_t> fitResult  = TrkFit(f125_fit,fx1,"fx1",1);
        chi2cc_gem = fitResult.first;
        integral_gem = fitResult.second;
      }
      double a0 = fx1.GetParameter(0);
      double a1 = fx1.GetParameter(1);
  
      char mmg1f125Title[80]; sprintf(mmg1f125Title,"MMG1-TRD:  Event=%lld Run=%d; z pos [time *8ns]; y pos [ch #]",jentry,RunNum);
      mmg1_f125_fit->SetTitle(mmg1f125Title);
      if (mmg1_f125_fit->GetEntries()!=0) {
        std::pair<Double_t, Double_t> fitResult  = TrkFit(mmg1_f125_fit,fx2,"fx2",1);
        chi2cc_mmg1 = fitResult.first;
        integral_mmg1 = fitResult.second;
      }
      double a0_mmg1 = fx2.GetParameter(0);
      double a1_mmg1 = fx2.GetParameter(1);
      
      //==============================================
      //        Second fa125 pulse loop
      //==============================================
      
      for (ULong64_t i=0; i<f125_pulse_count; i++) { //--- Fadc125 Pulse Loop
        
        float peak_amp = f125_pulse_peak_amp->at(i);
        float ped = f125_pulse_pedestal->at(i);
       	if (0 > ped || ped > 200 ) ped = 100;
       	float amp=peak_amp-ped;
       	if (amp<0) amp=0;
       	float time=f125_pulse_peak_time->at(i);
       	int fADCSlot = f125_pulse_slot->at(i);
       	int fADCChan = f125_pulse_channel->at(i);
       	
        int gemChan = GetGEMChan(fADCChan, fADCSlot);
        float gemChan_x = gemChan*0.4+3.2; // to [mm]
       	int mmg1Chan = GetMMG1Chan(fADCChan, fADCSlot, RunNum);
       	float mmg1Chan_x = mmg1Chan*0.4+3.2; // to [mm]
        
       	if (gemChan>-1 && amp>GEM_THRESH) {
          if (electron_tag && gem_el_amp_max<amp) {
            gem_el_amp_max=amp;
            gem_el_chan_max=gemChan_x;
          } else if (pion_tag && gem_pi_amp_max<amp) {
            gem_pi_amp_max=amp;
            gem_pi_chan_max=gemChan_x;
          }
        }
        if (amp>MMG1_THRESH && mmg1Chan>-1) {
          if (electron_tag && mmg1_el_amp_max<amp) {
            mmg1_el_amp_max=amp;
            mmg1_el_chan_max=mmg1Chan_x;
          } else if (pion_tag && mmg1_pi_amp_max<amp) {
            mmg1_pi_amp_max=amp;
            mmg1_pi_chan_max=mmg1Chan_x;
          }
        }
      } //--END second f125 pulse loop
      
      //==============================================
      //        Third fa125 pulse loop
      //==============================================
      
      for (ULong64_t i=0; i<f125_pulse_count; i++) {
        
      	float peak_amp = f125_pulse_peak_amp->at(i);
      	float ped = f125_pulse_pedestal->at(i);
      	if (0 > ped || ped > 200 ) ped = 100;
      	float amp=peak_amp-ped;
      	float time=f125_pulse_peak_time->at(i);
      	int fADCSlot = f125_pulse_slot->at(i);
      	int fADCChan = f125_pulse_channel->at(i);
      	
      	int gemChan = GetGEMChan(fADCChan, fADCSlot);
      	double gemChan_x = gemChan*0.4 + 3.2;
      	int mmg1Chan = GetMMG1Chan(fADCChan, fADCSlot, RunNum);
      	double mmg1Chan_x = mmg1Chan*0.4 + 3.2;
        
      	if (amp<0) amp=0;
      	
      	if (gemChan>-1 && amp>GEM_THRESH) {
          f125_fit->Fill(time,gemChan,amp);
          #if (USE_PULSE>0)
            //------------ for Clusters -------------
            int TimeWindowStart = 45;
            int time0=int(time)-TimeWindowStart;
            if ( 0 < time0 && time0 < 130 ) { // --- drop early and late hits ---
              hevtck->SetBinContent(time0,gemChan,amp/10);
              hevtk->SetBinContent(time0,gemChan,amp/10);
            }
          #endif
    	    if (gem_ampmax<amp) {
            gem_ampmax=amp;
            gem_timemax=time;
            gem_xchmax=gemChan;
          }
          if (55.<=time && time<=142. && (abs(mmg1_el_chan_max-(gemChan_x-1.35))<2.5 || abs(mmg1_pi_chan_max-(gemChan_x-1.35))<2.5)) {
            gem_pos_x[i] = gemChan_x;
            gem_amp_x[i] = amp;
            gem_idx_x++;
          } else { gem_pos_x[i]=-1000.; gem_amp_x[i]=-1000.; gem_idx_x++; }
          if (electron_tag && abs(mmg1_el_chan_max-(gemChan_x-1.35))<2.5) { //-- Within 2.5mm
            hgemPulseDiff_el->Fill(mmg1_el_chan_max-(gemChan_x-1.35));
      	    f125_el->Fill(amp);
            f125_el_amp2d->Fill(time,gemChan,amp);
      	    f125_el_clu2d->Fill(time,gemChan,1.);
      	    gem_xpos.push_back(gemChan);
      	    gem_dedx.push_back(amp);
      	    gem_zpos.push_back(time);
      	    gem_parID.push_back(1);
      	    gem_nhit++;
      	    gem_zHist->Fill(time, amp);
          } else if (pion_tag && abs(mmg1_pi_chan_max-(gemChan_x-1.35))<2.5) { //-- Within 2.5mm
            hgemPulseDiff_pi->Fill(mmg1_pi_chan_max-(gemChan_x-1.35));
            f125_pi->Fill(amp);
            f125_pi_amp2d->Fill(time,gemChan,amp);
            f125_pi_clu2d->Fill(time,gemChan,1.);
            gem_xpos.push_back(gemChan);
            gem_dedx.push_back(amp);
            gem_zpos.push_back(time);
            gem_parID.push_back(0);
            gem_nhit++;
            gem_zHist->Fill(time, amp);
          }
    	  }
    	  if (amp>MMG1_THRESH && mmg1Chan>-1) {
          mmg1_f125_fit->Fill(time,mmg1Chan,amp);
    	    if (mmg1_ampmax<amp) {
            mmg1_ampmax=amp;
            mmg1_timemax=time;
            mmg1_xchmax=mmg1Chan;
          }
          if (45.<=time && time<=200. && (abs(gem_el_chan_max-mmg1Chan_x-1.35)<2.5 || abs(gem_pi_chan_max-mmg1Chan_x-1.35)<2.5)) {
            mmg1_pos_x[i] = mmg1Chan_x;
            mmg1_amp_x[i] = amp;
            mmg1_idx_x++;
          } else { mmg1_pos_x[i]=-1000.; mmg1_amp_x[i]=-1000.; mmg1_idx_x++; }
          if (electron_tag && abs(gem_el_chan_max-mmg1Chan_x-1.35)<2.5) {//-- Within 2.5mm
            hmmg1PulseDiff_el->Fill(gem_el_chan_max-mmg1Chan_x-1.35);
      	    mmg1_f125_el_amp2d->Fill(time,mmg1Chan,amp);
      	    mmg1_f125_el->Fill(amp);
      	    mmg1_f125_el_clu2d->Fill(time,mmg1Chan,1.);
      	    mmg1_xpos.push_back(mmg1Chan);
      	    mmg1_dedx.push_back(amp);
      	    mmg1_zpos.push_back(time);
      	    mmg1_parID.push_back(1);
      	    mmg1_nhit++;
      	    mmg1_zHist->Fill(time, amp);
          } else if (pion_tag && abs(gem_pi_chan_max-mmg1Chan_x-1.35)<2.5) {//-- Within 2.5mm
            hmmg1PulseDiff_pi->Fill(gem_pi_chan_max-mmg1Chan_x-1.35);
            mmg1_f125_pi_amp2d->Fill(time,mmg1Chan,amp);
            mmg1_f125_pi->Fill(amp);
            mmg1_f125_pi_clu2d->Fill(time,mmg1Chan,1.);
            mmg1_xpos.push_back(mmg1Chan);
            mmg1_dedx.push_back(amp);
            mmg1_zpos.push_back(time);
            mmg1_parID.push_back(0);
            mmg1_nhit++;
            mmg1_zHist->Fill(time, amp);
          }
    	  }
    	} //--- end Fa125 Pulse Loop ---
      
      if (gem_el_chan_max>0.) hchan_g_el->Fill(gem_el_chan_max);
      if (gem_pi_chan_max>0.) hchan_g_pi->Fill(gem_pi_chan_max);
      if (mmg1_el_chan_max>0.) hchan_m_el->Fill(mmg1_el_chan_max);
      if (mmg1_pi_chan_max>0.) hchan_m_pi->Fill(mmg1_pi_chan_max);
      
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
      else if (pion_tag==1) {
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
      gem_amp_max=gem_ampmax;
      gem_time_max=gem_timemax;
      gem_xch_max=gem_xchmax;
      gem_chi2cc.push_back(chi2cc_gem);
      gem_integral.push_back(integral_gem);
      mmg1_amp_max=mmg1_ampmax;
      mmg1_time_max=mmg1_timemax;
      mmg1_xch_max=mmg1_xchmax;
      mmg1_chi2cc.push_back(chi2cc_mmg1);
      mmg1_integral.push_back(integral_mmg1);
      
      if (abs(gem_el_chan_max-(mmg1_xchmax*0.4+3.2)-1.35)<2.5 || abs(gem_pi_chan_max-(mmg1_xchmax*0.4+3.2)-1.35)<2.5 || abs(mmg1_el_chan_max-((gem_xchmax*0.4+3.2)-1.35))<2.5 || abs(mmg1_pi_chan_max-((gem_xchmax*0.4+3.2)-1.35))<2.5) gem_mmg1_max_xcorr->Fill((gem_xchmax*0.4+3.2), (mmg1_xchmax*0.4+3.2));
      
      for (int i=1; i<21; i++) {
        gem_zHist_vect.push_back(gem_zHist->GetBinContent(i));
        mmg1_zHist_vect.push_back(mmg1_zHist->GetBinContent(i));
      }
      
      if (gem_nhit>0.) hgem_nhits->Fill(gem_nhit);
      if (mmg1_nhit>0.) hmmg1_nhits->Fill(mmg1_nhit);
      if (gt1_nhit>0.) hgt1_nhits->Fill(gt1_nhit);
      if (gt2_nhit>0.) hgt2_nhits->Fill(gt2_nhit);
      if (gt3_nhit>0.) hgt3_nhits->Fill(gt3_nhit);
      
      for (ULong64_t i=0; i<gt_1_idx_x; i++) {
        if (gemtrkr_1_peak_pos_x[i]>1.) gem_gt1_xcorr->Fill((gem_xchmax*0.4+3.2), gemtrkr_1_peak_pos_x[i]);
      }
      for (ULong64_t i=0; i<gt_2_idx_x; i++) {
        if (gemtrkr_2_peak_pos_x[i]>1.) gem_gt2_xcorr->Fill((gem_xchmax*0.4+3.2), gemtrkr_2_peak_pos_x[i]);
      }
      for (ULong64_t i=0; i<gt_3_idx_x; i++) {
        if (gemtrkr_3_peak_pos_x[i]>1.) gem_gt3_xcorr->Fill((gem_xchmax*0.4+3.2), gemtrkr_3_peak_pos_x[i]);
      }
      
      for (ULong64_t j=0; j<gem_idx_x; j++) {
        if (gem_pos_x[j]>1.) {
          for (ULong64_t i=0; i<gem_idx_y; i++) {
            if (gem_peak_pos_y[i]>0.) hgem_xy->Fill(gem_pos_x[j], gem_peak_pos_y[i]);
          }
          for (ULong64_t i=0; i<mmg1_idx_x; i++) {
            if (mmg1_pos_x[i]>1.) gem_mmg1_xcorr->Fill(gem_pos_x[j], mmg1_pos_x[i]);
          }
/*          for (ULong64_t i=0; i<gt_1_idx_x; i++) {
            if (gemtrkr_1_peak_pos_x[i]>1.) gem_gt1_xcorr->Fill(gem_pos_x[j], gemtrkr_1_peak_pos_x[i]);
          }
          for (ULong64_t i=0; i<gt_2_idx_x; i++) {
            if (gemtrkr_2_peak_pos_x[i]>1.) gem_gt2_xcorr->Fill(gem_pos_x[j], gemtrkr_2_peak_pos_x[i]);
          }
          for (ULong64_t i=0; i<gt_3_idx_x; i++) {
            if (gemtrkr_3_peak_pos_x[i]>1.) gem_gt3_xcorr->Fill(gem_pos_x[j], gemtrkr_3_peak_pos_x[i]);
          }*/
        }
      }
      for (ULong64_t i=0; i<gem_idx_y; i++) {
        if (gem_peak_pos_y[i]>0.) {
          for (ULong64_t j=0; j<mmg1_idx_y; j++) {
            if (mmg1_peak_pos_y[j]>0.) gem_mmg1_ycorr->Fill(gem_peak_pos_y[i], mmg1_peak_pos_y[j]);
          }
        }
      }
      for (ULong64_t j=0; j<mmg1_idx_y; j++) {
        if (mmg1_peak_pos_y[j]>0.) {
          for (ULong64_t i=0; i<mmg1_idx_x; i++) {
            if (mmg1_pos_x[i]>1.) hmmg1_xy->Fill(mmg1_pos_x[i], mmg1_peak_pos_y[j]);
          }
        }
      }
      for (ULong64_t j=0; j<mmg1_idx_x; j++) {
        if (mmg1_pos_x[j]>1.) {
          //for (ULong64_t i=0; i<gt_1_idx_x; i++) {
            //if (gemtrkr_1_peak_pos_x[i]>1.) mmg1_gt1_xcorr->Fill(mmg1_pos_x[j], gemtrkr_1_peak_pos_x[i]);
          //}
          //for (ULong64_t i=0; i<gt_2_idx_x; i++) {
          //  if (gemtrkr_2_peak_pos_x[i]>1.) mmg1_gt2_xcorr->Fill(mmg1_pos_x[j], gemtrkr_2_peak_pos_x[i]);
         // }
         // for (ULong64_t i=0; i<gt_3_idx_x; i++) {
         //   if (gemtrkr_3_peak_pos_x[i]>1.) mmg1_gt3_xcorr->Fill(mmg1_pos_x[j], gemtrkr_3_peak_pos_x[i]);
         // }
        }
      }
      
      for (ULong64_t i=0; i<gt_1_idx_x; i++) {
        if (gemtrkr_1_peak_pos_x[i]>1.) mmg1_gt1_xcorr->Fill((mmg1_xchmax*0.4+3.2), gemtrkr_1_peak_pos_x[i]);
      }
      for (ULong64_t i=0; i<gt_2_idx_x; i++) {
        if (gemtrkr_2_peak_pos_x[i]>1.) mmg1_gt2_xcorr->Fill((mmg1_xchmax*0.4+3.2), gemtrkr_2_peak_pos_x[i]);
      }
      for (ULong64_t i=0; i<gt_3_idx_x; i++) {
        if (gemtrkr_3_peak_pos_x[i]>1.) mmg1_gt3_xcorr->Fill((mmg1_xchmax*0.4+3.2), gemtrkr_3_peak_pos_x[i]);
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
      mhevtc->Reset();
      mhevtf->Reset();
      hevt->Reset();
      hevtc->Reset();
      hevtf->Reset();
      
      for (ULong64_t i=0;i<f125_wraw_count; i++) { // --- fadc125 channels loop
        
        int fadc_window = f125_wraw_samples_count->at(i);
        int fADCSlot = f125_wraw_slot->at(i);
        int fADCChan = f125_wraw_channel->at(i);
        int gemChan = GetGEMChan(fADCChan, fADCSlot);
        int mmg1Chan = GetMMG1Chan(fADCChan, fADCSlot, RunNum);
        double DEDX_THR = GEM_THRESH, mDEDX_THR = MMG1_THRESH;
        int TimeWindowStart = 45;
        int TimeMin = 0;
        int TimeMax = 140;
        
        //--ped calculation for f125 raw data
        int nped = 0, ped = 100;
        double ped_sum = 0.;
        for (int si=TimeWindowStart-15; si<TimeWindowStart; si++) {
          int ped_samp = f125_wraw_samples->at(f125_wraw_samples_index->at(i)+si);
          ped_sum += ped_samp;
          nped++;
        }
        ped = ped_sum / nped;
        if (0. > ped || ped > 200 ) ped = 100;

        for (int si=0; si<fadc_window; si++) {
          int time=si;
          int adc = f125_wraw_samples->at(f125_wraw_samples_index->at(i)+si);
          adc = adc - ped;
          if (adc>4090) printf("!!!!!!!!!!!!!!!!!!!!!! ADC 125 overflow: %d \n",adc);
          if (adc>DEDX_THR) {
            if (gemChan>-1) {
              if (RunNum>5284.) { time-=(TimeWindowStart+35); } //--Second Xe Bottle
              else { time-=TimeWindowStart; }
              ///////////////if ( TimeMin > time || time > TimeMax ) continue; // --- drop early and late hits ---
              
              hevtc->SetBinContent(100-time,gemChan+1,adc/100.);
              hevt->SetBinContent(100-time,gemChan+1,adc/100.);
            }
          }
          if (adc > mDEDX_THR) {
            if (mmg1Chan>-1) {
              time-=(TimeWindowStart+70); //-35
              mhevtc->SetBinContent(100-time,mmg1Chan+1,adc/100.);
              mhevt->SetBinContent(100-time,mmg1Chan+1,adc/100.);
            }
          }
        } // --  end of samples loop
      } // -- end of fadc125 raw channels loop
      
      #ifdef SHOW_EVT_DISPLAY
          #if (USE_PULSE>0)
            c2->cd(1); hevtk->Draw("box");
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
        //--GEM-TRD
        float clust_Xmax[MAX_CLUST];
        float clust_Zmax[MAX_CLUST];
        float clust_Emax[MAX_CLUST];
        
        float clust_Xpos[MAX_CLUST];
        float clust_Zpos[MAX_CLUST];
        float clust_dEdx[MAX_CLUST];
        float clust_Size[MAX_CLUST];
        float clust_Width[MAX_CLUST][3];  // y1, y2, dy ; strips
        float clust_Length[MAX_CLUST][3]; // x1, x2, dx ; time
        float hits_Xpos[500];
        float hits_Zpos[500];
        float hits_dEdx[500];
        float hits_Size[MAX_CLUST];
        float hits_Width[MAX_CLUST];  // y1, y2, dy ; strips
        float hits_Length[MAX_CLUST]; // x1, x2, dx ; time
        //--MMG1TRD
        float mmg1_clust_Xmax[MAX_CLUST];
        float mmg1_clust_Zmax[MAX_CLUST];
        float mmg1_clust_Emax[MAX_CLUST];
        
        float mmg1_clust_Xpos[MAX_CLUST];
        float mmg1_clust_Zpos[MAX_CLUST];
        float mmg1_clust_dEdx[MAX_CLUST];
        float mmg1_clust_Size[MAX_CLUST];
        float mmg1_clust_Width[MAX_CLUST][3];  // y1, y2, dy ; strips
        float mmg1_clust_Length[MAX_CLUST][3]; // x1, x2, dx ; time
        float mmg1_hits_Xpos[500];
        float mmg1_hits_Zpos[500];
        float mmg1_hits_dEdx[500];
        float mmg1_hits_Size[MAX_CLUST];
        float mmg1_hits_Width[MAX_CLUST];  // y1, y2, dy ; strips
        float mmg1_hits_Length[MAX_CLUST]; // x1, x2, dx ; time
        
        for (int k=0; k<MAX_CLUST; k++) {
	        clust_Xpos[k]=0; clust_Zpos[k]=0; clust_dEdx[k]=0;  clust_Size[k]=0;
	        clust_Xmax[k]=0; clust_Zmax[k]=0; clust_Emax[k]=0;
          clust_Width[k][0]=999999;   	clust_Width[k][1]=-999999;   	clust_Width[k][2]=0;
          clust_Length[k][0]=999999;  	clust_Length[k][1]=-999999;  	clust_Length[k][2]=0;
          
          mmg1_clust_Xpos[k]=0; mmg1_clust_Zpos[k]=0; mmg1_clust_dEdx[k]=0;  mmg1_clust_Size[k]=0;
          mmg1_clust_Xmax[k]=0; mmg1_clust_Zmax[k]=0; mmg1_clust_Emax[k]=0;
          mmg1_clust_Width[k][0]=999999;     mmg1_clust_Width[k][1]=-999999;    mmg1_clust_Width[k][2]=0;
          mmg1_clust_Length[k][0]=999999;    mmg1_clust_Length[k][1]=-999999;   mmg1_clust_Length[k][2]=0;
        }
        int nclust=0, mmg1_nclust=0;
        #if (USE_PULSE>0)
          TH2F* hp = hevtk; // -- hevtk and hevtck should be same bin size
          TH2F* hpc = hevtck;
        #else
          TH2F* hp = hevt; // -- hevt and hevtc should be same bin size
          TH2F* hpc = hevtc;
          TH2F* hmp = mhevt;
          TH2F* hmpc = mhevtc;
        #endif
        //--GEM
        int nx=hp->GetNbinsX();    int ny=hp->GetNbinsY();
        double xmi=hp->GetXaxis()->GetBinLowEdge(1);     double xma=hp->GetXaxis()->GetBinUpEdge(nx);
        double ymi=hp->GetYaxis()->GetBinLowEdge(1);     double yma=hp->GetYaxis()->GetBinUpEdge(ny);
        double binx = (xma-xmi)/nx;      double biny = (yma-ymi)/ny;
        //--MMG1
        int nmx=hmp->GetNbinsX();    int nmy=hmp->GetNbinsY();
        double xmmi=hmp->GetXaxis()->GetBinLowEdge(1);   double xmma=hmp->GetXaxis()->GetBinUpEdge(nmx);
        double ymmi=hmp->GetYaxis()->GetBinLowEdge(1);   double ymma=hmp->GetYaxis()->GetBinUpEdge(nmy);
        double binmx = (xmma-xmmi)/nmx;      double binmy = (ymma-ymmi)/nmy;
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
        
        for (int iy=0; iy<ny; iy++) {  //-------------------- Clustering Loop (GEMTRD) ------------------------------------
          for (int ix=0; ix<nx; ix++) {
            double c1 = hpc->GetBinContent(ix+1,iy+1);                    // energy
            double x1=double(ix)/double(nx)*(xma-xmi)+xmi+binx/2.;    // drift time
            double y1=double(iy)/double(ny)*(yma-ymi)+ymi+biny/2.;    // X strip
            if (c1<THR2) continue;
            if (nclust==0) {
	            clust_Xpos[nclust]=y1; clust_Zpos[nclust]=x1;  clust_dEdx[nclust]=c1;  clust_Size[nclust]=1;
	            clust_Xmax[nclust]=y1; clust_Zmax[nclust]=x1;  clust_Emax[nclust]=c1;
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
	              if (c1>clust_Emax[k]) {
		              clust_Xmax[k]=y1;
		              clust_Zmax[k]=x1;
		              clust_Emax[k]=c1;
	              }
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
	            clust_Xpos[nclust]=y1; clust_Zpos[nclust]=x1;  clust_dEdx[nclust]=c1;  clust_Size[nclust]=1;
	            clust_Xmax[nclust]=y1; clust_Zmax[nclust]=x1;  clust_Emax[nclust]=c1;
              clust_Width[nclust][0]=y1;   	clust_Width[nclust][1]=y1;   	clust_Width[nclust][2]=0;
              clust_Length[nclust][0]=x1;  	clust_Length[nclust][1]=x1;  	clust_Length[nclust][2]=0;
              nclust++;
            }
          }
        } //---------------------- End Clustering Loop (GEMTRD) ------------------------------
        
        for (int iy=0; iy<nmy; iy++) {  //-------------------- Clustering Loop (MMG1TRD) ------------------------------------
          for (int ix=0; ix<nmx; ix++) {
            double c1 = hmpc->GetBinContent(ix+1,iy+1);                       // energy
            double x1=double(ix)/double(nmx)*(xmma-xmmi)+xmmi+binmx/2.;   // drift time
            double y1=double(iy)/double(nmy)*(ymma-ymmi)+ymmi+binmy/2.;        // X strip
            if (c1<THR2) continue;
            if (mmg1_nclust==0) {
              mmg1_clust_Xpos[mmg1_nclust]=y1; mmg1_clust_Zpos[mmg1_nclust]=x1;  mmg1_clust_dEdx[mmg1_nclust]=c1;  mmg1_clust_Size[mmg1_nclust]=1;
              mmg1_clust_Xmax[mmg1_nclust]=y1; mmg1_clust_Zmax[mmg1_nclust]=x1;  mmg1_clust_Emax[mmg1_nclust]=c1;
              mmg1_clust_Width[mmg1_nclust][0]=y1;    mmg1_clust_Width[mmg1_nclust][1]=y1;    mmg1_clust_Width[mmg1_nclust][2]=0;
              mmg1_clust_Length[mmg1_nclust][0]=x1;   mmg1_clust_Length[mmg1_nclust][1]=x1;   mmg1_clust_Length[mmg1_nclust][2]=0;
              mmg1_nclust++; continue;
            }
            int mmg1_added=0;
            for (int k=0; k<mmg1_nclust; k++) {
              double dist=sqrt(pow((y1-mmg1_clust_Xpos[k]),2.)+pow((x1-mmg1_clust_Zpos[k]),2.)); //--- dist hit to clusters
              if (dist<CL_DIST) {
                mmg1_clust_Xpos[k]=(y1*c1+mmg1_clust_Xpos[k]*mmg1_clust_dEdx[k])/(c1+mmg1_clust_dEdx[k]);  //--  new X pos
                mmg1_clust_Zpos[k]=(x1*c1+mmg1_clust_Zpos[k]*mmg1_clust_dEdx[k])/(c1+mmg1_clust_dEdx[k]);  //--  new Z pos
                if (c1>clust_Emax[k]) {
		              mmg1_clust_Xmax[k]=y1;
		              mmg1_clust_Zmax[k]=x1;
		              mmg1_clust_Emax[k]=c1;
	              }
                mmg1_clust_dEdx[k]=c1+mmg1_clust_dEdx[k];  // new dEdx
                mmg1_clust_Size[k]=1+mmg1_clust_Size[k];  // clust size in pixels
                if (y1<mmg1_clust_Width[k][0]) mmg1_clust_Width[k][0]=y1; if (y1>mmg1_clust_Width[k][1]) mmg1_clust_Width[k][1]=y1; mmg1_clust_Width[k][2]=mmg1_clust_Width[k][1]-mmg1_clust_Width[k][0];
                if (x1<mmg1_clust_Length[k][0]) mmg1_clust_Length[k][0]=x1;if (x1>mmg1_clust_Length[k][1]) mmg1_clust_Length[k][1]=x1; mmg1_clust_Length[k][2]=mmg1_clust_Length[k][1]-mmg1_clust_Length[k][0];
                hmpc->SetBinContent(ix,iy,k+1.);
                mmg1_added=1; break;
              }
            }
            if (mmg1_added==0) {
              if (mmg1_nclust+1>=MAX_CLUST) continue;
              mmg1_clust_Xpos[mmg1_nclust]=y1; mmg1_clust_Zpos[mmg1_nclust]=x1;  mmg1_clust_dEdx[mmg1_nclust]=c1;  mmg1_clust_Size[mmg1_nclust]=1;
              mmg1_clust_Xmax[mmg1_nclust]=y1; mmg1_clust_Zmax[mmg1_nclust]=x1;  mmg1_clust_Emax[mmg1_nclust]=c1;
              mmg1_clust_Width[mmg1_nclust][0]=y1;    mmg1_clust_Width[mmg1_nclust][1]=y1;    mmg1_clust_Width[mmg1_nclust][2]=0;
              mmg1_clust_Length[mmg1_nclust][0]=x1;   mmg1_clust_Length[mmg1_nclust][1]=x1;   mmg1_clust_Length[mmg1_nclust][2]=0;
              mmg1_nclust++;
            }
          }
        } //---------------------- End Clustering Loop (MMG1TRD) ------------------------------
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
        
        double maxClust_dEdx=0., totalClust_dEdx=0.;
        int ii=0;
        #ifdef VERBOSE
          printf("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
          printf("                Xpos   Ypos   Zpos       E    Width  Length   Size \n");
        #endif
        for (int k=0; k<nclust; k++) {
          #ifdef VERBOSE
            if (k<30) printf("%2d Clust(%2d): %6.1f %6.1f %8.1f %6.2f %6.2f %8.1f  ",k,k+1,clust_Xpos[k],clust_Zpos[k],clust_dEdx[k],clust_Width[k][2],clust_Length[k][2],clust_Size[k]);
          #endif
          //-------------  Cluster Filter (GEMTRD) -----------------
          if ((clust_Size[k] >= MinClustSize && zStart < clust_Zpos[k] && clust_Zpos[k] < zEnd && clust_Width[k][2]>MinClustWidth) || clust_Length[k][2]<MaxClustLength) {
            if (abs(mmg1_el_chan_max-(clust_Xpos[k]+3.2-1.35))<5. || abs(mmg1_pi_chan_max-(clust_Xpos[k]+3.2-1.35))<5.) {
              //FILL HIT DIFF HISTO
              if (electron_tag) hgemClusterDiff_el->Fill(mmg1_el_chan_max-(clust_Xpos[k]+3.2-1.35));
              else if (pion_tag) hgemClusterDiff_pi->Fill(mmg1_pi_chan_max-(clust_Xpos[k]+3.2-1.35));
        	    #if (USE_MAXPOS>0)
	              hits_Xpos[ii]=clust_Xmax[k];
	              hits_Zpos[ii]=clust_Zmax[k];
              #else
                hits_Xpos[ii]=clust_Xpos[k];
        	      hits_Zpos[ii]=clust_Zpos[k];
              #endif
        	    hits_dEdx[ii]=clust_dEdx[k];
              hits_Width[ii]=clust_Width[k][2];
              hits_Length[ii]=clust_Length[k][2];
        	    ii++;
        	    if (clust_dEdx[k]>maxClust_dEdx) maxClust_dEdx=clust_dEdx[k];
              totalClust_dEdx+=clust_dEdx[k];
              #ifdef VERBOSE
                if (k<30) printf("\n");
              #endif
            }
        	} else {
          #ifdef VERBOSE
            if (k<30) printf(" <--- skip \n");
          #endif
          }
        }
        int nhits=ii;
        clu_dedx_tot=totalClust_dEdx;
        
        double maxClust_m_dEdx=0., totalClust_m_dEdx=0.;
        int mmg1_ii=0;
        for (int k=0; k<mmg1_nclust; k++) {
          //-------------  Cluster Filter (MMG1TRD) -----------------
          if ((mmg1_clust_Size[k] >= MinClustSize && zStart < mmg1_clust_Zpos[k] && mmg1_clust_Zpos[k] < zEnd && mmg1_clust_Width[k][2]>MinClustWidth) || mmg1_clust_Length[k][2]<MaxClustLength) {
            if (abs(gem_el_chan_max-(mmg1_clust_Xpos[k]+3.2)-1.35)<5. || abs(gem_pi_chan_max-(mmg1_clust_Xpos[k]+3.2)-1.35)<5.) {
              //FILL HIT DIFF HISTO
              if (electron_tag) hmmg1ClusterDiff_el->Fill(gem_el_chan_max-(mmg1_clust_Xpos[k]+3.2)-1.35);
              else if (pion_tag) hmmg1ClusterDiff_pi->Fill(gem_pi_chan_max-(mmg1_clust_Xpos[k]+3.2)-1.35);
              mmg1_hits_Xpos[mmg1_ii]=mmg1_clust_Xpos[k];
              mmg1_hits_Zpos[mmg1_ii]=mmg1_clust_Zpos[k];
              mmg1_hits_dEdx[mmg1_ii]=mmg1_clust_dEdx[k];
              mmg1_hits_Width[mmg1_ii]=mmg1_clust_Width[k][2];
              mmg1_hits_Length[mmg1_ii]=mmg1_clust_Length[k][2];
              mmg1_ii++;
              if (mmg1_clust_dEdx[k]>maxClust_m_dEdx) maxClust_m_dEdx=mmg1_clust_dEdx[k];
              totalClust_m_dEdx+=mmg1_clust_dEdx[k];
            }
          }
        }
        int mmg1_nhits=mmg1_ii;
        mmg1_clu_dedx_tot=totalClust_m_dEdx;
        // ----------------------- end hist dist clustering ---------------------------------
        
        //=================================== Draw HITS and CLUST  ============================================
        #ifdef SHOW_EVTbyEVT
          char hevtTitle[80]; sprintf(hevtTitle,"GEM-TRD:  Event=%lld Run=%d e=%d #pi=%d; z pos [mm]; y pos [mm]",jentry,RunNum,electron_tag,pion_tag);
          hevt->SetTitle(hevtTitle);
          #if (USE_PULSE>0)
            hevtk->SetTitle(hevtTitle);
          #endif
          char mhevtTitle[80]; sprintf(mhevtTitle,"MMG1-TRD:  Event=%lld Run=%d e=%d #pi=%d; z pos [mm]; y pos [mm]",jentry,RunNum,electron_tag,pion_tag);
          mhevt->SetTitle(mhevtTitle);
          #ifdef VERBOSE
            printf("hits_SIZE=%d  Clust size = %d \n",nhits,nclust);
          #endif
          	c2->cd(1); gPad->Modified(); gPad->Update();
          	int COLMAP[]={1,2,3,4,6,5};
          	int pmt=22, pmt0 = 20; // PM type
            int max2draw = nclust;
            for (int i=0; i<max2draw; i++) {
              #if (USE_MAXPOS>0)
	              TMarker m = TMarker(clust_Zmax[i],clust_Xmax[i],pmt);
              #else
          	    TMarker m = TMarker(clust_Zpos[i], clust_Xpos[i], pmt);
              #endif
          	  int tcol=2;
          	  if (clust_Size[i]<MinClustSize) pmt=22; else pmt=pmt0;
          	  int mcol = COLMAP[tcol-1];   m.SetMarkerColor(mcol);   m.SetMarkerStyle(pmt);
          	  m.SetMarkerSize(0.7+clust_dEdx[i]/300);
          	  m.DrawClone();
          	}
            gPad->Modified(); gPad->Update();
            
            //--MMG1TRD
            c2->cd(6); gPad->Modified(); gPad->Update();
            int mCOLMAP[]={1,2,3,4,6,5};
            int mpmt=22, mpmt0 = 20; // PM type
            int mmax2draw = mmg1_nclust;
            for (int i=0; i<mmax2draw; i++) {
              #if (USE_MAXPOS>0)
	              TMarker m = TMarker(mmg1_clust_Zmax[i],mmg1_clust_Xmax[i], mpmt);
              #else
                TMarker m = TMarker(mmg1_clust_Zpos[i], mmg1_clust_Xpos[i], mpmt);
              #endif
              int tcol=2;
              if (mmg1_clust_Size[i]<MinClustSize) mpmt=22; else mpmt=mpmt0;
              int mcol = mCOLMAP[tcol-1];   m.SetMarkerColor(mcol);   m.SetMarkerStyle(mpmt);
              m.SetMarkerSize(0.7+mmg1_clust_dEdx[i]/300);
              m.DrawClone();
            }
            gPad->Modified(); gPad->Update();
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
          //-- MMG1
          std::vector<int> mmg1_tracks(mmg1_nhits, 0);
          std::vector<float> mmg1_Xcl;
          std::vector<float> mmg1_Zcl;
          mmg1_Xcl.clear();
          mmg1_Zcl.clear();
          for (int n=0; n<mmg1_nhits; n++) {
            mmg1_Xcl.push_back(mmg1_hits_Xpos[n]);
            mmg1_Zcl.push_back(mmg1_hits_Zpos[n]);
          }
          doPattern(mmg1_Xcl, mmg1_Zcl, mmg1_tracks);  //---- call GNN ---
          #ifdef VERBOSE
            printf("**> End Model simulation \n"); //===================================================
          #endif
          #ifdef SHOW_EVTbyEVT
              c2->cd(2); gPad->Modified(); gPad->Update();
              int COLMAP2[]={1,2,3,4,6,5};
              for(ULong64_t i=0; i<tracks.size(); i++) {
                #ifdef VERBOSE
                  if (i<30) printf("i=%d trk=%d |  %8.2f,%8.2f\n",i, tracks[i], Xcl[i], Zcl[i]);
              	#endif
                TMarker m = TMarker(hits_Zpos[i], hits_Xpos[i], 24);
              	int tcol = min(tracks[i], 6);
              	int mcol = COLMAP2[tcol-1];   m.SetMarkerColor(mcol);   m.SetMarkerStyle(41);     m.SetMarkerSize(1.5);
            	  m.DrawClone();  gPad->Modified(); gPad->Update();
              }
              //-- MMG1
              c2->cd(7); gPad->Modified(); gPad->Update();
              int mCOLMAP2[]={1,2,3,4,6,5};
              for(ULong64_t i=0; i<mmg1_tracks.size(); i++) {
                TMarker m = TMarker(mmg1_hits_Zpos[i], mmg1_hits_Xpos[i], 24);
                int tcol = min(mmg1_tracks[i], 6);
                int mcol = mCOLMAP2[tcol-1];   m.SetMarkerColor(mcol);   m.SetMarkerStyle(41);     m.SetMarkerSize(1.5);
                m.DrawClone();  gPad->Modified(); gPad->Update();
              }
              #ifdef VERBOSE
                printf("\n\n");
                printf("**> End Cluster Plot \n");
              #endif
          #endif
          //--------------------------------------------------
          //----           Track fitting                 -----
          //--------------------------------------------------
          
          #ifdef VERBOSE
            printf("==> GNN: tracks sort  : trk_siz=%ld \r\n", tracks.size());
          #endif
          //-----------------   tracks sorting -------------
          std::vector<std::vector<float>> TRACKS;
          TRACKS.resize(nhits);
          //std::vector<float> hit_coord(2, 0);
          std::vector<int>  TRACKS_N(nhits, 0);
          for (int i=0; i<nhits; i++)  { TRACKS_N[i] = 0;  }
          //std::vector<float> xz(2,0);
          for (int i2=0; i2<nhits; i2++) {
          	int num =  tracks[i2];
          	int num2 = std::max(0, std::min(num, nhits - 1));
            #ifdef VERBOSE
              if (i2<20) printf("==> lstm3:track sort i=%d  : num=%d(%d) x=%f z=%f \n", i2, num, num2,  Xcl[i2],Zcl[i2]);
          	#endif
            //xz[0]=Xcl[i2];
            //xz[1]=Zcl[i2];
          	TRACKS[num2].push_back(Xcl[i2]);
            TRACKS[num2].push_back(Zcl[i2]);
          	TRACKS_N[num2]++;
          }
          //-- MMG1
          std::vector<std::vector<float>> mmg1_TRACKS;
          mmg1_TRACKS.resize(mmg1_nhits);
          //std::vector<float> mmg1_hit_coord(2, 0);
          std::vector<int>  mmg1_TRACKS_N(mmg1_nhits, 0);
          for (int i=0; i<mmg1_nhits; i++)  { mmg1_TRACKS_N[i] = 0;  }
          //std::vector<float> mmg1_xz(2, 0);
          for (int i2=0; i2<mmg1_nhits; i2++) {
            int num =  mmg1_tracks[i2];
            int num2 = std::max(0, std::min(num, mmg1_nhits - 1));
            //mmg1_xz[0]=mmg1_Xcl[i2];
            //mmg1_xz[1]=mmg1_Zcl[i2];
            mmg1_TRACKS[num2].push_back(mmg1_Xcl[i2]);
            mmg1_TRACKS[num2].push_back(mmg1_Zcl[i2]);
            mmg1_TRACKS_N[num2]++;
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
          //------------------ end tracks sorting --------------------
          
          #if (USE_FIT==1)
            //-----------------------------------
            //---       linear fitting        ---
            //-----------------------------------
            static TMultiGraph *mg;
            if (mg != NULL ) delete mg;
            mg = new TMultiGraph();
            int NTRACKS=0;
            int MIN_HITS=2;
            Double_t p0, p1;
            
            for (int i2 = 1; i2 < nhits; i2++) {  //-- GEM tracks loop; zero track -> noise
              
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
              	gErrorIgnoreLevel = kBreak; // Suppress warning messages from empty fit data
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
            }  //-- end GEM tracks loop --
            
            //-- MMG1
            static TMultiGraph *mmg1_mg;
            if (mmg1_mg != NULL ) delete mmg1_mg;
            mmg1_mg = new TMultiGraph();
            int mmg1_NTRACKS=0;
            int mmg1_MIN_HITS=2;
            Double_t mp0, mp1;
            
            for (int i2=1; i2<mmg1_nhits; i2++) {  //-- MMG1 tracks loop; zero track -> noise
              
              if (mmg1_TRACKS_N[i2]<mmg1_MIN_HITS) continue;   //---- select 2 (x,z) or more hits on track ----
              std::vector<Double_t> mmg1_x;
              std::vector<Double_t> mmg1_y;
              for (int i3=0; i3<(int)mmg1_TRACKS[i2].size(); i3+=2) {
                mmg1_x.push_back(mmg1_TRACKS[i2].at(i3+1));
                mmg1_y.push_back(mmg1_TRACKS[i2].at(i3));
              }
              #ifdef SHOW_EVTbyEVT
                gErrorIgnoreLevel = kBreak; // Suppress warning messages from empty fit data
                TGraph *mmg1_g = new TGraph(mmg1_TRACKS_N[i2], &mmg1_x[0], &mmg1_y[0]);  mmg1_g->SetMarkerStyle(21); mmg1_g->SetMarkerColor(i2);
                TF1 *mmg1_f = new TF1("mmg1_f", "[1] * x + [0]");
                mmg1_g->Fit(mmg1_f,"Q");
                //  --- get fit parameters ---
                TF1 *mmg1_ffunc = mmg1_g->GetFunction("mmg1_f");
                mp0 = mmg1_ffunc->GetParameter(0);
                mp1 = mmg1_ffunc->GetParameter(1);
                Double_t mmg1_chi2x_nn = mmg1_ffunc->GetChisquare();
                Double_t mmg1_Ndfx_nn = mmg1_ffunc->GetNDF();
                double mmg1_chi2nn = mmg1_chi2x_nn/mmg1_Ndfx_nn;
                mmg1_mg->Add(mmg1_g, "p");
              #endif
              mmg1_NTRACKS++;
            }  //-- end MMG1 tracks loop --
            
            #ifdef SHOW_EVTbyEVT
              if (NTRACKS!=1) continue;  // --- skip event ----
              //if (nhits<3)    continue;  // --- skip event ----
              if (gem_trk_hit<1) continue;
                char mgTitle[80]; sprintf(mgTitle,"GEM ML-FPGA response, #Tracks=%d; z pos [mm]; y pos [mm]",NTRACKS);
                mg->SetTitle(mgTitle);
                c2->cd(3); mg->Draw("AP");
                mg->GetXaxis()->SetLimits(Xmin,Xmax);
                mg->SetMinimum(Ymin);
                mg->SetMaximum(Ymax);
                gPad->Modified(); gPad->Update();
                //-- MMG1
                char mmg1_mgTitle[80]; sprintf(mmg1_mgTitle,"MMG1 ML-FPGA response, #Tracks=%d; z pos [mm]; y pos [mm]",mmg1_NTRACKS);
                mmg1_mg->SetTitle(mmg1_mgTitle);
                c2->cd(8); mmg1_mg->Draw("AP");
                mmg1_mg->GetXaxis()->SetLimits(Xmin,Xmax);
                mmg1_mg->SetMinimum(Ymin);
                mmg1_mg->SetMaximum(Ymax);
                gPad->Modified(); gPad->Update();
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
            printf("All done, click middle of canvas ...\n");
            #ifdef VERBOSE
              printf(" a0=%f a1=%f (%f deg)  fx1(150)=%f chi2cc_gem=%f  \n",a0,a1,a1/3.1415*180.,fx1.Eval(150.),chi2cc_gem);
            #endif
            if (electron_tag || pion_tag) c2->cd(1); gPad->WaitPrimitive();
        #endif
      #endif   // --- End if USE_CLUST>0 ---
    #endif   //=======================  End Fa125 RAW process Loop  =====================================
    //============ END GEMTRD Pattern Recognition Tracking ==================
    
    if (NTRACKS==1) Count("singleTRK");
    if (NTRACKS>1) Count("multTRK");
    if (NTRACKS==1 && electron_tag) Count("snTRKel");
    if (NTRACKS==1 && pion_tag) Count("snTRKpi");
    if (mmg1_NTRACKS==1) Count("m_singleTRK");
    if (mmg1_NTRACKS>1) Count("m_multTRK");
    if (mmg1_NTRACKS==1 && electron_tag) Count("m_snTRKel");
    if (mmg1_NTRACKS==1 && pion_tag) Count("m_snTRKpi");
    if (electron_tag) {
      if (maxClust_dEdx!=0.) hClusterMaxdEdx_el->Fill(maxClust_dEdx);
      if (totalClust_dEdx!=0.) hClusterTotaldEdx_el->Fill(totalClust_dEdx);
      if (maxClust_m_dEdx!=0.) hmmg1ClusterMaxdEdx_el->Fill(maxClust_m_dEdx);
      if (totalClust_m_dEdx!=0.) hmmg1ClusterTotaldEdx_el->Fill(totalClust_m_dEdx);
    } else if (pion_tag) {
      if (maxClust_dEdx!=0.) hClusterMaxdEdx_pi->Fill(maxClust_dEdx);
      if (totalClust_dEdx!=0.) hClusterTotaldEdx_pi->Fill(totalClust_dEdx);
      if (maxClust_m_dEdx!=0.) hmmg1ClusterMaxdEdx_pi->Fill(maxClust_m_dEdx);
      if (totalClust_m_dEdx!=0.) hmmg1ClusterTotaldEdx_pi->Fill(totalClust_m_dEdx);
    }
    //=====================================================================================
    //===                Fill Root TTree Hits                                            ===
    //=====================================================================================
    
    #ifdef SAVE_TRACK_HITS
      gem_nclu=nhits;
      for (int n=0; n<nhits; n++) {
        clu_xpos.push_back(hits_Xpos[n]);
        clu_zpos.push_back(hits_Zpos[n]);
        clu_dedx.push_back(hits_dEdx[n]);
        clu_width.push_back(hits_Width[n]);
        if (hits_dEdx[n] > clu_dedx_max) {
          clu_dedx_max=hits_dEdx[n];
          clu_xpos_max=hits_Xpos[n];
          clu_zpos_max=hits_Zpos[n];
          clu_width_max=hits_Width[n];
        }
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
      
      mmg1_nclu=mmg1_nhits;
      for (int n=0; n<mmg1_nhits; n++) {
        mmg1_clu_xpos.push_back(mmg1_hits_Xpos[n]);
        mmg1_clu_zpos.push_back(mmg1_hits_Zpos[n]);
        mmg1_clu_dedx.push_back(mmg1_hits_dEdx[n]);
        mmg1_clu_width.push_back(mmg1_hits_Width[n]);
        if (mmg1_hits_dEdx[n] > mmg1_clu_dedx_max) {
          mmg1_clu_dedx_max=mmg1_hits_dEdx[n];
          mmg1_clu_xpos_max=mmg1_hits_Xpos[n];
          mmg1_clu_zpos_max=mmg1_hits_Zpos[n];
          mmg1_clu_width_max=mmg1_hits_Width[n];
        }
      }
      
      if (gem_nhit>0) EVENT_VECT_GEM->Fill();
      if (mmg1_nhit>0)EVENT_VECT_MMG1->Fill();
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
  TH1D *f125_drift = f125_el_amp2d->ProjectionX("f125_drift",100,135);
  TH1D *f125_drift_c = new TH1D(*f125_drift); f125_drift_c->SetStats(0);
  double gemDriftScale = 1./f125_drift_c->GetEntries();
  f125_drift_c->Scale(gemDriftScale);
  TH1D *mmg1_f125_drift = mmg1_f125_el_amp2d->ProjectionX("mmg1_f125_drift",95,130);
  TH1D *mmg1_f125_drift_c = new TH1D(*mmg1_f125_drift); mmg1_f125_drift_c->SetStats(0);
  double mmg1DriftScale = 1./mmg1_f125_drift_c->GetEntries();
  mmg1_f125_drift_c->Scale(mmg1DriftScale);
  f125_drift_c->SetLineColor(4);  HistList->Add(f125_drift_c);
  mmg1_f125_drift_c->SetLineColor(2); HistList->Add(mmg1_f125_drift_c);
  f125_drift_c->GetXaxis()->SetTitle("Time Response (8ns)");
  f125_drift_c->GetYaxis()->SetTitle("1. / nEntries");
  f125_drift_c->SetTitle("Drift Time Distribution (Electron Response)");
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
  #if ANALYZE_MERGED
    char rootFileName[256]; sprintf(rootFileName, "RootOutput/cern24/merged/Run_%06d_%0dEntries_Output.root",RunNum,nEntries);
  #else
    char rootFileName[256]; sprintf(rootFileName, "RootOutput/cern24/Run_%06d_Output.root",RunNum);
  #endif
  fOut = new TFile(rootFileName, "RECREATE");
  fOut->cd();
  cout<<"Writing Output File: "<<rootFileName<<endl;
  HistList->Write("HistDQM", TObject::kSingleKey);
  c4->Write();
  fOut->Close();
  delete fOut;
  
  //=====================================================================================
  //===                 S A V E   T R A C K   H I T   T T R E E S                    ====
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
  //===                 P L O T     H I S T O G R A M S                               ===
  //=====================================================================================
  #if ANALYZE_MERGED
    const char *OutputDir="RootOutput/cern24/merged";
  #else
    const char *OutputDir="RootOutput/cern24";
  #endif
  #ifdef SAVE_PDF
    char ctit[120];
    #if ANALYZE_MERGED
      sprintf(G_DIR,"%s/Run_%06d_%06dEntries",OutputDir,RunNum,nEntries);
    #else
      sprintf(G_DIR,"%s/Run_%06d",OutputDir,RunNum);
    #endif
    sprintf(ctit,"File=%s",G_DIR);
    bool COMPACT=false;
    TCanvas *cc;
    int nxd=3;
    int nyd=5;
    char pdfname[120];  sprintf(pdfname,"%s_evdisp.pdf",G_DIR);  //c0->Print(pdfname);
    
    //---------------------  page 1 --------------------
    htitle(" Count ");   // if (!COMPACT) cc=NextPlot(0,0);
    nxd=2; nyd=4;
    cc=NextPlot(nxd,nyd);  gPad->SetLogy(); hcount->Draw();
    cc=NextPlot(nxd,nyd);  hgem_nhits->Draw();
    cc=NextPlot(nxd,nyd);  hmmg1_nhits->Draw();
    cc=NextPlot(nxd,nyd);  hgt1_nhits->Draw();
    cc=NextPlot(nxd,nyd);  hgt2_nhits->Draw();
    cc=NextPlot(nxd,nyd);  hgt3_nhits->Draw();
    
    //---------------------  page 2a --------------------
    htitle(" TRD Correlations  ");   if (!COMPACT) cc=NextPlot(0,0);
    cc=NextPlot(nxd,nyd);  hgemtrkr_1_gem->Draw("colz");
    cc=NextPlot(nxd,nyd);  hgemtrkr_1_mmg1->Draw("colz");
    cc=NextPlot(nxd,nyd);  gem_mmg1_xcorr->Draw("colz");
    cc=NextPlot(nxd,nyd);  gem_mmg1_max_xcorr->Draw("colz");
    cc=NextPlot(nxd,nyd);  gem_mmg1_ycorr->Draw("colz");
    cc=NextPlot(nxd,nyd);  hgem_xy->Draw("colz");
    cc=NextPlot(nxd,nyd);  hmmg1_xy->Draw("colz");
    
    htitle(" X Correlations  ");   if (!COMPACT) cc=NextPlot(0,0);
    cc=NextPlot(nxd,nyd);  gem_gt1_xcorr->Draw("colz");
    cc=NextPlot(nxd,nyd);  gem_gt2_xcorr->Draw("colz");
    cc=NextPlot(nxd,nyd);  gem_gt3_xcorr->Draw("colz");
    cc=NextPlot(nxd,nyd);  mmg1_gt1_xcorr->Draw("colz");
    cc=NextPlot(nxd,nyd);  mmg1_gt2_xcorr->Draw("colz");
    cc=NextPlot(nxd,nyd);  mmg1_gt3_xcorr->Draw("colz");
    
    //---------------------  page 2a --------------------
    htitle(" SRS GEM-TRKR 1  ");   if (!COMPACT) cc=NextPlot(0,0);
    //nxd=2; nyd=4;
    cc=NextPlot(nxd,nyd); hgemtrkr_1_peak_xy->Draw("colz");
    cc=NextPlot(nxd,nyd); hgemtrkr_1_max_xy->Draw("colz");
    cc=NextPlot(nxd,nyd); hgemtrkr_1_peak_x->Draw();
    cc=NextPlot(nxd,nyd); hgemtrkr_1_peak_y->Draw();
    cc=NextPlot(nxd,nyd); hgemtrkr_1_peak_x_height->Draw();
    cc=NextPlot(nxd,nyd); hgemtrkr_1_peak_y_height->Draw();
    cc=NextPlot(nxd,nyd); hgemtrkr_1_max_xch->Draw();
    cc=NextPlot(nxd,nyd); hgemtrkr_1_max_xamp->Draw();
    
    htitle(" SRS GEM-TRKR 2  ");   if (!COMPACT) cc=NextPlot(0,0);
    cc=NextPlot(nxd,nyd); hgemtrkr_2_peak_xy->Draw("colz");
    cc=NextPlot(nxd,nyd); hgemtrkr_2_max_xy->Draw("colz");
    cc=NextPlot(nxd,nyd); hgemtrkr_2_peak_x->Draw();
    cc=NextPlot(nxd,nyd); hgemtrkr_2_peak_y->Draw();
    cc=NextPlot(nxd,nyd); hgemtrkr_2_peak_x_height->Draw();
    cc=NextPlot(nxd,nyd); hgemtrkr_2_peak_y_height->Draw();
    cc=NextPlot(nxd,nyd); hgemtrkr_2_max_xch->Draw();
    cc=NextPlot(nxd,nyd); hgemtrkr_2_max_xamp->Draw();
    
    htitle(" SRS GEM-TRKR 3  ");   if (!COMPACT) cc=NextPlot(0,0);
    //nxd=2; nyd=3;
    cc=NextPlot(nxd,nyd); hgemtrkr_3_peak_xy->Draw("colz");
    cc=NextPlot(nxd,nyd); hgemtrkr_3_max_xy->Draw("colz");
    cc=NextPlot(nxd,nyd); hgemtrkr_3_peak_x->Draw();
    cc=NextPlot(nxd,nyd); hgemtrkr_3_peak_y->Draw();
    cc=NextPlot(nxd,nyd); hgemtrkr_3_peak_x_height->Draw();
    cc=NextPlot(nxd,nyd); hgemtrkr_3_peak_y_height->Draw();
    cc=NextPlot(nxd,nyd); hgemtrkr_3_max_xch->Draw();
    cc=NextPlot(nxd,nyd); hgemtrkr_3_max_xamp->Draw();
    
    htitle(" SRS GEM / MMG1  ");   if (!COMPACT) cc=NextPlot(0,0);
    //nxd=2; nyd=3;
    cc=NextPlot(nxd,nyd); mmg1_peak_y->Draw();
    cc=NextPlot(nxd,nyd); hmmg1_peak_y_height->Draw();
    cc=NextPlot(nxd,nyd); hmmg1_peak_y_height_el->Draw();
    cc=NextPlot(nxd,nyd); hmmg1_peak_y_height_pi->Draw();
    cc=NextPlot(nxd,nyd); gem_peak_y->Draw();
    cc=NextPlot(nxd,nyd); hgem_peak_y_height->Draw();
    cc=NextPlot(nxd,nyd); hgem_peak_y_height_el->Draw();
    cc=NextPlot(nxd,nyd); hgem_peak_y_height_pi->Draw();
    //TBox fbox(xbc1,ybc1,xbc2,ybc2);  //---- draw box cut ---
    //fbox.Draw("same");
    //fbox.SetLineColor(kRed);
    //fbox.SetFillStyle(0);
    //fbox.SetLineWidth(1);
    
    //htitle(" SRS GEM / MMG1  ");   if (!COMPACT) cc=NextPlot(0,0);
    //nxd=2; nyd=4;
    //cc=NextPlot(nxd,nyd); gem_mmg1_doubleY->Draw("colz");
    //cc=NextPlot(nxd,nyd); gem_mmg1_doubleX->Draw("colz");
    
    //--------------------- new page --------------------
    htitle(" fADC125 Max Hit Distributions with Tracking Restriction");   if (!COMPACT) cc=NextPlot(0,0);
    //nxd=2; nyd=3;
    cc=NextPlot(nxd,nyd);  hchan_g_el->Draw();
    cc=NextPlot(nxd,nyd);  hchan_g_pi->Draw();
    cc=NextPlot(nxd,nyd);  hchan_m_el->Draw();
    cc=NextPlot(nxd,nyd);  hchan_m_pi->Draw();
    
    //--------------------- new page --------------------
    htitle(" fADC125 Raw (Clustering) Distributions ");   if (!COMPACT) cc=NextPlot(0,0);
    //nxd=2; nyd=4;
    cc=NextPlot(nxd,nyd);  hClusterMaxdEdx_el->Draw();
    cc=NextPlot(nxd,nyd);  hClusterMaxdEdx_pi->Draw();
    cc=NextPlot(nxd,nyd);  hClusterTotaldEdx_el->Draw();
    cc=NextPlot(nxd,nyd);  hClusterTotaldEdx_pi->Draw();
    cc=NextPlot(nxd,nyd);  hmmg1ClusterMaxdEdx_el->Draw();
    cc=NextPlot(nxd,nyd);  hmmg1ClusterMaxdEdx_pi->Draw();
    cc=NextPlot(nxd,nyd);  hmmg1ClusterTotaldEdx_el->Draw();
    cc=NextPlot(nxd,nyd);  hmmg1ClusterTotaldEdx_pi->Draw();
    
    //--------------------- new page --------------------
    htitle(" fADC125 Raw (Clustering) Track Differences ");   if (!COMPACT) cc=NextPlot(0,0);
    cc=NextPlot(nxd,nyd);  hgemPulseDiff_el->Draw();
    cc=NextPlot(nxd,nyd);  hgemPulseDiff_pi->Draw();
    cc=NextPlot(nxd,nyd);  hmmg1PulseDiff_el->Draw();
    cc=NextPlot(nxd,nyd);  hmmg1PulseDiff_pi->Draw();
    cc=NextPlot(nxd,nyd);  hgemClusterDiff_el->Draw();
    cc=NextPlot(nxd,nyd);  hgemClusterDiff_pi->Draw();
    cc=NextPlot(nxd,nyd);  hmmg1ClusterDiff_el->Draw();
    cc=NextPlot(nxd,nyd);  hmmg1ClusterDiff_pi->Draw();
    
   //---------------------  page 3 --------------------
    htitle("  TRD (fa125) Amp Distributions ");    if (!COMPACT) cc=NextPlot(0,0);
    //nxd=2; nyd=3;
    cc=NextPlot(nxd,nyd);   gPad->SetLogy();  f125_el->Draw();
    cc=NextPlot(nxd,nyd);   gPad->SetLogy();  mmg1_f125_el->Draw();
    cc=NextPlot(nxd,nyd);   gPad->SetLogy();  f125_el_max->Draw();
    cc=NextPlot(nxd,nyd);   gPad->SetLogy();  f125_pi_max->Draw();
    cc=NextPlot(nxd,nyd);   gPad->SetLogy();  mmg1_f125_el_max->Draw();
    cc=NextPlot(nxd,nyd);   gPad->SetLogy();  mmg1_f125_pi_max->Draw();
    
   //---------------------  page 3a --------------------
    htitle("  GEM-TRD (fa125) Amp 2D");    if (!COMPACT) cc=NextPlot(0,0);
    //nxd=2; nyd=3;
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
    //nxd=2; nyd=3;
    cc=NextPlot(nxd,nyd);   mmg1_f125_el_amp2d->Draw("colz");
    cc=NextPlot(nxd,nyd);   mmg1_f125_pi_amp2d->Draw("colz");
    cc=NextPlot(nxd,nyd);   mmg1_f125_el_clu2d->Draw("colz");
    cc=NextPlot(nxd,nyd);   mmg1_f125_pi_clu2d->Draw("colz");
   
   //---------------------  page 3a --------------------
    htitle("  External Tracking");    if (!COMPACT) cc=NextPlot(0,0);
    //nxd=2; nyd=4;
    cc=NextPlot(nxd,nyd);   f125_el_tracker_hits->Draw();
    cc=NextPlot(nxd,nyd);   mmg1_f125_el_tracker_hits->Draw();
    cc=NextPlot(nxd,nyd);   f125_el_tracker_eff->Draw();
    cc=NextPlot(nxd,nyd);   mmg1_f125_el_tracker_eff->Draw();
    cc=NextPlot(nxd,nyd);   f125_pi_tracker_hits->Draw();
    cc=NextPlot(nxd,nyd);   mmg1_f125_pi_tracker_hits->Draw();
    cc=NextPlot(nxd,nyd);   f125_pi_tracker_eff->Draw();
    cc=NextPlot(nxd,nyd);   mmg1_f125_pi_tracker_eff->Draw();
    
    htitle("  External Tracking");    if (!COMPACT) cc=NextPlot(0,0);
    nxd=2; nyd=3;
    cc=NextPlot(nxd,nyd);   if (f125_el_tracker_hits->GetEntries()!=0) { f125_el_tracker_hits->Sumw2(); f125_el_tracker_eff->Sumw2();  f125_el_tracker_eff->Divide(f125_el_tracker_hits);  f125_el_tracker_eff->Draw(); }
    cc=NextPlot(nxd,nyd);   if (mmg1_f125_el_tracker_hits->GetEntries()!=0) { mmg1_f125_el_tracker_hits->Sumw2(); mmg1_f125_el_tracker_eff->Sumw2(); mmg1_f125_el_tracker_eff->Divide(mmg1_f125_el_tracker_hits);  mmg1_f125_el_tracker_eff->Draw(); }
    cc=NextPlot(nxd,nyd);   if (f125_pi_tracker_hits->GetEntries()!=0) { f125_pi_tracker_hits->Sumw2(); f125_pi_tracker_eff->Sumw2();  f125_pi_tracker_eff->Divide(f125_pi_tracker_hits);  f125_pi_tracker_eff->Draw(); }
    cc=NextPlot(nxd,nyd);   if (mmg1_f125_pi_tracker_hits->GetEntries()!=0) { mmg1_f125_pi_tracker_hits->Sumw2(); mmg1_f125_pi_tracker_eff->Sumw2(); mmg1_f125_pi_tracker_eff->Divide(mmg1_f125_pi_tracker_hits);  mmg1_f125_pi_tracker_eff->Draw(); }
    
    htitle("  External Tracking");    if (!COMPACT) cc=NextPlot(0,0);
    cc=NextPlot(nxd,nyd);   gem_residuals->Draw();
    cc=NextPlot(nxd,nyd);   gem_residualscorr->Draw();
    cc=NextPlot(nxd,nyd);   mmg1_residuals->Draw();
    cc=NextPlot(nxd,nyd);   mmg1_residualscorr->Draw();
    cc=NextPlot(nxd,nyd);   gem_residual_ch->Draw("colz");
    cc=NextPlot(nxd,nyd);   mmg1_residual_ch->Draw("colz");
    
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
