# cern24_TestBeam

## Repository for CERN July 2024 TRD Test Beam

################################  
### Single Run Analysis Workflow

Workflow on JLab Gluon compute nodes:   
Log in to Gluon nodes with JLab computing account (must have 2FA)
```
ssh -XY [$USERNAME]@scilogin.jlab.org  
Password: [Pin+OTP]  
ssh hallgw  
Password: [Pin+OTP]  
ssh gluon[100-150]  
Password: [JLab CUE]  
```
After cloning repo, make directory links   
```
./make_gluon_links.sh  
```
And use executable files for run analysis   
```
./trdclass_cern24.sh [$RUNNUMBER] [$MAXNUMBEROFEVENTS] [$FIRSTEVENT]  
./trd_mlp_cern.sh [$RUNNUMBER]
```
With this workflow, raw .evio data files have already been processed into .root files that now live in the `ROOT/` directory. These .root data files are analyzed with the `trdclass_cern24.C` analysis macro. Output from there is saved in the `RootOutput/cern24` directory. The output from this in the form of a .root TTree file is passed on to `trd_mlp_cern.C` where a rejection factor calculation is done for different particle efficiencies. Output from this NN macro is saved in the `mlpOutput/cern24` directory.  

###############################  
### Multi-Run Analysis Workflow

Workflow on JLab Gluon compute nodes:   
Log in to Gluon nodes with JLab computing account (must have 2FA)
```
ssh -XY [$USERNAME]@scilogin.jlab.org  
Password: [Pin+OTP]  
ssh hallgw  
Password: [Pin+OTP]  
ssh gluon[100-150]  
Password: [JLab CUE]  
```
After cloning repo, make directory links   
```
./make_gluon_links.sh  
```
Make directory to store multi-run merged files   
```
mkdir ROOT_MERGED/  
```
Use ROOT macro found [here](https://github.com/lkasper16/ROOT_macros/commit/076639ec3d75df0e5a0ff28e1b8d2efbba1501b6) to pass a .txt file containing a list of run data files desired to be merged into one root TChain file with the following naming convention: *eventsChain[RUNNUMBER]\_[NENTRIES]Entries\_[NTREES]Trees.root*   
```
cd ROOT_MERGED/
root -l -q 'multiFileMerge.C(“firstBottleDoubleFleece.txt”)'  
```
Once the new root file is generated, open the analysis macro and set the `ANALYZE_MERGED` flag to 1   
```
cd ../  
vi trdclass_cern24.C  
> #define ANALYZE_MERGED 1  
```
Then, use the executable file for run analysis   
```
./trdclass_cern24Merged.sh [$RUNNUMBER] [$NENTRIES] [$NTREES] [$MAXNUMBEROFEVENTS] [$FIRSTEVENT]  
```
To then execute the trd_mlp_cern24.C macro over the output and calculate rejection factors, first change the same flag to 1 inside this macro   
```
vi trd_mlp_cern24.C
> #define ANALYZE_MERGED 1  
```
... And use the executable file for processing   
```
./trd_mlp_cernMerged.sh [$RUNNUMBER] [$NENTRIES]  
```
With this workflow, raw .evio data files have already been processed into .root files that now live in the `ROOT/` directory. These .root data files are now merged together into one .root file that lives in the `ROOT_MERGED/` directory. This file is then analyzed with the `trdclass_cern24.C` analysis macro. Output from there is saved in the `RootOutput/cern24/merged/` directory. The output from this in the form of a .root TTree file is passed on to `trd_mlp_cern.C` where a rejection factor calculation is done for different particle efficiencies. Output from this NN macro is saved in the `mlpOutput/cern24/merged` directory.  

See JLab's [2FA documentation](https://jlab.servicenowservices.com/sp?id=kb_article_view&sysparm_article=KB0012313&sys_kb_id=a8caee091b990910a552ed3ce54bcbe3&spa=1.)  
See JLab analysis code repository: [trd_root](https://github.com/JeffersonLab/trd_root/tree/main)  
See JLab JANA repository: [JANA4ML4FPGA](https://github.com/JeffersonLab/JANA4ML4FPGA/tree/main)  

