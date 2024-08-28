# cern24_TestBeam

## Repository for CERN July 2024 TRD Test Beam

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
```
With this workflow, raw .evio data files have already been processed into .root files that now live in the `ROOT/` directory. These .root data files are analyzed with the `trdclass_cern24.C` analysis macro. Output from there is saved in the `RootOutput/cern24` directory.  

See JLab's [2FA documentation](https://jlab.servicenowservices.com/sp?id=kb_article_view&sysparm_article=KB0012313&sys_kb_id=a8caee091b990910a552ed3ce54bcbe3&spa=1.)  
See JLab analysis code repository: [trd_root](https://github.com/JeffersonLab/trd_root/tree/main)  
See JLab JANA repository: [JANA4ML4FPGA](https://github.com/JeffersonLab/JANA4ML4FPGA/tree/main)  

