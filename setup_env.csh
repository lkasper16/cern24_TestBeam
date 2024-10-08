
setenv ML4FPGA_TOP_DIR /gluonfs1/home/hdtrdops/soft_ml4fpga

# Start conda environment
source $ML4FPGA_TOP_DIR/miniconda/etc/profile.d/conda.csh
conda activate ml4fpga

# The path where edpm stores its JSon database and creates env files
setenv EDPM_DATA_PATH $ML4FPGA_TOP_DIR/edpm_data

# This tells EDPM not to generate source thisroot.sh
setenv ROOT_INSTALLED_BY_CONDA 1

# source environment generated by EDPM
# means ROOT and others
source $EDPM_DATA_PATH/env.csh

# =============================
# jana4ml4fpga
# =============================

# Make sure JANA_PLUGIN_PATH is set
if ( ! $?JANA_PLUGIN_PATH ) then
    setenv JANA_PLUGIN_PATH "$ML4FPGA_TOP_DIR/jana4ml4fpga/jana4ml4fpga-main/install/plugins"
else
    setenv JANA_PLUGIN_PATH "$ML4FPGA_TOP_DIR/jana4ml4fpga/jana4ml4fpga-main/install/plugins":${JANA_PLUGIN_PATH}
endif

# Make sure PATH is set
if ( ! $?PATH ) then
    setenv PATH "$ML4FPGA_TOP_DIR/jana4ml4fpga/jana4ml4fpga-main/install/bin"
else
    setenv PATH "$ML4FPGA_TOP_DIR/jana4ml4fpga/jana4ml4fpga-main/install/bin":${PATH}
endif

# Make sure LD_LIBRARY_PATH is set
if ( ! $?LD_LIBRARY_PATH ) then
    setenv LD_LIBRARY_PATH "$ML4FPGA_TOP_DIR/jana4ml4fpga/jana4ml4fpga-main/install/lib"
else
    setenv LD_LIBRARY_PATH "$ML4FPGA_TOP_DIR/jana4ml4fpga/jana4ml4fpga-main/install/lib":${LD_LIBRARY_PATH}
endif