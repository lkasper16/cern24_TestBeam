
# =============================
# clhep
# =============================
setenv CLHEP "/home/hdtrdops/jana4ml4fpga/jana4ml4fpga/miniconda/envs/ml4fpga"
setenv CLHEP_BASE_DIR "/home/hdtrdops/jana4ml4fpga/jana4ml4fpga/miniconda/envs/ml4fpga"
setenv CLHEP_INCLUDE_DIR "/home/hdtrdops/jana4ml4fpga/jana4ml4fpga/miniconda/envs/ml4fpga/include"
setenv CLHEP_LIB_DIR "/home/hdtrdops/jana4ml4fpga/jana4ml4fpga/miniconda/envs/ml4fpga/lib"

# Make sure PATH is set
if ( ! $?PATH ) then
    setenv PATH "/home/hdtrdops/jana4ml4fpga/jana4ml4fpga/miniconda/envs/ml4fpga/bin"
else
    setenv PATH "/home/hdtrdops/jana4ml4fpga/jana4ml4fpga/miniconda/envs/ml4fpga/bin":${PATH}
endif

# Make sure LD_LIBRARY_PATH is set
if ( ! $?LD_LIBRARY_PATH ) then
    setenv LD_LIBRARY_PATH "/home/hdtrdops/jana4ml4fpga/jana4ml4fpga/miniconda/envs/ml4fpga/lib"
else
    setenv LD_LIBRARY_PATH "/home/hdtrdops/jana4ml4fpga/jana4ml4fpga/miniconda/envs/ml4fpga/lib":${LD_LIBRARY_PATH}
endif

# =============================
# hepmc3
# =============================

# Make sure PATH is set
if ( ! $?PATH ) then
    setenv PATH "/home/hdtrdops/jana4ml4fpga/jana4ml4fpga/miniconda/envs/ml4fpga/bin"
else
    setenv PATH "/home/hdtrdops/jana4ml4fpga/jana4ml4fpga/miniconda/envs/ml4fpga/bin":${PATH}
endif
setenv HEPMC3_DIR "/home/hdtrdops/jana4ml4fpga/jana4ml4fpga/miniconda/envs/ml4fpga"

# Make sure CMAKE_PREFIX_PATH is set
if ( ! $?CMAKE_PREFIX_PATH ) then
    setenv CMAKE_PREFIX_PATH "/home/hdtrdops/jana4ml4fpga/jana4ml4fpga/miniconda/envs/ml4fpga/share/HepMC3/cmake"
else
    setenv CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH}:"/home/hdtrdops/jana4ml4fpga/jana4ml4fpga/miniconda/envs/ml4fpga/share/HepMC3/cmake"
endif

# =============================
# jana2
# =============================
setenv JANA_HOME "/home/hdtrdops/jana4ml4fpga/jana2/jana2-v2.1.0"

# Make sure JANA_PLUGIN_PATH is set
if ( ! $?JANA_PLUGIN_PATH ) then
    setenv JANA_PLUGIN_PATH "$JANA_HOME/plugins"
else
    setenv JANA_PLUGIN_PATH ${JANA_PLUGIN_PATH}:"$JANA_HOME/plugins"
endif

# Make sure PATH is set
if ( ! $?PATH ) then
    setenv PATH "$JANA_HOME/bin"
else
    setenv PATH "$JANA_HOME/bin":${PATH}
endif

# Make sure LD_LIBRARY_PATH is set
if ( ! $?LD_LIBRARY_PATH ) then
    setenv LD_LIBRARY_PATH "/home/hdtrdops/jana4ml4fpga/jana2/jana2-v2.1.0/lib"
else
    setenv LD_LIBRARY_PATH "/home/hdtrdops/jana4ml4fpga/jana2/jana2-v2.1.0/lib":${LD_LIBRARY_PATH}
endif

# Make sure CMAKE_PREFIX_PATH is set
if ( ! $?CMAKE_PREFIX_PATH ) then
    setenv CMAKE_PREFIX_PATH "/home/hdtrdops/jana4ml4fpga/jana2/jana2-v2.1.0/lib/cmake/JANA"
else
    setenv CMAKE_PREFIX_PATH "/home/hdtrdops/jana4ml4fpga/jana2/jana2-v2.1.0/lib/cmake/JANA":${CMAKE_PREFIX_PATH}
endif

# =============================
# jana4ml4fpga
# =============================
setenv jana4ml4fpga_HOME "/home/hdtrdops/jana4ml4fpga/jana4ml4fpga/jana4ml4fpga-main"

# Make sure JANA_PLUGIN_PATH is set
if ( ! $?JANA_PLUGIN_PATH ) then
    setenv JANA_PLUGIN_PATH "/home/hdtrdops/jana4ml4fpga/jana4ml4fpga/jana4ml4fpga-main/plugins"
else
    setenv JANA_PLUGIN_PATH "/home/hdtrdops/jana4ml4fpga/jana4ml4fpga/jana4ml4fpga-main/plugins":${JANA_PLUGIN_PATH}
endif

# Make sure PATH is set
if ( ! $?PATH ) then
    setenv PATH "/home/hdtrdops/jana4ml4fpga/jana4ml4fpga/jana4ml4fpga-main/bin"
else
    setenv PATH "/home/hdtrdops/jana4ml4fpga/jana4ml4fpga/jana4ml4fpga-main/bin":${PATH}
endif

# Make sure LD_LIBRARY_PATH is set
if ( ! $?LD_LIBRARY_PATH ) then
    setenv LD_LIBRARY_PATH "/home/hdtrdops/jana4ml4fpga/jana4ml4fpga/jana4ml4fpga-main/lib"
else
    setenv LD_LIBRARY_PATH "/home/hdtrdops/jana4ml4fpga/jana4ml4fpga/jana4ml4fpga-main/lib":${LD_LIBRARY_PATH}
endif

# =============================
# root
# =============================

# Make sure CMAKE_PREFIX_PATH is set
if ( ! $?CMAKE_PREFIX_PATH ) then
    setenv CMAKE_PREFIX_PATH "/home/hdtrdops/jana4ml4fpga/jana4ml4fpga/miniconda/envs/ml4fpga/cmake/"
else
    setenv CMAKE_PREFIX_PATH "/home/hdtrdops/jana4ml4fpga/jana4ml4fpga/miniconda/envs/ml4fpga/cmake/":${CMAKE_PREFIX_PATH}
endif

if (! $?ROOT_INSTALLED_BY_CONDA) then
   if ( -f /home/hdtrdops/jana4ml4fpga/jana4ml4fpga/miniconda/envs/ml4fpga/bin/thisroot.csh ) then
      source /home/hdtrdops/jana4ml4fpga/jana4ml4fpga/miniconda/envs/ml4fpga/bin/thisroot.csh
   endif
endif

