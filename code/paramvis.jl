# ------------ GENERIC MODULES -------------------------------------------------
using PyPlot
import CSV
import JLD

# ------------ FLOW CODES ------------------------------------------------------

# GeometricTools https://github.com/byuflowlab/GeometricTools
import GeometricTools
gt = GeometricTools

# WFVisual https://github.com/byuflowlab/WFVisual.jl
# import WFVisual
reload("WFVisual")
wfv=WFVisual


# ------------ GLOBAL VARIABLES ------------------------------------------------
global module_path; module_path,_ = splitdir(@__FILE__);   # Path to this module
global extdrive_path = "/media/edoalvar/MyExtDrive/simulationdata5/" # Path to
                                                           # external drive

# ------------ HEADERS ---------------------------------------------------------
for module_name in ["functions", "generate_vtks"]
  include("paramvis_"*module_name*".jl")
end


function main()
    generate_vtks(; data_path="data/direct4",
                                save_path="temps/direct4_00",
                                prompt=true,
                                previs=false, paraview=false)
end
