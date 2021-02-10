from HFSSdrawpy import Body, Entity, Modeler
from HFSSdrawpy.interfaces.hfss_modeler import *
from HFSSdrawpy.libraries.base_elements import *
from HFSSdrawpy.utils import Vector, parse_entry

"""
testing script for waveport assignment and network data export

ATTENTION; script creates a new project in HFSS and analysis

generates a shielded microstrip line on sapphire including ground plane 
for terminal assignment

"""

###########################################################setup new project

# DrivenModal or DrivenTerminal
simtype = "DrivenModal"

desktop = get_desktop()
project = desktop.new_project()
project.make_active()
if simtype == "DrivenTerminal":
    design = project.new_dt_design("test")
if simtype == "DrivenModal":
    design = project.new_dm_design("test")

setup_name = "Setup"

modeler = "hfss"
pm = Modeler(modeler)

pm.is_litho = False
pm.is_hfss = True

chip_body = Body(pm, "chip")

# #############################################################Drawing

track = pm.set_variable("0.32mm")
sub_h = pm.set_variable("0.43mm")
cover_H = pm.set_variable("0.57mm")
MSL_length = pm.set_variable("1mm")
width = pm.set_variable("3mm")

# define substrate + MSL + GND
chip_subs = chip_body.box([-width / 2, 0, -sub_h], [width, MSL_length, sub_h], name="chip_subs")
chip_subs.assign_material("sapphire")

MSL = chip_body.rect([-track / 2, 0, 0], [track, MSL_length, 0], name="MSL")
MSL.assign_perfect_E("_perfE")

GND = chip_body.rect([-width / 2, 0, -sub_h], [width, MSL_length, 0], name="GND")
GND.assign_perfect_E("_perfE")

# define vacuum
cover = chip_body.box([-width / 2, 0, 0], [width, MSL_length, cover_H], name="air_top")
cover.assign_material("vacuum")

# ###########################################################define ports
port1 = chip_body.rect([-width / 2, 0, -sub_h], [width, 0, cover_H + sub_h], name="1")
port1.assign_waveport(Nmodes=1)
if simtype == "DrivenTerminal":
    port1.assign_terminal_auto(GND)


port2 = chip_body.rect([-width / 2, MSL_length, -sub_h], [width, 0, cover_H + sub_h], name="2")
port2.assign_waveport(Nmodes=1)
if simtype == "DrivenTerminal":
    port2.assign_terminal_auto(GND)

############################################################# SOLVE
dm_setup = design.create_dm_setup(freq_ghz=5, name=setup_name)
dm_sweep = dm_setup.insert_sweep(0.1, 10, 1001, name="Sweep", type="Interpolating")
analysis = dm_setup.analyze()

# #############################################################export
path = project.get_path()
file = "test.s2p"
solutions = dm_setup.get_solutions()
solutions.export_network_data(sweep="Sweep", efile=path + file)
print("FILE IS LOCATED IN " + path + file)
