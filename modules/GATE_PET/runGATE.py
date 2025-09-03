#GATE Module for PET SCANS

import opengate as gate
'''
#The Simulation object
sim = gate.Simulation()

sim.g4_verbose = False
sim.visu = False
sim.random_seed = 'auto'
sim.number_of_threads = 1

gate.help_on_user_info(gate.Simulation)
gate.help_on_user_info(gate.Simulation.number_of_threads)
'''

''' 
#units
cm = gate.g4_units.cm
eV = gate.g4_units.eV
MeV = gate.g4_units.MeV
x = 32 * cm
energy = 150 * MeV
print(f'The energy is {energy/eV} eV')
'''

#Components of a simulation

if __name__ == "__main__":
    sim = gate.Simulation()
    sim.output_dir = "/Users/peteryazdi/Desktop/BC_Cancer/TDT/modules/GATE_PET/output"

    #Phantom:
    wb = sim.add_volume("Box", name="waterbox")
    # configure the volume ...
    cm = gate.g4_units.cm
    wb.size = [10 * cm, 5 * cm, 10 * cm]
    wb.material = "G4_WATER" # so has similar properties to human tissue
    # ...

    stats = sim.add_actor("SimulationStatisticsActor", "Stats")

    #Loading physics: QGSP_BERT_EMV (hadronic + EM interactions).
    #Spawning primaries: 1000 240 MeV protons from a point at the origin with isotropic directions.
    #Tracking each particle through the geometry and applying the physics.
    #Actors: you only added SimulationStatisticsActor, so it prints run stats; 
    # nothing is saved unless you give actors filenames (dose, phase space, etc.).
    source = sim.add_source("GenericSource", name="Default")
    source.particle = "proton"
    source.energy.type = "mono"
    source.energy.mono = 240 * gate.g4_units.MeV
    source.n = 1000  # emit exactly 1000 particles ( keep small if vis on) 
    source.position.type = "point"
    source.position.translation = [0, 0, 0]
    source.direction.type = "iso"
    # ...
    
    #dose actor
    dose = sim.add_actor("DoseActor", "DoseInWater")
    dose.attached_to = wb
    dose.dose.active = True
    dose.dose.output_filename = "dose_in_water.mhd"
    
    stats.output_filename = "stats.txt"
    
    #visualtion:
    sim.visu = True
    sim.visu_type = "vrml"
    sim.visu_filename = "scene.wrl"
    
    #JSON save
    sim.store_json_archive = True
    sim.json_archive_filename = "sim_config.json"  # saved under sim.output_dir

    sim.run()
    
    

def runGATE():
    return