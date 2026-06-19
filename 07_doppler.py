import openmc
import openmc.deplete
import matplotlib.pyplot as plt
import numpy as np

openmc.config['chain_file'] = '/home/thomas/nuclear_data/chain_endfb80_pwr.xml'

batches = 100
inactive = 20
particles = 1000

pitch_value = 1.26
reseau = 17
guide_positions = [(2, 5), (2, 8), (2, 11), (3,3), (3,13), (5,2), (5,5), (5,8), (5,11), (5,14), (8,2), (8,5), (8,8), (8,11), (8,14), (11,2), (11,5), (11,8), (11,11), (11,14), (13,3), (13,13), (14,5), (14,8), (14,11)]
energy_bins = np.logspace(-3, 7, 500)
density_to_calculate = (0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)
Volume = np.pi * 0.41**2 * 1  # volume de combustible dans une cellule (cm³)
power_per_cell = 170.98  # puissance par cellule (W)
time_steps = [0.5, 1, 1.5, 2, 5, 10, 15, 20, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 200]

temp_initial = 294
temp_steps = [294, 600, 900, 1200, 2500]

use_only_cellule = True
calculate_density = False
iterate_temperature = True

settings = openmc.Settings()
resultats_k = []
results = []

def fuel_UO2_material(enrichment: float, density: float) -> openmc.Material:
    UO2 = openmc.Material(material_id=1, name="UO2")  
    UO2.set_density('g/cm3', density) 
    U_mass = enrichment*235+(1-enrichment)*238
    O_mass = 2*16
    UO2_mass = U_mass + O_mass
    U_ratio = U_mass / UO2_mass
    O_ratio = O_mass / UO2_mass
    UO2.add_nuclide("U235", U_ratio*enrichment, 'wo')
    UO2.add_nuclide("U238", U_ratio*(1-enrichment), 'wo')  
    UO2.add_nuclide("O16", O_ratio, 'wo')
    return UO2

def cladding_Zr_material(density: float) -> openmc.Material:
    Cladding_material = openmc.Material(name="Cladding")
    Cladding_material.set_density('g/cm3', density)
    Cladding_material.add_element('Zr', 1.0)
    return Cladding_material

def water_material() -> openmc.Material:
    H2O = openmc.Material(name="H2O")
    H2O.set_density('g/cm3', 1.0)
    H2O.add_nuclide('H1', 2.0, 'ao')
    H2O.add_nuclide('O16', 1.0, 'ao')
    H2O.add_s_alpha_beta('c_H_in_H2O')
    return H2O

def water_variable_material(density: float) -> openmc.Material:
    H2O = openmc.Material(name="H2O")
    H2O.set_density('g/cm3', density)
    H2O.add_nuclide('H1', 2.0, 'ao')
    H2O.add_nuclide('O16', 1.0, 'ao')
    H2O.add_s_alpha_beta('c_H_in_H2O')
    return H2O

def void_material() -> openmc.Material:
    Void = openmc.Material(name="Void")
    Void.set_density('g/cm3', 0.0)
    return Void

def parameter_study(settings: openmc.Settings) -> None:
    settings.batches = batches
    settings.inactive = inactive
    settings.particles = particles
    settings.source = openmc.IndependentSource(space=openmc.stats.Point((0, 0, 0)))
    settings.temperature = {'method': 'interpolation'}

def fuel_and_cladding_cylinder(r_cladding:float, r_fuel:float) -> tuple[openmc.Cylinder, openmc.Cylinder]:
    Cladding_geometry = openmc.Cylinder(r=r_cladding)
    Fuel_geometry = openmc.Cylinder(r=r_fuel)
    return Fuel_geometry, Cladding_geometry

def fuel_and_cladding_cell(Fuel_geometry:openmc.Cylinder , Cladding_geometry:openmc.Cylinder, fuel_material:openmc.Material, cladding_material:openmc.Material) -> tuple[openmc.Cell, openmc.Cell]:
    Fuel_cell = openmc.Cell(name="fuel cell", fill=fuel_material, region=-Fuel_geometry)
    Cladding_cell = openmc.Cell(name="cladding cell", fill=cladding_material, region=+Fuel_geometry & -Cladding_geometry)
    return Fuel_cell, Cladding_cell

def build_tallies():
    # tally nappe de puissance (ton mesh existant)
    mesh = openmc.RegularMesh()
    mesh.dimension = [reseau, reseau]
    mesh.lower_left = [-reseau * pitch_value / 2] * 2
    mesh.upper_right = [reseau * pitch_value / 2] * 2
    fission_tally = openmc.Tally(name='fission_map')
    fission_tally.filters = [openmc.MeshFilter(mesh)]
    fission_tally.scores = ['fission']

    # tally spectre
    spec_tally = openmc.Tally(name='spectrum')
    spec_tally.filters = [openmc.EnergyFilter(energy_bins)]
    spec_tally.scores = ['flux']

    return openmc.Tallies([fission_tally, spec_tally]) 

def build_model(density: float, temperature: float):
    UO2 = fuel_UO2_material(0.04, 10.5)
    UO2.temperature = temperature
    Cladding_material = cladding_Zr_material(6.5)
    Cladding_material.temperature = temp_initial
    H2O = water_variable_material(density)
    H2O.temperature = temp_initial

    m = openmc.Materials([UO2])
    m.append(Cladding_material)
    m.append(H2O)

    Fuel_geometry, Cladding_geometry = fuel_and_cladding_cylinder(0.475, 0.41)
    Fuel_cell, Cladding_cell = fuel_and_cladding_cell(Fuel_geometry, Cladding_geometry, UO2, Cladding_material)

    Water_cell = openmc.Cell(name="water cell", fill=H2O, region=+Cladding_geometry)
    water_cell_center = openmc.Cell(name="water center", fill=H2O) 
    water_universe = openmc.Universe(cells=[water_cell_center])

    parameter_study(settings)

    if use_only_cellule:
            print("Calcul avec une seule cellule (sans réseau)")
            box = openmc.model.RectangularPrism(pitch_value, pitch_value, boundary_type='reflective')
            Water_cell = openmc.Cell(fill=H2O, region=+Cladding_geometry & -box)
            universe = openmc.Universe(cells=[Fuel_cell, Cladding_cell, Water_cell])
            geometry = openmc.Geometry(universe)
            model = openmc.Model(geometry, m, settings)
            return model

    print("Calcul en reseau")
    pin_universe = openmc.Universe(cells=[Fuel_cell, Cladding_cell, Water_cell])
    lattice_array = np.full((reseau, reseau), pin_universe)

    for (i, j) in guide_positions:
        lattice_array[i, j] = water_universe

    grid = openmc.RectLattice()

    pitch = (pitch_value, pitch_value)
    grid.lower_left = (-reseau * pitch_value/2, -reseau * pitch_value/2)
    grid.pitch = pitch
    grid.universes = lattice_array
    box = openmc.model.RectangularPrism(reseau * pitch_value, reseau * pitch_value, boundary_type='reflective')
    root_cell = openmc.Cell(fill=grid, region=-box)   # remplie par le RÉSEAU
    root_universe = openmc.Universe(cells=[root_cell])
    geometry = openmc.Geometry(root_universe)
    root_universe.plot(width=(reseau * pitch_value, reseau * pitch_value), pixels=(100*reseau, 100*reseau))
    plt.savefig('geom.png')
    tallies = build_tallies()
    model = openmc.Model(geometry, m, settings, tallies=tallies)
    return model

if calculate_density:
    density = density_to_calculate
else:
    density = (1.0,)

for d in density:
    print(f"Calcul pour une densité de {d:.2f} g/cm³")

    if iterate_temperature:
        temp_values = temp_steps
    else:
        temp_values = (temp_initial,)

    for temp in temp_values:
        print(f"Calcul pour une température de {temp} K")

        model_to_calculate = build_model(d, temp)
      
        sp_path = model_to_calculate.run()

        with openmc.StatePoint(sp_path) as sp:
            k = sp.keff.n
            ecart = sp.keff.s
            resultats_k.append([d, k, ecart])
            results.append([temp, k, ecart])
        
if iterate_temperature:
    data = np.array(results)
    T_arr, k_arr, s_arr = data[:,0], data[:,1], data[:,2]

    # Tracé k vs température
    plt.figure(figsize=(8,6))
    plt.errorbar(T_arr, k_arr, yerr=s_arr, fmt='o-', capsize=3)
    plt.xlabel("Température du combustible (K)")
    plt.ylabel("k∞")
    plt.title("Effet Doppler — cellule REP UO₂ 4%\n(seul le combustible chauffe)")
    plt.grid(True, alpha=0.3)
    plt.savefig('doppler.png', dpi=150)

    # Coefficient Doppler en pcm/K, entre 294 K et 900 K (plage opérationnelle)
    rho = (k_arr - 1) / k_arr * 1e5    # pcm
    # coefficient = Δρ / ΔT
    i1 = np.where(T_arr == 294)[0][0]
    i2 = np.where(T_arr == 900)[0][0]
    alpha_doppler = (rho[i2] - rho[i1]) / (900 - 294)
    print(f"Coefficient Doppler (294→900 K) : {alpha_doppler:.2f} pcm/K")

if calculate_density:
    data = np.array(resultats_k)
    d_arr, k_arr, s_arr = data[:,0], data[:,1], data[:,2]
    rho = (k_arr - 1) / k_arr * 1e5   # réactivité en pcm
    delta_rho = rho[d_arr==1.0][0] - rho[d_arr==0.01][0]
    print(f"Δρ (pleine eau → quasi-vide) = {delta_rho:.0f} pcm")

    plt.figure(figsize=(8,6))
    plt.errorbar(d_arr, k_arr, yerr=s_arr, fmt='o-', capsize=3)
    plt.xlabel("Densité de l'eau (g/cm³)")
    plt.ylabel("k∞")
    plt.title("Coefficient de vide — cellule REP UO₂ 4%")
    plt.axhline(1.0, ls='--', color='grey', label='Seuil critique')
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.savefig('void_coefficient.png', dpi=150)