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

calcul_energy_spectrum = False
use_only_cellule = True
calculate_density = False
depletion = True

settings = openmc.Settings()
resultats_k = []

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

def water_variable_material(density: float) -> openmc.Material:
    H2O = openmc.Material(name="H2O")
    H2O.set_density('g/cm3', density)
    H2O.add_nuclide('H1', 2.0, 'ao')
    H2O.add_nuclide('O16', 1.0, 'ao')
    H2O.add_s_alpha_beta('c_H_in_H2O')
    return H2O

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

def build_model(density: float):
    UO2 = fuel_UO2_material(0.04, 10.5)
    if depletion:
        UO2.depletable = True
        UO2.volume = Volume
    Cladding_material = cladding_Zr_material(6.5)
    H2O = water_variable_material(density)

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

    model_to_calculate = build_model(d)

    if depletion:
        op = openmc.deplete.CoupledOperator(model_to_calculate)
        integrator = openmc.deplete.PredictorIntegrator(
            op, time_steps, power_per_cell, timestep_units='d')
        integrator.integrate()
        break   # ou structure le code pour ne pas tomber dans le bloc statique
        
    sp_path = model_to_calculate.run()

    with openmc.StatePoint(sp_path) as sp:
        if calcul_energy_spectrum:
            tal = sp.get_tally(name='fission_map')
            fission_raw = tal.get_values(scores=['fission']).reshape(17, 17)
            ecart_type = tal.get_values(scores=['fission'], value='std_dev').reshape(17, 17)
            t = sp.get_tally(name='spectrum')
            flux = t.get_values(scores=['flux']).flatten()
        k = sp.keff.n
        ecart = sp.keff.s
        resultats_k.append([d, k, ecart])
        
        if calcul_energy_spectrum:
            E = energy_bins                      # 501 bornes pour 500 bins
            E_mid = np.sqrt(E[:-1] * E[1:])      # milieu géométrique de chaque bin (log)
            dlnE = np.diff(np.log(E))            # largeur de léthargie de chaque bin
            flux_per_lethargy = flux / dlnE      # LA correction clé

            plt.figure(figsize=(9, 6))
            plt.step(E_mid, flux_per_lethargy, where='mid')
            plt.xscale('log')
            plt.xlabel('Énergie (eV)')
            plt.ylabel('Flux par unité de léthargie (u.a.)')
            title = f"Spectre neutronique — Assemblage REP 17×17\nUO₂ 4%, densite : {d:.2f}, OpenMC ENDF/B-VIII.0"
            plt.title(title)
            plt.grid(True, which='both', alpha=0.3)
            plt.savefig('spectrum.png', dpi=150)

            # Normalisation : facteur de forme (1.0 = puissance moyenne)
            fission_map = fission_raw / fission_raw[fission_raw > 0].mean()
            print (f"Fxy crayon maximum: {fission_map.max()}")
            print (f"Position crayon: {np.unravel_index(fission_map.argmax(), fission_map.shape)}")
            print (f"Écart type relatif: {ecart_type[fission_raw > 0].mean() / fission_raw[fission_raw > 0].mean():.2%}")
            # Heatmap
            fig, ax = plt.subplots(figsize=(8, 8))
            im = ax.imshow(fission_map, cmap='hot', origin='upper',
                    vmin=0, vmax=fission_map.max())
            plt.colorbar(im, ax=ax, label='Facteur de puissance relatif')
            ax.set_title('Nappe de puissance — Assemblage REP 17×17\n'
                    'UO₂ 4%, pitch 1.26 cm, OpenMC ENDF/B-VIII.0')
            ax.set_xlabel('Position j')
            ax.set_ylabel('Position i')
            title_power_map = f"power_map_{d:.2f}.png"
            plt.savefig(title_power_map, dpi=150, bbox_inches='tight')

if depletion:
    results = openmc.deplete.Results("depletion_results.h5")
    _, u235  = results.get_atoms("1", "U235")
    _, pu239 = results.get_atoms("1", "Pu239")
    time, keff = results.get_keff()
    # time en secondes, keff[:,0] = valeurs, keff[:,1] = incertitudes
    
    # Convertir le temps en burnup (GWj/t)
    masse_U_tonnes = Volume * 10.5 * 0.881 * 1e-6   # masse U en tonnes
    burnup = power_per_cell * (time / 86400) / 1e9 / masse_U_tonnes  # GWj/t
    # (time/86400 : secondes → jours ; power en W → GW : /1e9)
    
    plt.figure(figsize=(9,6))
    plt.errorbar(burnup, keff[:,0], yerr=keff[:,1], fmt='o-', capsize=3)
    plt.xlabel("Burnup (GWj/t)")
    plt.ylabel("k∞")
    plt.title("Évolution de k∞ avec le burnup — cellule REP UO₂ 4%")
    plt.axhline(1.0, ls='--', color='grey')
    plt.grid(True, alpha=0.3)
    plt.savefig('burnup_keff.png', dpi=150)

    # --- Courbe d'évolution isotopique ---
    plt.figure(figsize=(9, 6))
    plt.plot(burnup, u235,  'o-', label='U-235')
    plt.plot(burnup, pu239, 's-', label='Pu-239')
    plt.xlabel("Burnup (GWj/t)")
    plt.ylabel("Nombre d'atomes")
    plt.title("Évolution isotopique avec le burnup — cellule REP UO₂ 4%")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.savefig('burnup_isotopes.png', dpi=150)

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