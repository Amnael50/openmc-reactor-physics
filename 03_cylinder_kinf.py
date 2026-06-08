import openmc
import matplotlib.pyplot as plt

UO2 = openmc.Material(name="UO2")  
UO2.set_density('g/cm3', 10.5) 
U_mass = 0.04*235+0.96*238
O_mass = 2*16
UO2_mass = U_mass + O_mass
U_ratio = U_mass / UO2_mass
O_ratio = O_mass / UO2_mass
UO2.add_nuclide("U235", U_ratio*0.04, 'wo')
UO2.add_nuclide("U238", U_ratio*0.96, 'wo')  
UO2.add_nuclide("O16", O_ratio, 'wo')
m = openmc.Materials([UO2])

Cladding_material = openmc.Material(name="Cladding")
Cladding_material.set_density('g/cm3', 6.5)
Cladding_material.add_element('Zr', 1.0)
m.append(Cladding_material)

H2O = openmc.Material(name="H2O")
H2O.set_density('g/cm3', 1.0)
H2O.add_nuclide('H1', 2.0, 'ao')
H2O.add_nuclide('O16', 1.0, 'ao')
H2O.add_s_alpha_beta('c_H_in_H2O')
m.append(H2O)

Cladding_geometry = openmc.Cylinder(r=0.475)
Fuel_geometry = openmc.Cylinder(r=0.41)

dimensions = (1.0, 1.26, 1.3, 1.35, 1.4, 1.45, 1.5, 1.75, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0)

resultats = []
settings = openmc.Settings()
settings.batches = 100
settings.inactive = 20
settings.particles = 3000
settings.source = openmc.IndependentSource(space=openmc.stats.Point((0, 0, 0)))

for dimension in dimensions:
    box = openmc.model.RectangularPrism(dimension, dimension, boundary_type='reflective')
    Fuel_cell = openmc.Cell(name="fuel cell", fill=UO2, region=-Fuel_geometry)
    Cladding_cell = openmc.Cell(name="cladding cell", fill=Cladding_material, region=+Fuel_geometry & -Cladding_geometry)
    Water_cell = openmc.Cell(name="water cell", fill=H2O, region=+Cladding_geometry&-box)
    universe = openmc.Universe(cells=[Fuel_cell, Cladding_cell, Water_cell])
    geometry = openmc.Geometry(universe)

    model = openmc.Model(geometry, m, settings)
    sp_path = model.run(output=False)
    print(f"Calcul de la dimension: {dimension} cm")
    with openmc.StatePoint(sp_path) as sp:
        k = sp.keff.n
        ecart = sp.keff.s
        resultats.append([dimension, k, ecart])

print(resultats)
pitchs, ks, sigmas = list(zip(*resultats))
plt.errorbar(pitchs, ks, yerr=sigmas, fmt='o', capsize=3, label='k∞ ± σ (~180 pcm)')
plt.xlabel("Pitch (cm)")
plt.ylabel("k-inf")
plt.title("k-inf vs Pitch")
plt.axhline(1.0, linestyle='--', label='Seuil de criticité')
plt.axvline(1.26, linestyle='-', color='red', label='REP réel')
plt.legend()
plt.grid(True)
plt.savefig('moderation_curve.png', dpi=150)
plt.show()

