import openmc
import numpy as np
import matplotlib.pyplot as plt

batches = 100
inactive = 20
particles = 1000

r_fuel, r_clad, H = 0.41, 0.475, 365
pitch = 1.26
reseau = 17
insertion = [1.0, 0.90, 0.75, 0.60, 0.50, 0.35, 0.25, 0.10, 0.00]
guide_positions = [(2, 5), (2, 8), (2, 11), (3,3), (3,13), (5,2), (5,5), (5,8), (5,11), (5,14), (8,2), (8,5), (8,8), (8,11), (8,14), (11,2), (11,5), (11,8), (11,11), (11,14), (13,3), (13,13), (14,5), (14,8), (14,11)]

resultats = []

z0 = openmc.ZPlane(z0=0.0, boundary_type='vacuum')
z1 = openmc.ZPlane(z0=H,   boundary_type='vacuum')

def create_fuel_UO2_material(enrichment: float, density: float) -> openmc.Material:
    UO2 = openmc.Material(name="UO2")  
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

def create_cladding_Zr_material(density: float) -> openmc.Material:
    Cladding_material = openmc.Material(name="Cladding")
    Cladding_material.set_density('g/cm3', density)
    Cladding_material.add_element('Zr', 1.0)
    return Cladding_material

def create_control_rod_material():
    b4c = openmc.Material(name="B4C")
    b4c.set_density('g/cm3', 2.52)
    b4c.add_nuclide('B10', 4 * 0.199, 'ao')   # B naturel : ~19.9% B10
    b4c.add_nuclide('B11', 4 * 0.801, 'ao')
    b4c.add_element('C', 1.0, 'ao')
    return b4c

def create_water_material() -> openmc.Material:
    H2O = openmc.Material(name="H2O")
    H2O.set_density('g/cm3', 1.0)
    H2O.add_nuclide('H1', 2.0, 'ao')
    H2O.add_nuclide('O16', 1.0, 'ao')
    H2O.add_s_alpha_beta('c_H_in_H2O')
    return H2O

def create_water_variable_material(density: float) -> openmc.Material:
    H2O = openmc.Material(name="H2O")
    H2O.set_density('g/cm3', density)
    H2O.add_nuclide('H1', 2.0, 'ao')
    H2O.add_nuclide('O16', 1.0, 'ao')
    H2O.add_s_alpha_beta('c_H_in_H2O')
    return H2O

def create_void_material() -> openmc.Material:
    Void = openmc.Material(name="Void")
    Void.set_density('g/cm3', 0.0)
    return Void

def parameter_study(settings: openmc.Settings) -> None:
    settings.batches = batches
    settings.inactive = inactive
    settings.particles = particles
    settings.source = openmc.IndependentSource(space=openmc.stats.Point((0, 0, 0)))

def make_fuel_rod_universe(r_cladding:float, r_fuel:float, fuel_material:openmc.Material, cladding_material:openmc.Material, water_material:openmc.Material) -> openmc.Universe:
    Cladding_geometry = openmc.ZCylinder(r=r_cladding)
    Fuel_geometry = openmc.ZCylinder(r=r_fuel)
    Fuel_cell = openmc.Cell(name="fuel cell", fill=fuel_material, region=-Fuel_geometry & +z0 & -z1)
    Cladding_cell = openmc.Cell(name="cladding cell", fill=cladding_material, region=+Fuel_geometry & -Cladding_geometry & +z0 & -z1)
    water_cell  = openmc.Cell(name="water_cell", fill=water_material, region=+Cladding_geometry & +z0 & -z1)
    return openmc.Universe(cells=[Fuel_cell, Cladding_cell, water_cell])

def make_guide_rod_universe(fill_rod_mat: openmc.Material, fill_water_mat: openmc.Material, z_insert: float) -> openmc.Universe:
    z_mid = openmc.ZPlane(z0=z_insert)
    water_part = openmc.Cell(fill=fill_water_mat, region=+z0 & -z_mid)   # bas : eau
    b4c_part   = openmc.Cell(fill=fill_rod_mat,   region=+z_mid & -z1)   # haut : B4C
    return openmc.Universe(cells=[water_part, b4c_part])

def build(rod_inserted: bool, z_insert:float):
    fuel_material = create_fuel_UO2_material(0.04, 10.5)
    cladding_material = create_cladding_Zr_material(6.55)
    water_material = create_water_variable_material(1.0)
    control_rod_material = create_control_rod_material()
    materials = openmc.Materials([fuel_material, cladding_material, water_material, control_rod_material])
    pin = make_fuel_rod_universe(r_clad, r_fuel, fuel_material, cladding_material, water_material)
    guide = make_guide_rod_universe(control_rod_material, water_material, z_insert)

    arr = np.full((reseau, reseau), pin)
    for (i, j) in guide_positions:
        arr[i, j] = guide

    lat = openmc.RectLattice()
    lat.pitch = (pitch, pitch)
    lat.lower_left = (-reseau*pitch/2, -reseau*pitch/2)
    lat.universes = arr

    box = openmc.model.RectangularPrism(reseau*pitch, reseau*pitch, boundary_type='reflective')
    root = openmc.Cell(fill=lat, region=-box & +z0 & -z1)
    geom = openmc.Geometry([root])

    s = openmc.Settings()
    s.batches, s.inactive, s.particles = 100, 20, 2000
    s.source = openmc.IndependentSource(space=openmc.stats.Point((0,0,H/2)))

    return openmc.Model(geom, materials, s)

# --- Étape 2 : grappe tout-ou-rien ---
def run_k(rod_inserted, z_insert):
    model = build(rod_inserted, (1 - z_insert) * H)
    sp_path = model.run(output=False)
    with openmc.StatePoint(sp_path) as sp:
        k = sp.keff.n
        ecart = sp.keff.s
        resultats.append([z_insert, k, ecart])
        return k, ecart

for z_insert in insertion:
    run_k(rod_inserted=True, z_insert=z_insert)

data = np.array(resultats)
insertion, k_arr, s_arr = data[:,0], data[:,1], data[:,2]
k_extract = k_arr.max()
k_insert = k_arr.min()

for z, k, s in resultats:
    print(f"insertion {z:.2f} : k = {k:.5f} ± {s:.5f}")

curve = []
for k in k_arr:
        curve.append((k_extract-k)*100000)

m = build(True, (1-0.5)*H)
m.geometry.root_universe.plot(basis='xz', width=(reseau*pitch, H), pixels=(400,600))
plt.savefig('rod_xz.png', dpi=150)

plt.figure(figsize=(8,6))
plt.errorbar(insertion, curve, yerr=s_arr*1e5, fmt='o-', capsize=3)
plt.xlabel("Fraction d'insertion de la grappe")
plt.ylabel("Anti-réactivité intégrée (pcm)")
plt.title("Courbe d'insertion de grappe — assemblage REP 17×17")
plt.grid(True, alpha=0.3)
plt.savefig('control_rod_curve.png', dpi=150)

