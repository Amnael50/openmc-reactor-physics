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
Fuel_cell = openmc.Cell(name="fuel cell", fill=UO2, region=-Fuel_geometry)
Cladding_cell = openmc.Cell(name="cladding cell", fill=Cladding_material, region=+Fuel_geometry & -Cladding_geometry)
Water_cell = openmc.Cell(name="water cell", fill=H2O, region=+Cladding_geometry)
water_cell_center = openmc.Cell(name="water center", fill=H2O) 
water_universe = openmc.Universe(cells=[water_cell_center])
pin_universe = openmc.Universe(cells=[Fuel_cell, Cladding_cell, Water_cell])
grid = openmc.RectLattice()
pitch = (1.26, 1.26)
grid.lower_left = (-3 * 1.26/2, -3 * 1.26/2)
grid.pitch = pitch
grid.universes = [[pin_universe, pin_universe, pin_universe],
                  [pin_universe, water_universe, pin_universe],
                  [pin_universe, pin_universe, pin_universe]]
box = openmc.model.RectangularPrism(3.78, 3.78, boundary_type='reflective')
root_cell = openmc.Cell(fill=grid, region=-box)   # remplie par le RÉSEAU
root_universe = openmc.Universe(cells=[root_cell])
geometry = openmc.Geometry(root_universe)
root_universe.plot(width=(3.78, 3.78), pixels=(300, 300))
plt.savefig('geom.png')
settings = openmc.Settings()
settings.batches = 100
settings.inactive = 20
settings.particles = 5000
settings.source = openmc.IndependentSource(space=openmc.stats.Point((0, 0, 0)))

model = openmc.Model(geometry, m, settings)
model.run()