import openmc

fuel_material = openmc.Material(name="U235-20%")  
fuel_material.set_density('g/cm3', 19.0) 
fuel_material.add_nuclide("U235", 0.20, 'wo')
fuel_material.add_nuclide("U238", 0.80, 'wo')  
m = openmc.Materials([fuel_material])

fuel_geometry = openmc.Sphere(r=16, boundary_type='vacuum')

fuel_cell = openmc.Cell(name="fuel cell", fill=fuel_material, region=-fuel_geometry)
universe = openmc.Universe(cells=[fuel_cell])
geometry = openmc.Geometry(universe)

settings = openmc.Settings()
settings.batches = 100
settings.inactive = 20
settings.particles = 5000
settings.source = openmc.Source(space=openmc.stats.Point((0, 0, 0)))

model = openmc.Model(geometry, m, settings)
model.run()
