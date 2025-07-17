import math
def compute_stage_mass(radius_tank,fuel_height,oxidizer_height,outer_diameter_body,height_body,n_engines):
    # Constants
    g0 = 9.81                             # gravity (m/s²)
    pressure_bar = 4
    pressure_pa = pressure_bar * 1e5

    # Tank geometry
    #radius_tank = 5.5/2                   # radius of tanks (m)
    #oxidizer_height = 41.16              # m
    #fuel_height = 22.88                  # m

    # Body tube
    outer_diameter_body = 6            # m
    #height_body = 67.13                   # m (entire stage)
    wall_fraction = 0.0008                # wall thickness as % of diameter (1.5%)

    # Material
    materials = {
    "Al-Li 2195": (2800, 500e6),
    "Stainless Steel 304": (7930, 1000e6),
    "Titanium Ti-6Al-4V": (4430, 880e6),
    "Carbon fibre composite": (1600, 600e6)
    }

    # Fuel and oxidizer tank 
    material_name = "Stainless Steel 304"
    density, yield_strength = materials[material_name]
    safety_factor = 2.0
    allowable_stress = yield_strength / safety_factor

    # Tank thickness
    tank_thickness = (pressure_pa * radius_tank) / allowable_stress

    def compute_tank_mass(height):
        cylinder_area = 2 * math.pi * radius_tank * height
        cylinder_volume = cylinder_area * tank_thickness
        dome_area = 4 * math.pi * radius_tank**2
        dome_volume = dome_area * tank_thickness
        total_volume = cylinder_volume + dome_volume
        return total_volume * density

    # Body tube thickness
    wall_thickness_body = outer_diameter_body * wall_fraction
    inner_radius_body = outer_diameter_body / 2 - wall_thickness_body
    outer_radius_body = outer_diameter_body / 2

    body_cylinder_area = 2 * math.pi * outer_radius_body * height_body
    body_volume = body_cylinder_area * wall_thickness_body
    body_mass = body_volume * materials["Stainless Steel 304"][0]

    # Forward fairing
    forward_fairing_length = 4.0  # meters
    cone_area_forward = math.pi * outer_diameter_body / 2 * (
    (outer_diameter_body**2 / 4 + forward_fairing_length**2) ** 0.5)
    fairing_thickness = 0.01
    forward_fairing_mass = cone_area_forward * fairing_thickness * materials["Carbon fibre composite"][0]

    # Aft fairing
    aft_fairing_length = 4.0
    cone_area_aft = math.pi * outer_diameter_body / 2 * (
    (outer_diameter_body**2 / 4 + aft_fairing_length**2) ** 0.5)
    aft_fairing_mass = cone_area_aft * fairing_thickness * materials["Carbon fibre composite"][0]

    # Engine mount structure
    # Approximated as 2% of tank + fairing mass
    oxidizer_mass = compute_tank_mass(oxidizer_height)
    fuel_mass = compute_tank_mass(fuel_height)
    engine_structure_mass = 0.1 * (oxidizer_mass + fuel_mass + body_mass)

    # Total structural mass
    total_structural_mass = (
    oxidizer_mass +
    fuel_mass +
    body_mass +
    forward_fairing_mass +
    aft_fairing_mass +
    engine_structure_mass +
    n_engines*1600
    )

    # Output
    """print(f"Material: {material_name}")
    print(f"Tank Wall Thickness: {tank_thickness*1000:.2f} mm")
    print(f"Body Tube Thickness: {wall_thickness_body*1000:.2f} mm\n")

    print(f"Oxidizer Tank Mass: {oxidizer_mass:.2f} kg")
    print(f"Fuel Tank Mass: {fuel_mass:.2f} kg")
    print(f"Body Tube Mass: {body_mass:.2f} kg")
    print(f"Forward Fairing Mass: {forward_fairing_mass:.2f} kg")
    print(f"Aft Fairing Mass: {aft_fairing_mass:.2f} kg")
    print(f"Engine Mount Structure Mass: {engine_structure_mass:.2f} kg\n")

    print(f"Total Stage 1 Structural Mass: {total_structural_mass:.2f} kg")"""

    return total_structural_mass*1

