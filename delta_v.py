from __init__ import *

# -----------------------------------------------
# Class: delta_v
# Purpose: Calculates total delta-v required to reach orbit,
#          accounting for gravitational, drag, and optional boostback losses.
# -----------------------------------------------
class delta_v:
    def __init__(self, thrust_by_weight_initial, Isp, target_orbit, initial_rocket_mass, boostback):
        self.boostback = boostback
        self.Isp = Isp
        self.initial_rocekt_mass = initial_rocket_mass
        self.rocket_mass = initial_rocket_mass
        self.thrust_by_weight_initial = thrust_by_weight_initial
        self.thrust = thrust_by_weight_initial * initial_rocket_mass * g
        self.target_orbit = target_orbit

        # Compute orbital and initial velocity
        self.orbit_velocity = (g * R_earth**2 / (R_earth + self.target_orbit))**0.5
        self.initial_velocity = R_earth * np.cos(latitude / 180 * np.pi) * w_earth

        # Initialize losses
        self.steering_velocity_loss = 0
        self.velocity = 0
        self.flight_time = 0

        # Run all computations
        self.calculate_flight_time()
        self.calculate_drag_delta_v()
        self.calculate_gravitation_delta_v()
        if self.boostback:
            self.calculate_boostback_delta_v()
        self.calculate_total_delta_v()

    def calculate_flight_time(self):
        # Simulate ascent until orbital velocity is achieved
        while self.velocity < self.orbit_velocity:
            self.velocity += (self.thrust / self.rocket_mass - g) * dt
            m_dot = self.thrust / self.Isp / g
            self.rocket_mass -= m_dot * dt
            self.flight_time += dt
            if self.rocket_mass < 0:
                print("Insufficient Fuel")
                break

    def calculate_gravitation_delta_v(self):
        self.gravitational_delta_v = self.flight_time * g

    def calculate_drag_delta_v(self):
        self.drag_delta_v = 50.0  # Placeholder for drag model

    def calculate_boostback_delta_v(self):
        self.boostback_delta_v = 1200.0  # Placeholder boostback cost

    def calculate_total_delta_v(self):
        # Add losses and orbital requirement
        self.total_delta_v = -self.initial_velocity + self.gravitational_delta_v + self.orbit_velocity + self.drag_delta_v

    def display_breakdown(self, as_table=True):
        print("\n###################################   DELTA-V BREAKDOWN   ##############################\n")
        required_delta_v = round(self.total_delta_v + self.boostback_delta_v, 2) if self.boostback else round(self.total_delta_v, 2)

        delta_v_items = {
            "Required Delta-v (m/s)": required_delta_v,
            "Initial Velocity (m/s)": round(self.initial_velocity, 2),
            "Orbital Velocity (m/s)": round(self.orbit_velocity, 2),
            "Gravitational Loss (m/s)": round(self.gravitational_delta_v, 2),
            "Drag Loss (m/s)": round(self.drag_delta_v, 2),
        }

        if self.boostback:
            delta_v_items["Boostback Delta-v (m/s)"] = round(self.boostback_delta_v, 2)

        if not as_table:
            for key, value in delta_v_items.items():
                print(f"{key} = {value}")
        else:
            import pandas as pd
            df = pd.DataFrame.from_dict(delta_v_items, orient='index', columns=["Delta-v (m/s)"])
            print(df.to_markdown())
        return delta_v_items

# -----------------------------------------------
# Class: delta_v_2
# Purpose: Computes delta-v achieved with a fixed amount of fuel
# -----------------------------------------------
class delta_v_2:
    def __init__(self, thrust_by_weight_initial, Isp, initial_rocket_mass, total_fuel_available):
        self.Isp = Isp
        self.initial_rocekt_mass = initial_rocket_mass
        self.rocket_mass = initial_rocket_mass
        self.thrust_by_weight_initial = thrust_by_weight_initial
        self.thrust = thrust_by_weight_initial * initial_rocket_mass * g
        self.initial_fuel_weight = total_fuel_available

        self.calculate_flight_time()
        self.calculate_drag_delta_v()
        self.calculate_gravitation_delta_v()
        self.calculate_total_delta_v()

    def calculate_flight_time(self):
        current_fuel = self.initial_fuel_weight
        self.flight_time = 0
        self.current_velocity = 0
        while current_fuel > 0:
            self.current_velocity += (self.thrust / self.rocket_mass - g) * dt
            m_dot = self.thrust / self.Isp / g
            self.rocket_mass -= m_dot * dt
            current_fuel -= m_dot * dt
            self.flight_time += dt

    def calculate_gravitation_delta_v(self):
        self.gravitational_delta_v = self.flight_time * g

    def calculate_drag_delta_v(self):
        self.drag_delta_v = 50  # Placeholder

    def calculate_total_delta_v(self):
        self.total_delta_v = self.gravitational_delta_v + self.current_velocity + self.drag_delta_v

# -----------------------------------------------
# Class: delta_v_3
# Purpose: Delta-v required to orbit using a physics-based model
#          with thrust limiting, drag, gravity, and steering losses
# -----------------------------------------------
class delta_v_3:
    def __init__(self, thrust_by_weight_initial, Isp, target_orbit, initial_rocket_mass,
                 thrust_by_weight_limit=1000, steering_angle_offset=0, initial_orientation=0):
        self.Isp = Isp
        self.initial_rocket_mass = initial_rocket_mass
        self.rocket_mass = initial_rocket_mass
        self.thrust_by_weight_initial = thrust_by_weight_initial
        self.thrust_by_weight_limit = thrust_by_weight_limit
        self.thrust = thrust_by_weight_initial * initial_rocket_mass * g
        self.target_orbit = target_orbit
        self.steering_angle_offset = steering_angle_offset
        self.orbit_velocity = (g * R_earth**2 / (R_earth + self.target_orbit))**0.5

        self.steering_velocity_loss = 0
        self.gravitational_delta_v = 0
        self.drag_delta_v = 0
        self.orientation = initial_orientation
        self.velocity = 0
        self.total_delta_v = 0
        self.flight_time = 0

        self.calculate_total_delta_v()

    def calculate_drag(self):
        # Basic drag model
        return self.velocity**2 * 0.5 * np.pi * 3**2 * 0.00001

    def calculate_total_delta_v(self):
        while self.velocity < self.orbit_velocity:
            # Enforce thrust limit
            if self.thrust / self.rocket_mass / g > self.thrust_by_weight_limit:
                self.thrust = self.thrust_by_weight_limit * self.rocket_mass * g

            # Update physics terms
            self.velocity += ((self.thrust - self.calculate_drag()) / self.rocket_mass - g * np.cos(self.orientation)) * dt
            m_dot = self.thrust / self.Isp / g

            # Accumulate losses
            self.gravitational_delta_v += g * dt
            self.drag_delta_v += self.calculate_drag() * dt / self.rocket_mass
            self.steering_velocity_loss += self.thrust * self.velocity * (1 - np.cos(self.orientation)) * dt

            # Update mass
            self.rocket_mass -= m_dot * dt
            self.flight_time += dt
            if self.rocket_mass < 0:
                print("Insufficient Fuel")
                break

        self.total_delta_v = self.velocity + self.gravitational_delta_v + self.drag_delta_v + self.steering_velocity_loss

# -----------------------------------------------
# Class: delta_v_4
# Purpose: Variant of delta_v_3 with fuel limit (instead of velocity goal)
# -----------------------------------------------
class delta_v_4:
    def __init__(self, thrust_by_weight_initial, Isp, initial_rocket_mass, fuel_available,
                 thrust_by_weight_limit=1000, steering_angle_offset=0, initial_orientation=0):
        self.Isp = Isp
        self.initial_rocket_mass = initial_rocket_mass
        self.rocket_mass = initial_rocket_mass
        self.thrust_by_weight_initial = thrust_by_weight_initial
        self.thrust_by_weight_limit = thrust_by_weight_limit
        self.thrust = thrust_by_weight_initial * initial_rocket_mass * g
        self.fuel_available = fuel_available
        self.steering_angle_offset = steering_angle_offset

        self.steering_velocity_loss = 0
        self.gravitational_delta_v = 0
        self.drag_delta_v = 0
        self.orientation = initial_orientation
        self.velocity = 0
        self.total_delta_v = 0
        self.flight_time = 0

        self.calculate_total_delta_v()

    def calculate_drag(self):
        return self.velocity**2 * 0.5 * np.pi * 3**2 * 0.00001

    def calculate_total_delta_v(self):
        while self.fuel_available > 0:
            if self.thrust / self.rocket_mass / g > self.thrust_by_weight_limit:
                self.thrust = self.thrust_by_weight_limit * self.rocket_mass * g

            self.velocity += ((self.thrust - self.calculate_drag()) / self.rocket_mass - g * np.cos(self.orientation)) * dt
            m_dot = self.thrust / self.Isp / g

            self.gravitational_delta_v += g * dt
            self.drag_delta_v += self.calculate_drag() * dt / self.rocket_mass
            self.steering_velocity_loss += self.thrust * self.velocity * (1 - np.cos(self.orientation)) * dt

            self.rocket_mass -= m_dot * dt
            self.fuel_available -= m_dot * dt
            self.flight_time += dt
            if self.rocket_mass < 0:
                print("Insufficient Fuel")
                break

        self.total_delta_v = self.velocity + self.gravitational_delta_v + self.drag_delta_v + self.steering_velocity_loss
