from __init__ import *
# -----------------------------------------------
# Class: delta_v
# Purpose: Calculates total delta-v required to reach orbit,
#          accounting for gravitational, drag, and optional boostback losses.
# -----------------------------------------------
class delta_v:
    def __init__(self, thrust_by_weight_initial, Isp, target_orbit, initial_rocket_mass, boostback):
        import path_planning
        self.boostback = boostback
        self.Isp = Isp
        self.initial_rocekt_mass = initial_rocket_mass
        self.rocket_mass = initial_rocket_mass
        self.thrust_by_weight_initial = thrust_by_weight_initial
        self.thrust = thrust_by_weight_initial * initial_rocket_mass * g
        self.target_orbit = target_orbit

        # Compute orbital and initial velocity
        self.orbit_velocity = (g * R_earth**2 / (R_earth + self.target_orbit*1e3))**0.5
        self.initial_velocity = R_earth * np.cos(latitude) * w_earth

        if self.boostback:
            self.calculate_boostback_delta_v()
        #self.calculate_total_delta_v()
        self.path=path_planning.path_planning(thrust_by_weight_initial,target_orbit*1e3,Isp)
        print(f"Calculating delta v requirements")
        self.path.path_planner(n=10,n_epochs=8,dynamic_n=True,check_history=True)
        self.optimal_orbital_parameters=self.path.mins
        self.gravitational_delta_v=np.abs(self.path.gravity_delta_v)
        self.drag_delta_v=self.path.drag_delta_v
        self.total_delta_v=self.path.total_delta_v
        print(f"TOTAL DELTA V REQUIRED: {self.total_delta_v}\n")


    def calculate_boostback_delta_v(self):
        self.boostback_delta_v = [1500.0,2000.0]  # Placeholder boostback cost

    def calculate_total_delta_v(self):
        # Add losses and orbital requirement
        self.total_delta_v = -self.initial_velocity + self.gravitational_delta_v + self.orbit_velocity + self.drag_delta_v

    def display_breakdown(self, as_table=True):
        out_str="\n###################################   DELTA-V BREAKDOWN   ##############################\n\n"
        #("\n###################################   DELTA-V BREAKDOWN   ##############################\n")
        #required_delta_v = round(self.total_delta_v + np.sum(self.boostback_delta_v), 2) if self.boostback else round(self.total_delta_v, 2)
        required_delta_v = self.total_delta_v
        q=0
        for a in self.boostback:
            if a:
                required_delta_v+=self.boostback_delta_v[q]
            q+=1
        delta_v_items = {
            "Required Delta-v (m/s)": required_delta_v,
            "Initial Velocity (m/s)": round(self.initial_velocity, 2),
            "Orbital Velocity (m/s)": round(self.orbit_velocity, 2),
            "Gravitational Loss (m/s)": round(self.gravitational_delta_v, 2),
            "Drag Loss (m/s)": round(self.drag_delta_v, 2),
        }
        q=1
        for a in self.boostback:
            if a:
                delta_v_items[f"Boostback Delta-v for stage {q} (m/s)"] = round(self.boostback_delta_v[q-1], 2)
            q+=1
        if not as_table:
            for key, value in delta_v_items.items():
                out_str+=f"{key} = {value}"
                #print(f"{key} = {value}")
        else:
            import pandas as pd
            df = pd.DataFrame.from_dict(delta_v_items, orient='index', columns=["Delta-v (m/s)"])
            out_str+=df.to_markdown()+"\n"
            #print(df.to_markdown())
        return out_str

class delta_v_obsolete:
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
        self.initial_velocity = R_earth * np.cos(latitude) * w_earth

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
                break

    def calculate_gravitation_delta_v(self):
        self.gravitational_delta_v = self.flight_time * g

    def calculate_drag_delta_v(self):
        self.drag_delta_v = 50.0  # Placeholder for drag model

    def calculate_boostback_delta_v(self):
        self.boostback_delta_v = [1500.0,2000.0]  # Placeholder boostback cost

    def calculate_total_delta_v(self):
        # Add losses and orbital requirement
        self.total_delta_v = -self.initial_velocity + self.gravitational_delta_v + self.orbit_velocity + self.drag_delta_v

    def display_breakdown(self, as_table=True):
        out_str="\n###################################   DELTA-V BREAKDOWN   ##############################\n\n"
        #("\n###################################   DELTA-V BREAKDOWN   ##############################\n")
        #required_delta_v = round(self.total_delta_v + np.sum(self.boostback_delta_v), 2) if self.boostback else round(self.total_delta_v, 2)
        required_delta_v = self.total_delta_v
        q=0
        for a in self.boostback:
            if a:
                required_delta_v+=self.boostback_delta_v[q]
            q+=1
        delta_v_items = {
            "Required Delta-v (m/s)": required_delta_v,
            "Initial Velocity (m/s)": round(self.initial_velocity, 2),
            "Orbital Velocity (m/s)": round(self.orbit_velocity, 2),
            "Gravitational Loss (m/s)": round(self.gravitational_delta_v, 2),
            "Drag Loss (m/s)": round(self.drag_delta_v, 2),
        }
        q=1
        for a in self.boostback:
            if a:
                delta_v_items[f"Boostback Delta-v for stage {q} (m/s)"] = round(self.boostback_delta_v[q-1], 2)
            q+=1
        if not as_table:
            for key, value in delta_v_items.items():
                out_str+=f"{key} = {value}"
                #print(f"{key} = {value}")
        else:
            import pandas as pd
            df = pd.DataFrame.from_dict(delta_v_items, orient='index', columns=["Delta-v (m/s)"])
            out_str+=df.to_markdown()+"\n"
            #print(df.to_markdown())
        return out_str