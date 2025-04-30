from __init__ import *

class rocket:

    def __init__(self,payload):
        self.number_of_stages=0
        self.total_mass=payload
        self.total_fuel_mass=0
        self.position=np.array([0.0,0.0,0.0])
        self.velocity=np.array([0.0,0.0,0.0])
        self.acceleration=np.array([0.0,0.0,0.0])
        self.stages={}
        self.stage_masses=[]
        self.step_fuel_masses=[]
        self.step_masses=[]
        self.thrust_by_weight=[]
        self.Isps=[]
        self.payload=payload

    def add_rocket_data(self,data):
        self.data=data

    def add_stage(self,step_propellant_mass,stage_mass,step_mass,Isp):
        #adds a stage at to the rocket
        self.number_of_stages+=1
        self.stage_masses.append(stage_mass)
        self.step_fuel_masses.append(step_propellant_mass)
        self.step_masses.append(step_mass)
        self.total_fuel_mass+=step_propellant_mass
        self.Isps.append(Isp)
        self.total_mass+=step_mass
        '''
        self.total_fuel_mass+=data['fuel_mass']
        self.total_mass+=data['mass']
        self.stages.update({stage:data})
        '''
        return
        
    def display_states(self,position=True,velocity=True,acceleration=True,masses=True,number_of_stages=True):
        if position:
            print("Position:",self.position)
        if velocity:
            print("Velocity:",self.velocity)
        if acceleration:
            print("Acceleration:",self.acceleration)
        if masses:
            print("Total Mass:",self.total_mass)
            print("Total Fuel Mass:",self.total_fuel_mass)
        if number_of_stages:
            print("Number of stages:",self.number_of_stages)
        

    def drag_force(self, Cd=0.5, A=113):
        rho = 1.225 * np.exp(-self.position[0]/8500)  # kg/m³
        F_drag = 0.5 * Cd * rho * self.velocity[0]**2 * A
        return F_drag
    
    def simulate_burn(self):
        time=0
        acc=np.array([])
        vel=np.array([])
        pos=np.array([])
        print("SIMULATING ROCKET LIFTOFF")
        for i in range(0,self.number_of_stages):
            thrust=self.thrust_by_weight[i]*self.total_mass*g
            while self.stage_fuel_masses[i]>0:
                thrust_by_weight=thrust/self.total_mass/g
                #if thrust_by_weight>4:
                #    thrust=4*self.total_mass*g
                m_dot=thrust/self.Isps[i]/g
                acceleration=(thrust-self.drag_force())/self.total_mass-g
                self.acceleration=np.array([acceleration,0.0,0.0])
                self.velocity+=self.acceleration*dt
                self.position+=self.velocity*dt+0.5*self.acceleration*dt**2
                '''
                acc=np.append(acc,np.array([thrust/self.stage_masses[i],0,0]))
                vel=np.append(acc,np.array([thrust/self.stage_masses[i],0,0]))
                '''
                self.step_fuel_masses[i]-=m_dot*dt
                self.stage_masses[i]-=m_dot*dt
                self.total_mass-=m_dot*dt
                self.total_fuel_mass-=m_dot*dt
                time+=dt
                #print(time,self.position[0],self.stage_masses[i],thrust/self.stage_masses[i]/g)
            self.total_mass-=self.step_masses[i]
            print("STAGE ",i+1," BURNOUT AT ",self.position[0]/1000," KM ALTITUDE", " AT ",time," SECONDS")

        
    def simulate_burn_2(self):
        time=0
        acc=np.array([])
        vel=np.array([])
        pos=np.array([])
        print("SIMULATING ROCKET LIFTOFF")
        for i in range(0,self.number_of_stages):
            thrust=self.thrust_by_weight[i]*self.total_mass*g
            while self.step_fuel_masses[i]>0:
                #thrust_by_weight=thrust/self.total_mass/g
                #if thrust_by_weight>4:
                    #thrust=4*self.total_mass*g
                m_dot=thrust/self.Isps[i]/g
                acceleration=(thrust-self.drag_force())/self.total_mass-g
                self.acceleration=np.array([acceleration,0.0,0.0])
                self.velocity+=self.acceleration*dt
                self.position+=self.velocity*dt#+0.5*self.acceleration*dt**2
                '''
                acc=np.append(acc,np.array([thrust/self.stage_masses[i],0,0]))
                vel=np.append(acc,np.array([thrust/self.stage_masses[i],0,0]))
                '''
                self.step_fuel_masses[i]-=m_dot*dt
                self.stage_masses[i]-=m_dot*dt
                self.total_mass-=m_dot*dt
                self.total_fuel_mass-=m_dot*dt
                self.step_masses-=m_dot*dt
                time+=dt
                #if time%10<0.01:
                    #print(time,self.position[0],self.stage_masses[i],thrust/self.total_mass/g)
            self.total_mass-=self.step_masses[i]
            print("STAGE ",i+1," BURNOUT AT ",round(self.position[0]/1000,0)," KM ALTITUDE", " AT ",round(time,2)," SECONDS")
            #self.display_states()

    def size_rocket(self):
        self.rocket_height=0
        self.step_height=np.zeros(self.number_of_stages)
        self.step_engine_number=np.zeros(self.number_of_stages)
        self.step_oxidiser_tank_height=np.zeros(self.number_of_stages)
        self.intertank_height=np.ones(self.number_of_stages)+1.5
        self.step_fuel_tank_height=np.zeros(self.number_of_stages)
        for i in range(0,self.number_of_stages):
            oxidiser_mass=self.step_fuel_masses[i]*self.data["o/f"]/(self.data["o/f"]+1)
            fuel_mass=self.step_fuel_masses[i]-oxidiser_mass
            self.step_oxidiser_tank_height[i]=(oxidiser_mass/self.data["oxidiser density"]-4/3*np.pi*(self.data["tank diameters"]/2)**3)/np.pi/(self.data["tank diameters"]/2)**2
            self.step_fuel_tank_height[i]=(fuel_mass/self.data["fuel density"]-4/3*np.pi*(self.data["tank diameters"]/2)**3)/np.pi/(self.data["tank diameters"]/2)**2
            self.step_height[i]=self.step_fuel_tank_height[i]+self.step_oxidiser_tank_height[i]+self.data["engine height"]
            self.step_height[i]+=self.intertank_height[i]
            self.rocket_height+=self.step_height[i]
            self.step_engine_number[i]=self.stage_masses[i]*self.data["thrust by weight"][i]*g/self.data["engine thrust"][i]
            print(self.step_height[i])
        self.rocket_height+=5 #interstage
        print("Rocket Height:",self.rocket_height)
        print("Number Of Engine:",self.step_engine_number)


    def tabulate_rocket_stages(self):
        #print()
        print()
        print("###################################   ROCKET DIMENSIONS   ##############################")
        print()
        #print()
        tank_radius = self.data["tank diameters"] / 2
        tank_area = pi * tank_radius**2
        hemisphere_volume = (4/3) * pi * tank_radius**3

        results = {
            "Stage": [],
            "Fuel Mass (kg)": [],
            "Oxidiser Mass (kg)": [],
            "Fuel Tank Height (m)": [],
            "Oxidiser Tank Height (m)": [],
            "Intertank Height (m)": [],
            "Engine Count": [],
            "Total Stage Height (m)": [],
        }

        self.rocket_height = 0
        self.step_height = np.zeros(self.number_of_stages)
        self.step_engine_number = np.zeros(self.number_of_stages)
        self.step_oxidiser_tank_height = np.zeros(self.number_of_stages)
        self.step_fuel_tank_height = np.zeros(self.number_of_stages)
        self.intertank_height = np.ones(self.number_of_stages) * 0

        for i in range(self.number_of_stages):
            # Mass breakdown

            oxf = self.data["o/f"]
            ox_density = self.data["oxidiser density"]
            fuel_density = self.data["fuel density"]

            fuel_mass = self.step_fuel_masses[i] / (oxf + 1)
            oxidiser_mass = self.step_fuel_masses[i] - fuel_mass

            # Tank heights
            fuel_tank_height = (fuel_mass / fuel_density - hemisphere_volume) / tank_area + tank_radius*2
            oxidiser_tank_height = (oxidiser_mass / ox_density - hemisphere_volume) / tank_area + tank_radius*2

            # Total stage height
            total_height = fuel_tank_height + oxidiser_tank_height + self.data["engine height"] + self.intertank_height[i]

            # Engine count
            required_thrust = self.stage_masses[i] * self.data["thrust by weight"][i] * g
            engine_count = ceil(required_thrust / self.data["engine thrust"][i])

        # Store
            results["Stage"].append(f"Stage {i+1}")
            results["Fuel Mass (kg)"].append(round(fuel_mass, 2))
            results["Oxidiser Mass (kg)"].append(round(oxidiser_mass, 2))
            results["Fuel Tank Height (m)"].append(round(fuel_tank_height, 2))
            results["Oxidiser Tank Height (m)"].append(round(oxidiser_tank_height, 2))
            results["Intertank Height (m)"].append(round(self.intertank_height[i], 2))
            results["Engine Count"].append(engine_count)
            results["Total Stage Height (m)"].append(round(total_height, 2))

            # Update rocket height
            self.rocket_height += total_height

        self.rocket_height += 2  # Interstage

        df = pd.DataFrame(results)
        print(df.to_markdown())
        print(f"\nTotal Rocket Height: {self.rocket_height:.2f} m")
        print("Note that the total rocket height does not include the height of the payload fairing")
        return df
    
    def simulate_trajectory(self):
        import trajectory
        #print("Simulating trajectory")
        data={'stage masses':np.array(self.stage_masses)+self.payload,
              'propellant masses':self.step_fuel_masses,
              'thrust by weight':[1.9,0.8],
              'isp':self.Isps}
        traj=trajectory.trajectory2(data)
        traj.model()
        traj.plotter()
        traj.plot_altitudes()