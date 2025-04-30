from __init__ import *
class trajectory:
    
    def __init__(self,data):
        self.position=[np.array([0,0,R_earth])]
        self.velocity=[np.array([0,0,0])]
        self.acceleration=[np.array([0,0,0])]
        self.gravity=[np.array([0,0,-g])]
        self.orientation=[np.array([0,0,1.0])]
        self.time=[0.0]
        self.data=data
        self.mass=[2e9]
        self.betas=[0.0]

    def model(self):
        i=0
        self.thrust_by_weight=self.data['thrust by weight'][0]
        self.thrust_magnitude=self.thrust_by_weight*self.mass[-1]*g
        while np.linalg.norm(self.position[-1])<(R_earth+300e3) and np.linalg.norm(self.position[-1])>(R_earth-50):
            if self.thrust_magnitude/self.mass[-1]/g>1.8:
                self.thrust_by_weight=1.8
                self.thrust_magnitude=self.thrust_by_weight*self.mass[-1]*g
            theta = self.betas[-1]
            Rx = np.array([
                        [1, 0, 0],
                        [0, np.cos(theta), -np.sin(theta)],
                        [0, np.sin(theta),  np.cos(theta)]
                        ])
            theta_rad=self.betas[-1]
            Ry = np.array([
                            [ np.cos(theta_rad), 0, np.sin(theta_rad)],
                            [ 0,                1, 0               ],
                            [-np.sin(theta_rad), 0, np.cos(theta_rad)]
                        ])
            self.thrust=Ry@(self.thrust_magnitude*(self.orientation[-1]))
            self.acceleration+=[(self.thrust)/self.mass[-1]+self.gravity[-1]]
            self.velocity+=[self.velocity[-1]+self.acceleration[-1]*dt]
            self.position+=[self.position[-1]+self.velocity[-1]*dt+self.acceleration[-1]*dt*dt/2] 
            self.time+=[self.time[-1]+dt]
            self.mass+=[self.mass[-1]-self.thrust_magnitude/self.data['average isp']/g*dt]
            self.orientation+=[self.velocity[-1]/np.linalg.norm(self.velocity[-1])]
            self.gravity+=[-g*self.position[-1]/np.linalg.norm(self.position[-1])]
            if self.time[-1]>40 and self.time[-1]<80:
                self.betas[-1]=0.1
            else:
                self.betas[-1]=0.0
            i+=1
            #print(self.position[-1][2])
        #print(self.time[-1])

    def model_sherical(self):
        i = 0
        r = R_earth  # Earth radius at surface (starting point)
        theta = 0.0 # azimuthal angle
        phi = 0.0   # polar angle (straight up)

        r_dot = 0.0
        theta_dot = 0.0
        phi_dot = 0.0

        self.thrust_by_weight = self.data['thrust by weight'][0]
        self.thrust_magnitude = self.thrust_by_weight * self.mass[-1] * g

        self.position = [[r, theta, phi]]  # now storing r, θ, φ instead of x,y,z
        self.velocity = [[r_dot, theta_dot, phi_dot]]
        self.acceleration = [[0.0, 0.0, 0.0]]

        while r < (6371e3 + 300e3) and r > 6371e3 - 50:
            if self.thrust_magnitude / self.mass[-1] / g > 4:
                self.thrust_by_weight = 4
                self.thrust_magnitude = self.thrust_by_weight * self.mass[-1] * g

            if i >= len(self.betas):
                beta = self.betas[-1]
            else:
                beta = self.betas[i]

            # Thrust split into radial and angular components
            thrust_r = self.thrust_magnitude * np.cos(beta)
            thrust_phi = self.thrust_magnitude * np.sin(beta)

            # Acceleration in r and φ (ignoring θ dynamics for now = pure vertical plane)
            r_ddot = thrust_r / self.mass[-1] + self.gravity[2] * np.cos(phi)
            phi_ddot = thrust_phi / (self.mass[-1] * r)

            # Update velocities
            r_dot += r_ddot * dt
            phi_dot += phi_ddot * dt

            # Update positions
            r += r_dot * dt + 0.5 * r_ddot * dt * dt
            phi += phi_dot * dt + 0.5 * phi_ddot * dt * dt

            # Save new states
            self.acceleration += [[r_ddot, 0, phi_ddot]]
            self.velocity += [[r_dot, 0, phi_dot]]
            self.position += [[r, theta, phi]]
            self.mass += [self.mass[-1] - self.thrust_magnitude / self.data['average isp'] / g * dt]
            self.time += [self.time[-1] + dt]

            #############################################################
            if self.time>40 and self.time<80:
                self.betas[-1]=0.1

            i += 1
            print(self.position[-1][0]-R_earth)
        print("Flight time:", self.time[-1])


    def plotter(self, position=True, velocity=False, acceleration=False,show_earth=True):
        if position:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
            x = np.array(self.position).T[0] / 1e3
            y = np.array(self.position).T[1] / 1e3
            z = np.array(self.position).T[2] / 1e3
            ax.plot3D(x, y, z)
            if show_earth:
            # --- Add Earth surface ---
                R_earth = 6371  # Earth's radius in km
                u = np.linspace(0, 2 * np.pi, 100)
                v = np.linspace(0, np.pi, 100)
                x_earth = R_earth * np.outer(np.cos(u), np.sin(v))
                y_earth = R_earth * np.outer(np.sin(u), np.sin(v))
                z_earth = R_earth * np.outer(np.ones(np.size(u)), np.cos(v))

                ax.plot_surface(x_earth, y_earth, z_earth, color='lightblue', alpha=0.5, edgecolor='none')

            # --- Equal aspect ratio ---
            max_range = np.array([x.max()-x.min(), y.max()-y.min(), z.max()-z.min()]).max() / 2.0
            mid_x = (x.max() + x.min()) * 0.5
            mid_y = (y.max() + y.min()) * 0.5
            mid_z = (z.max() + z.min()) * 0.5
            ax.set_xlim(mid_x - max_range, mid_x + max_range)
            ax.set_ylim(mid_y - max_range, mid_y + max_range)
            ax.set_zlim(mid_z - max_range, mid_z + max_range)

            output_file = "output files/TRC files/planned trajectory.png"
            fig.savefig(output_file)
            plt.show()

        if velocity:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
            x = np.array(self.velocity).T[0]
            y = np.array(self.velocity).T[1]
            z = np.array(self.velocity).T[2]
            ax.plot3D(x, y, z)
            if show_earth:
            # --- Add Earth surface ---
                R_earth = 6371  # Earth's radius in km
                u = np.linspace(0, 2 * np.pi, 100)
                v = np.linspace(0, np.pi, 100)
                x_earth = R_earth * np.outer(np.cos(u), np.sin(v))
                y_earth = R_earth * np.outer(np.sin(u), np.sin(v))
                z_earth = R_earth * np.outer(np.ones(np.size(u)), np.cos(v))

                ax.plot_surface(x_earth, y_earth, z_earth, color='lightblue', alpha=0.5, edgecolor='none')

            # --- Equal aspect ratio ---
            max_range = np.array([x.max()-x.min(), y.max()-y.min(), z.max()-z.min()]).max() / 2.0
            mid_x = (x.max() + x.min()) * 0.5
            mid_y = (y.max() + y.min()) * 0.5
            mid_z = (z.max() + z.min()) * 0.5
            ax.set_xlim(mid_x - max_range, mid_x + max_range)
            ax.set_ylim(mid_y - max_range, mid_y + max_range)
            ax.set_zlim(mid_z - max_range, mid_z + max_range)

            output_file = "output files/TRC files/planned trajectory.png"
            fig.savefig(output_file)
            plt.show()

    def plotter_spherical(self, position=True, velocity=False, acceleration=False):
        if position:
            # Extract r, theta, phi from stored positions
            r = np.array(self.position).T[0]
            theta = np.array(self.position).T[1]
            phi = np.array(self.position).T[2]

            # Convert to Cartesian coordinates for plotting
            x = r * np.sin(phi) * np.cos(theta)
            y = r * np.sin(phi) * np.sin(theta)
            z = r * np.cos(phi)

            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')

            # Plot the trajectory
            ax.plot3D(x / 1e3, y / 1e3, z / 1e3)  # plot in km

            ax.set_xlabel('X (km)')
            ax.set_ylabel('Y (km)')
            ax.set_zlabel('Z (km)')
            ax.set_title('Trajectory in 3D Space (converted from Spherical)')

            output_file = "output files/TRC files/planned_trajectory_spherical.png"
            fig.savefig(output_file)
            plt.show()

        if velocity:
            # Optional: plot velocity vectors (would need similar conversion if storing angular rates)
            pass


    def dyn_plotter(self):
        # Create animation
        self.positions=np.array(self.position)/1e3
        self.fig = plt.figure()
        self.ax = self.fig.add_subplot(111, projection='3d')
        self.ax.set_xlim(np.min(self.positions[:,0]), np.max(self.positions[:,0]))
        self.ax.set_ylim(np.min(self.positions[:,1]), np.max(self.positions[:,1]))
        self.ax.set_zlim(np.min(self.positions[:,2]), np.max(self.positions[:,2]))
        self.line, = self.ax.plot([], [], [], lw=2, color='blue')
        self.point, = self.ax.plot([], [], [], 'ro')
        ani = FuncAnimation(self.fig, self.update, frames=len(self.position), init_func=self.init,
                    interval=dt, blit=True)
        plt.show()

    def init(self):
        self.line.set_data([], [])
        self.line.set_3d_properties([])
        self.point.set_data([], [])
        self.point.set_3d_properties([])
        return self.line, self.point

    # Update function
    def update(self,frame):
        x = self.positions[:frame, 0]
        y = self.positions[:frame, 1]
        z = self.positions[:frame, 2]
        self.line.set_data(x, y)
        self.line.set_3d_properties(z)
    
        self.point.set_data(x[-1:], y[-1:])
        self.point.set_3d_properties(z[-1:])
        return self.line, self.point

def objective(betas_flat):
        betas = list(betas_flat)
        traj = trajectory(data)
        traj.betas = betas + [betas[-1]] * (len(traj.position) - len(betas))  # pad if needed
        traj.model()
        return traj.time[-1] if traj.position[-1][2] > 0 else 1e6  # penalty if flight fails

def optimize_betas(data, num_steps=100):
        initial_guess = np.zeros(num_steps)  # start with all zeros
        bounds = [(0, 0.1)] * num_steps  # limit tilt angles
        print("optimising betas")
        result = differential_evolution(objective, bounds,
                                    maxiter=10, popsize=1)
        if result.success:
            print(f"Optimal flight time: {result.fun:.2f} s")
            return result.x
        else:
            print("Optimization failed:", result.message)
            return initial_guess

###########################################################################################################################

class trajectory2:
    
    def __init__(self,data):
        #self.position=[np.array([0,0,R_earth])
        self.position=[np.array([R_earth*np.cos(latitude),0,R_earth*np.sin(latitude)])]
        self.rotation_velocity=np.cross(np.array([0,0,w_earth]),self.position)[0]
        self.velocity=[np.array([0,0,0])]
        #self.velocity[-1]=self.rotation_velocity
        self.acceleration=[np.array([0,0,0])]
        #self.gravity=np.array([0,0,-g])
        self.gravity=-g*np.array([np.cos(latitude),0,np.sin(latitude)])
        self.orientation=[np.array([np.cos(latitude),0,np.sin(latitude)])]
        self.time=[0.0]
        self.data=data
        self.betas=[0.0]
        self.stage_masses=data['stage masses']
        self.stage_propellant_masses=data['propellant masses']
        self.thrust_by_weights=data['thrust by weight']
        self.number_of_stages=len(self.stage_masses)
        self.burnouts=[]

    def model(self,post_burnout=True):
        i=0
        for current_stage in range(0,self.number_of_stages):

            self.thrust_by_weight=self.data['thrust by weight'][current_stage]
            self.thrust_magnitude=self.thrust_by_weight*self.stage_masses[current_stage]*g
            while self.stage_propellant_masses[current_stage]>0:
                if self.thrust_magnitude/self.stage_masses[current_stage]/g>4:
                    self.thrust_by_weight=4
                    self.thrust_magnitude=self.thrust_by_weight*self.stage_masses[current_stage]*g
                theta = -self.betas[-1]
                Rx = np.array([
                            [1, 0, 0],
                            [0, np.cos(theta), -np.sin(theta)],
                            [0, np.sin(theta),  np.cos(theta)]
                            ])
                theta_rad=-self.betas[-1]
                Ry = np.array([
                                [ np.cos(theta_rad), 0, np.sin(theta_rad)],
                                [ 0,                1, 0               ],
                                [-np.sin(theta_rad), 0, np.cos(theta_rad)]
                            ])
                self.thrust=Rx@(self.thrust_magnitude*(self.orientation[-1]))
                self.acceleration+=[(self.thrust)/self.stage_masses[current_stage]+self.gravity]
                self.velocity+=[self.velocity[-1]+self.acceleration[-1]*dt]
                self.position+=[self.position[-1]+self.velocity[-1]*dt+self.acceleration[-1]*dt*dt/2+self.rotation_velocity*dt] 
                self.time+=[self.time[-1]+dt]
                #self.mass+=[self.mass[-1]-self.thrust_magnitude/self.data['average isp']/g*dt]
                self.stage_masses[current_stage]-=self.thrust_magnitude/self.data['isp'][current_stage]/g*dt
                self.stage_propellant_masses[current_stage]-=self.thrust_magnitude/self.data['isp'][current_stage]/g*dt
                self.orientation+=[self.velocity[-1]/np.linalg.norm(self.velocity[-1])]
                self.gravity=-g*self.position[-1]/np.linalg.norm(self.position[-1])
                if self.time[-1]>15 and self.time[-1]<200:
                    self.betas[-1]=0.25
                else:
                    self.betas[-1]=0.0
                i+=1
                #if i%10==0:
                    #print(self.data['isp'][current_stage]*g*dt)
            print(f"Stage {current_stage+1} burnout at {round((np.linalg.norm(self.position[-1])-R_earth)/1e3,2)} km at {round(self.time[-1]/60,2)} minutes")
            self.burnouts+=[self.position[-1]]
            #print(self.position[-1][2])
        #print(np.linalg.norm(self.velocity[-1]))
        self.burnouts=np.array(self.burnouts)
        #print(self.acceleration[1],self.velocity[1])
        print("Note that the burnout altitudes are not very accurate as they depend on the actual ascent trajectory which is not yet fixed")
        if post_burnout:
            while i<240*60*10:
                self.acceleration+=[self.gravity]
                self.velocity+=[self.velocity[-1]+self.acceleration[-1]*dt]
                self.position+=[self.position[-1]+self.velocity[-1]*dt+self.acceleration[-1]*dt*dt/2+self.rotation_velocity*dt] 
                self.time+=[self.time[-1]+dt]
                self.orientation+=[self.velocity[-1]/np.linalg.norm(self.velocity[-1])]
                self.gravity=-g*self.position[-1]/np.linalg.norm(self.position[-1])
                i+=1
    
    def plotter(self, position=True, velocity=False, acceleration=False,show_earth=True):
        if position:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
            x = np.array(self.position).T[0] / 1e3
            y = np.array(self.position).T[1] / 1e3
            z = np.array(self.position).T[2] / 1e3
            ax.plot3D(x, y, z)

            x_burnout = self.burnouts[:, 0] / 1e3  # in km
            y_burnout = self.burnouts[:, 1] / 1e3
            z_burnout = self.burnouts[:, 2] / 1e3
            ax.scatter(x_burnout, y_burnout, z_burnout, color='red', s=20, label='Burnout Points')


            if show_earth:
            # --- Add Earth surface ---
                R_earth = 6371  # Earth's radius in km
                u = np.linspace(0, 2 * np.pi, 100)
                v = np.linspace(0, np.pi, 100)
                x_earth = R_earth * np.outer(np.cos(u), np.sin(v))
                y_earth = R_earth * np.outer(np.sin(u), np.sin(v))
                z_earth = R_earth * np.outer(np.ones(np.size(u)), np.cos(v))

                ax.plot_surface(x_earth, y_earth, z_earth, color='lightblue', alpha=0.5, edgecolor='none')

            # --- Equal aspect ratio ---
            max_range = np.array([x.max()-x.min(), y.max()-y.min(), z.max()-z.min()]).max() / 2.0
            mid_x = (x.max() + x.min()) * 0.5
            mid_y = (y.max() + y.min()) * 0.5
            mid_z = (z.max() + z.min()) * 0.5
            ax.set_xlim(mid_x - max_range, mid_x + max_range)
            ax.set_ylim(mid_y - max_range, mid_y + max_range)
            ax.set_zlim(mid_z - max_range, mid_z + max_range)

            ax.set_xlabel('X (km)')
            ax.set_ylabel('Y (km)')
            ax.set_zlabel('Z (km)')            

            output_file = "output files/TRC files/planned trajectory.png"
            fig.savefig(output_file)
            plt.show()

        if velocity:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
            x = np.array(self.velocity).T[0]
            y = np.array(self.velocity).T[1]
            z = np.array(self.velocity).T[2]
            ax.plot3D(x, y, z)
            if show_earth:
            # --- Add Earth surface ---
                R_earth = 6371  # Earth's radius in km
                u = np.linspace(0, 2 * np.pi, 100)
                v = np.linspace(0, np.pi, 100)
                x_earth = R_earth * np.outer(np.cos(u), np.sin(v))
                y_earth = R_earth * np.outer(np.sin(u), np.sin(v))
                z_earth = R_earth * np.outer(np.ones(np.size(u)), np.cos(v))

                ax.plot_surface(x_earth, y_earth, z_earth, color='lightblue', alpha=0.5, edgecolor='none')

            # --- Equal aspect ratio ---
            max_range = np.array([x.max()-x.min(), y.max()-y.min(), z.max()-z.min()]).max() / 2.0
            mid_x = (x.max() + x.min()) * 0.5
            mid_y = (y.max() + y.min()) * 0.5
            mid_z = (z.max() + z.min()) * 0.5
            ax.set_xlim(mid_x - max_range, mid_x + max_range)
            ax.set_ylim(mid_y - max_range, mid_y + max_range)
            ax.set_zlim(mid_z - max_range, mid_z + max_range)

            ax.set_xlabel('X (km)')
            ax.set_ylabel('Y (km)')
            ax.set_zlabel('Z (km)')

            output_file = "output files/TRC files/planned trajectory.png"
            fig.savefig(output_file)
            plt.show()


    def dyn_plotter(self,show_earth=True):
        # Create animation
        x = np.array(self.position).T[0] / 1e3
        y = np.array(self.position).T[1] / 1e3
        z = np.array(self.position).T[2] / 1e3
        self.positions=np.array(self.position)/1e3
        self.fig = plt.figure()
        self.ax = self.fig.add_subplot(111, projection='3d')
        self.ax.set_xlim(np.min(self.positions[:,0]), np.max(self.positions[:,0]))
        self.ax.set_ylim(np.min(self.positions[:,1]), np.max(self.positions[:,1]))
        self.ax.set_zlim(np.min(self.positions[:,2]), np.max(self.positions[:,2]))
        if show_earth:
            # --- Add Earth surface ---
            R_earth = 6371  # Earth's radius in km
            u = np.linspace(0, 2 * np.pi, 100)
            v = np.linspace(0, np.pi, 100)
            x_earth = R_earth * np.outer(np.cos(u), np.sin(v))
            y_earth = R_earth * np.outer(np.sin(u), np.sin(v))
            z_earth = R_earth * np.outer(np.ones(np.size(u)), np.cos(v))
            self.ax.plot_surface(x_earth, y_earth, z_earth, color='lightblue', alpha=0.5, edgecolor='none')
            max_range = np.array([x.max()-x.min(), y.max()-y.min(), z.max()-z.min()]).max() / 2.0
            mid_x = (x.max() + x.min()) * 0.5
            mid_y = (y.max() + y.min()) * 0.5
            mid_z = (z.max() + z.min()) * 0.5
            self.ax.set_xlim(mid_x - max_range, mid_x + max_range)
            self.ax.set_ylim(mid_y - max_range, mid_y + max_range)
            self.ax.set_zlim(mid_z - max_range, mid_z + max_range)
        
        self.line, = self.ax.plot([], [], [], lw=2, color='blue')
        self.point, = self.ax.plot([], [], [], 'ro')
        ani = FuncAnimation(self.fig, self.update, frames=len(self.position), init_func=self.init,
                    interval=dt, blit=True)
        

        x_burnout = self.burnouts[:, 0] / 1e3  # in km
        y_burnout = self.burnouts[:, 1] / 1e3
        z_burnout = self.burnouts[:, 2] / 1e3
        self.ax.scatter(x_burnout, y_burnout, z_burnout, color='red', s=20, label='Burnout Points')

        self.ax.set_xlabel('X (km)')
        self.ax.set_ylabel('Y (km)')
        self.ax.set_zlabel('Z (km)')

        plt.show()

    def plot_altitudes(self):
        plt.figure()
        self.altitudes=[(np.linalg.norm(self.position[i])-R_earth)/1e3 for i in range(len(self.position))]
        plt.plot(self.time,self.altitudes)
        plt.title("Altitude vs Time")
        plt.show()
        plt.figure()
        self.velocities=[(np.linalg.norm(self.velocity[i])) for i in range(len(self.velocity))]
        plt.plot(self.time,self.velocities)
        plt.title("Velocity vs Time")
        plt.show()

    def init(self):
        self.line.set_data([], [])
        self.line.set_3d_properties([])
        self.point.set_data([], [])
        self.point.set_3d_properties([])
        return self.line, self.point

    # Update function
    def update(self,frame):
        x = self.positions[:frame, 0]
        y = self.positions[:frame, 1]
        z = self.positions[:frame, 2]
        self.line.set_data(x, y)
        self.line.set_3d_properties(z)
    
        self.point.set_data(x[-1:], y[-1:])
        self.point.set_3d_properties(z[-1:])
        return self.line, self.point

def objective(betas_flat):
        betas = list(betas_flat)
        traj = trajectory(data)
        traj.betas = betas + [betas[-1]] * (len(traj.position) - len(betas))  # pad if needed
        traj.model()
        return traj.time[-1] if traj.position[-1][2] > 0 else 1e6  # penalty if flight fails

def optimize_betas(data, num_steps=100):
        initial_guess = np.zeros(num_steps)  # start with all zeros
        bounds = [(0, 0.1)] * num_steps  # limit tilt angles
        print("optimising betas")
        result = differential_evolution(objective, bounds,
                                    maxiter=10, popsize=1)
        if result.success:
            print(f"Optimal flight time: {result.fun:.2f} s")
            return result.x
        else:
            print("Optimization failed:", result.message)
            return initial_guess


'''

input_path = "input files/" + "TRC" + ".json"
with open(input_path) as file:
    data = json.load(file)
print(data['thrust by weight'][0])
#optimal_betas = optimize_betas(data, num_steps=5000)
traj=trajectory(data)
#traj.betas = list(optimal_betas)
traj.model()
print(traj.time[-1])
#traj.plotter()                  




input_path = "input files/" + "TRC" + ".json"
with open(input_path) as file:
    data = json.load(file)
print(data['thrust by weight'][0])
#optimal_betas = optimize_betas(data, num_steps=5000)
traj=trajectory(data)
#traj.betas = list(optimal_betas)
traj.model()
traj.plotter()                  

'''