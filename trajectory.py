# class trajectory is being used for the simulation

from __init__ import *
class trajectory:
    
    def __init__(self,data):
        #self.position=[np.array([0,0,R_earth])
        self.position=[np.array([R_earth*np.cos(latitude),0,R_earth*np.sin(latitude)])]
        self.omega_earth=np.array([0,0,w_earth])
        self.rotation_velocity=np.cross(np.array([0,0,w_earth]),self.position)[0]
        self.velocity=[np.array([0,0,0])]
        #self.velocity[-1]=self.rotation_velocity
        self.acceleration=[np.array([0,0,0])]
        self.aoa=[0]
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
        self.velocity_inertial_frame=[self.velocity[-1]+np.cross(self.omega_earth,self.position[-1])]
        self.acceleration_inertial_frame=[self.acceleration[-1]+2*np.cross(self.omega_earth,self.velocity[-1])+np.cross(self.omega_earth,np.cross(self.omega_earth,self.position[-1]))]
        self.position_inertial_frame=[self.position[-1]]
        self.burnouts_inertial_frame=[]
        self.burnout_times=[]

    def model(self,post_burnout=True,dt=dt,vector_end_time=200,beta_max=0.2):
        i=0
        out_str=f"\nRocket liftoff from {np.round(self.position[0]/1000,2)} km. Earth rotates along the Z-axis"
        for current_stage in range(0,self.number_of_stages):
            delta_v_achieved=0
            self.thrust_by_weight=self.data['thrust by weight'][current_stage]
            self.thrust_magnitude=self.thrust_by_weight*self.stage_masses[current_stage]*g

            while self.stage_propellant_masses[current_stage]>0:
                self.betas+=[0]
                if self.thrust_magnitude/self.stage_masses[current_stage]/g>4:
                    self.thrust_by_weight=4
                    self.thrust_magnitude=self.thrust_by_weight*self.stage_masses[current_stage]*g
                
                theta = -self.betas[i]
                Rx = np.array([
                            [1, 0, 0],
                            [0, np.cos(theta), -np.sin(theta)],
                            [0, np.sin(theta),  np.cos(theta)]
                            ])
                
                theta_rad=-self.betas[i]
                
                Ry = np.array([
                                [ np.cos(theta_rad), 0, np.sin(theta_rad)],
                                [ 0,                1, 0               ],
                                [-np.sin(theta_rad), 0, np.cos(theta_rad)]
                            ])
                
                theta_earth=w_earth*self.time[-1]
                Rz = np.array([
                                    [np.cos(theta_earth), -np.sin(theta_earth), 0],
                                    [np.sin(theta_earth),  np.cos(theta_earth), 0],
                                    [0,                    0,                   1]
                            ])


                self.thrust=Rz@(self.thrust_magnitude*(self.orientation[-1]))
                coriolis=-2*np.cross(self.omega_earth,self.velocity[-1])
                centrifugal=-np.cross(self.omega_earth,np.cross(self.omega_earth,self.position[-1]))

                thrust_inertial_frame=Rz@self.thrust

                self.acceleration_inertial_frame+=[(thrust_inertial_frame)/self.stage_masses[current_stage]]

                self.acceleration+=[self.acceleration_inertial_frame[-1]+self.gravity+coriolis+centrifugal]
                if np.linalg.norm(self.velocity[-1])>0:
                    delta_v_achieved+=np.dot((self.thrust)/self.stage_masses[current_stage]*dt,self.velocity[-1])/np.linalg.norm(self.velocity[-1])
                
                self.acceleration_inertial_frame+=[self.acceleration[-1]+2*np.cross(self.omega_earth,self.velocity[-1])+np.cross(self.omega_earth,np.cross(self.omega_earth,self.position[-1]))]
                
                self.velocity+=[self.velocity[-1]+self.acceleration[-1]*dt]
                self.position+=[self.position[-1]+self.velocity[-1]*dt+self.acceleration[-1]*dt*dt/2]#+self.rotation_velocity*dt] 
                self.time+=[self.time[-1]+dt]
                self.orientation+=[self.velocity[-1]/np.linalg.norm(self.velocity[-1])]
                self.gravity=-g*self.position[-1]/np.linalg.norm(self.position[-1])

                self.velocity_inertial_frame+=[self.velocity[-1]+np.cross(self.omega_earth,self.position[-1])]
                
                self.position_inertial_frame+=[self.position_inertial_frame[-1]+self.velocity_inertial_frame[-1]*dt+self.acceleration_inertial_frame[-1]*dt*dt/2]

                self.stage_masses[current_stage]-=self.thrust_magnitude/self.data['isp'][current_stage]/g*dt
                self.stage_propellant_masses[current_stage]-=self.thrust_magnitude/self.data['isp'][current_stage]/g*dt
                self.orientation+=[self.velocity[-1]/np.linalg.norm(self.velocity[-1])]
                self.gravity=-g*self.position[-1]/np.linalg.norm(self.position[-1])

                if self.time[-1]>15 and self.time[-1]<vector_end_time:
                    self.betas[i+1]=beta_max
                    #print(i,self.betas[i])
                else:
                    self.betas[i+1]=0
                i+=1
                
            #print()
            out_str+=f"Stage {current_stage+1} burnout:\nat {round((np.linalg.norm(self.position[-1])-R_earth)/1e3,2)} km at {round(self.time[-1]/60,2)} minutes\n"
            #print(f"Stage {current_stage+1} burnout:\nat {round((np.linalg.norm(self.position[-1])-R_earth)/1e3,2)} km at {round(self.time[-1]/60,2)} minutes")
            self.burnouts+=[self.position[-1]]
            self.burnouts_inertial_frame+=[self.position_inertial_frame[-1]]
            self.burnout_times+=[self.time[-1]]
            
            theta_=w_earth*self.time[-1]
            Rz = np.array([
                                [ np.cos(theta_), -np.sin(theta_), 0],
                                [np.sin(theta_), np.cos(theta_), 0],
                                [0,             0,                1]
                            ])
            
            out_str+=f"downrange position from launch site (in earth centered earth fixed frame) = {np.round((self.position[-1]-self.position[0])/1e3,2)} km\n"
            out_str+=f"velocity vector at burnout (in earth centered earth fixed frame) = {np.round(self.velocity[-1],2)} m/s\n"
            out_str+=f"acceleration vector at burnout (in earth centered earth fixed frame) = {np.round(self.acceleration[-1],2)} m/s^2\n"
            out_str+=f"direction vector at burnout (in earth centered earth fixed frame) = {np.round(self.orientation[-1],4)}\n"
            out_str+=f"downrange position from launch site (in earth centered inertial frame) = {np.round((self.position_inertial_frame[-1]-Rz@self.position_inertial_frame[0])/1e3,2)} km\n"
            out_str+=f"position of launch site at burnout (in earth centered inertial frame) = {np.round((Rz@self.position_inertial_frame[0])/1e3,2)} km\n"
            out_str+=f"velocity vector at burnout (in earth centered inertial frame) = {np.round(self.velocity_inertial_frame[-1],2)} m/s\n"
            out_str+=f"acceleration vector at burnout (in earth centered inertial frame) = {np.round(self.acceleration_inertial_frame[-1],2)} m/s\n"
            out_str+=f"direction vector at burnout (in earth centered inertial frame) = {np.round(Rz@self.orientation[-1],4)}\n"
            out_str+=f"delta-v achieved = {delta_v_achieved} m/s\n"
            out_str+="\n"
            
            """print(f"downrange position from launch site (in earth centered earth fixed frame) = {np.round((self.position[-1]-self.position[0])/1e3,2)} km")
            print(f"velocity vector at burnout (in earth centered earth fixed frame) = {np.round(self.velocity[-1],2)} m/s")
            print(f"acceleration vector at burnout (in earth centered earth fixed frame) = {np.round(self.acceleration[-1],2)} m/s^2")
            print(f"direction vector at burnout (in earth centered earth fixed frame) = {np.round(self.orientation[-1],4)}")
            print(f"downrange position from launch site (in earth centered inertial frame) = {np.round((self.position_inertial_frame[-1]-Rz@self.position_inertial_frame[0])/1e3,2)} km")
            print(f"position of launch site at burnout (in earth centered inertial frame) = {np.round((Rz@self.position_inertial_frame[0])/1e3,2)} km")
            print(f"velocity vector at burnout (in earth centered inertial frame) = {np.round(self.velocity_inertial_frame[-1],2)} m/s")
            print(f"acceleration vector at burnout (in earth centered inertial frame) = {np.round(self.acceleration_inertial_frame[-1],2)} m/s")
            print(f"direction vector at burnout (in earth centered inertial frame) = {np.round(Rz@self.orientation[-1],4)}")
            print()"""
            #print(f"direction vector at burnout (in earth centered inertial frame) = {np.round(self.orientation_inertial_frame[-1],2)} m/s")
            #print(self.position[-1][2])
        #print(np.linalg.norm(self.velocity[-1]))
        self.burnouts=np.array(self.burnouts)
        self.burnouts_inertial_frame=np.array(self.burnouts_inertial_frame)
        #print(self.acceleration[1],self.velocity[1])
        out_str+="Note that the burnout altitudes are not very accurate as they depend on the actual ascent trajectory which is not yet fixed\n"
        #print("Note that the burnout altitudes are not very accurate as they depend on the actual ascent trajectory which is not yet fixed")

        if post_burnout:
            while self.time[-1]<total_time:
                coriolis=-2*np.cross(self.omega_earth,self.velocity[-1])
                centrifugal=-np.cross(self.omega_earth,np.cross(self.omega_earth,self.position[-1]))
                self.acceleration+=[self.gravity+coriolis+centrifugal]
                self.velocity+=[self.velocity[-1]+self.acceleration[-1]*dt]
                self.position+=[self.position[-1]+self.velocity[-1]*dt+self.acceleration[-1]*dt*dt/2] 
                self.time+=[self.time[-1]+dt]
                self.orientation+=[self.velocity[-1]/np.linalg.norm(self.velocity[-1])]
                self.gravity=-g*self.position[-1]/np.linalg.norm(self.position[-1])


                self.velocity_inertial_frame+=[self.velocity[-1]+np.cross(self.omega_earth,self.position[-1])]
                self.acceleration_inertial_frame+=[self.acceleration[-1]+2*np.cross(self.omega_earth,self.velocity[-1])+np.cross(self.omega_earth,np.cross(self.omega_earth,self.position[-1]))]
                self.position_inertial_frame+=[self.position_inertial_frame[-1]+self.velocity_inertial_frame[-1]*dt+self.acceleration_inertial_frame[-1]*dt*dt/2]

                i+=1
        return out_str
    
    def plotter_inertial_frame(self, position=True, velocity=False, show_earth=True):
        if position:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
            x = np.array(self.position_inertial_frame).T[0] / 1e3
            y = np.array(self.position_inertial_frame).T[1] / 1e3
            z = np.array(self.position_inertial_frame).T[2] / 1e3
            ax.plot3D(x, y, z)

            x_burnout = self.burnouts_inertial_frame[:, 0] / 1e3  # in km
            y_burnout = self.burnouts_inertial_frame[:, 1] / 1e3
            z_burnout = self.burnouts_inertial_frame[:, 2] / 1e3
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

            
            plt.show()

    def plotter(self,show_earth=True):
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

            output_file = f"output files/{self.data["rocket name"]} files/3d trajectory.png"
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

    def plot_altitudes_inertial_frame(self,show_plots=False):
        plt.figure()
        self.altitudes=[(np.linalg.norm(self.position_inertial_frame[i])-R_earth)/1e3 for i in range(len(self.position_inertial_frame))]
        self.altitudes_burnout_inertial_frame=[(np.linalg.norm(self.burnouts_inertial_frame[i])-R_earth)/1e3 for i in range(len(self.burnouts_inertial_frame))]
        plt.plot(self.time,self.altitudes)
        plt.scatter(self.burnout_times,self.altitudes_burnout_inertial_frame)
        plt.title("Altitude vs Time")
        output_file = f"output files/{self.data["rocket name"]} files/altitude vs time in inertial frame.png"
        plt.savefig(output_file)
        if show_plots:
            plt.show()
        plt.figure()
        self.velocities=[(np.linalg.norm(self.velocity_inertial_frame[i])) for i in range(len(self.velocity_inertial_frame))]
        plt.plot(self.time,self.velocities)
        plt.title("Velocity vs Time")
        output_file = f"output files/{self.data["rocket name"]} files/velocity vs time in inertial frame.png"
        plt.savefig(output_file)
        if show_plots:
            plt.show()

    def plot_altitudes(self,show_plots=False):
        plt.figure()
        self.altitudes=[(np.linalg.norm(self.position[i])-R_earth)/1e3 for i in range(len(self.position))]
        plt.plot(self.time,self.altitudes)
        plt.title("Altitude vs Time")
        if show_plots:
            plt.show()
        plt.figure()
        self.velocities=[(np.linalg.norm(self.velocity[i])) for i in range(len(self.velocity))]
        plt.plot(self.time,self.velocities)
        plt.title("Velocity vs Time")
        if show_plots:
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

    def cost_function(self, betas_array):
        self.__init__(self.data)  # reset the simulation
        self.betas = list(betas_array)
        self.model()
        final_altitude = []
        for i in range(1,201):
            final_altitude+=[np.linalg.norm(self.position[-i]) - R_earth]
        final_altitude=np.array(final_altitude)
        
        print("Current loss = ",np.sum(abs(final_altitude - self.data["target orbit"]))," Current final altitude = ",np.linalg.norm(self.position[-1]) - R_earth)
        return np.sum(abs(final_altitude - self.data["target orbit"]*1e3))

    def optimize_betas(self, steps=total_time/dt):
        initial_guess = [0.0] * (int(steps))
        bounds = [(0, 0.1)] * (int(steps))
        print("Beta array length = ",len(initial_guess))
        result = minimize(self.cost_function, initial_guess, bounds=bounds, method="L-BFGS-B")
        self.betas = list(result.x)
        return result.x#, result.fun
    
def optimize_betas_2(data,t1=20,t2=100,n_ts=20,beta_max=0.2,n_betas=20):
        print("Optimising vector end time")
        plt.figure()
        min_loss=1e100
        ans1=0
        ans2=0
        count=1
        for beta in np.linspace(0.05,beta_max,n_betas):
            for t in np.linspace(t1,t2,n_ts):
                traj=trajectory(copy.deepcopy(data))  # reset the simulation
                traj.model(vector_end_time=t,beta_max=beta)
                traj.altitudes=[(np.linalg.norm(traj.position_inertial_frame[i])-R_earth)/1e3 for i in range(len(traj.position_inertial_frame))]
                traj.altitudes_burnout_inertial_frame=[(np.linalg.norm(traj.burnouts_inertial_frame[i])-R_earth)/1e3 for i in range(len(traj.burnouts_inertial_frame))]
                plt.plot(traj.time,traj.altitudes)
                plt.scatter(traj.burnout_times,traj.altitudes_burnout_inertial_frame)
                plt.title("Altitude vs Time")
                final_altitude = []
                for i in range(1,2001):
                    final_altitude+=[np.linalg.norm(traj.position[-i]) - R_earth]
                final_altitude=np.array(final_altitude)
                current_loss=np.average(np.abs(final_altitude - traj.data["target orbit"]))/1e3
                print("Current loss = ",current_loss," Current final altitude = ",int(np.linalg.norm(traj.position[-1]) - R_earth)/1e3)
                if current_loss<min_loss:
                    ans1=t
                    ans2=beta
                    min_loss=current_loss
                print(np.round(count/n_betas/n_ts*100,2),"% completed")
                count+=1
        plt.show() 
        return ans1,ans2

###########################################################################################################################

class trajectory2():
    
    def __init__(self,data):
        #self.position=[np.array([0,0,R_earth])
        self.position=[np.array([R_earth*np.cos(latitude),0,R_earth*np.sin(latitude)])]
        self.omega_earth=np.array([0,0,w_earth])
        self.rotation_velocity=np.cross(np.array([0,0,w_earth]),self.position)[0]
        self.velocity=[np.array([0,0,0])]
        #self.velocity[-1]=self.rotation_velocity
        self.acceleration=[np.array([0,0,0])]
        self.aoa=[0]
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
        self.velocity_inertial_frame=[self.velocity[-1]+np.cross(self.omega_earth,self.position[-1])]
        self.acceleration_inertial_frame=[self.acceleration[-1]+2*np.cross(self.omega_earth,self.velocity[-1])+np.cross(self.omega_earth,np.cross(self.omega_earth,self.position[-1]))]
        self.position_inertial_frame=[self.position[-1]]
        self.burnouts_inertial_frame=[]
        self.burnout_times=[]

    def model(self,post_burnout=True,dt=dt,vector_end_time=200,beta_max=0.2,total_time=total_time):
        i=0
        out_str=f"\nRocket liftoff from {np.round(self.position[0]/1000,2)} km. Earth rotates along the Z-axis\n"
        for current_stage in range(0,self.number_of_stages):
            delta_v_achieved=0
            self.thrust_by_weight=self.data['thrust by weight'][current_stage]
            self.thrust_magnitude=self.thrust_by_weight*self.stage_masses[current_stage]*g

            while self.stage_propellant_masses[current_stage]>0:
                self.betas+=[0]
                if self.thrust_magnitude/self.stage_masses[current_stage]/g>4:
                    self.thrust_by_weight=4
                    self.thrust_magnitude=self.thrust_by_weight*self.stage_masses[current_stage]*g
                
                theta = self.betas[i]
                Rx = np.array([
                            [1, 0, 0],
                            [0, np.cos(theta), -np.sin(theta)],
                            [0, np.sin(theta),  np.cos(theta)]
                            ])
                
                theta_rad = self.betas[i]
                
                Ry = np.array([
                                [ np.cos(theta_rad), 0, np.sin(theta_rad)],
                                [ 0,                1, 0               ],
                                [-np.sin(theta_rad), 0, np.cos(theta_rad)]
                            ])
                
                Rz = np.array([
                            [np.cos(theta), -np.sin(theta), 0],
                            [np.sin(theta),  np.cos(theta), 0],
                            [0            , 0             , 1]
                            ])
                
                self.thrust=Rz@(self.thrust_magnitude*(self.orientation[-1]))
                coriolis=-2*np.cross(self.omega_earth,self.velocity[-1])
                centrifugal=-np.cross(self.omega_earth,np.cross(self.omega_earth,self.position[-1]))

                self.acceleration+=[(self.thrust)/self.stage_masses[current_stage]+self.gravity+coriolis+centrifugal]
                if np.linalg.norm(self.velocity[-1])>0:
                    delta_v_achieved+=np.dot((self.thrust)/self.stage_masses[current_stage]*dt,self.velocity[-1])/np.linalg.norm(self.velocity[-1])
                
                self.acceleration_inertial_frame+=[self.acceleration[-1]+2*np.cross(self.omega_earth,self.velocity[-1])+np.cross(self.omega_earth,np.cross(self.omega_earth,self.position[-1]))]
                
                self.velocity+=[self.velocity[-1]+self.acceleration[-1]*dt]
                self.position+=[self.position[-1]+self.velocity[-1]*dt+self.acceleration[-1]*dt*dt/2]#+self.rotation_velocity*dt] 
                self.time+=[self.time[-1]+dt]
                self.orientation+=[self.velocity[-1]/np.linalg.norm(self.velocity[-1])]
                self.gravity=-g*self.position[-1]/np.linalg.norm(self.position[-1])

                self.velocity_inertial_frame+=[self.velocity[-1]+np.cross(self.omega_earth,self.position[-1])]
                
                theta_=w_earth*self.time[-1]
                Rz = np.array([
                                    [ np.cos(theta_), -np.sin(theta_), 0],
                                    [np.sin(theta_), np.cos(theta_), 0],
                                    [0,             0,                1]
                                ])

                #self.position_inertial_frame+=[self.position_inertial_frame[-1]+self.velocity_inertial_frame[-1]*dt+self.acceleration_inertial_frame[-1]*dt*dt/2]
                self.position_inertial_frame+=[Rz@self.position[-1]]
                self.stage_masses[current_stage]-=self.thrust_magnitude/self.data['isp'][current_stage]/g*dt
                self.stage_propellant_masses[current_stage]-=self.thrust_magnitude/self.data['isp'][current_stage]/g*dt
                self.orientation+=[self.velocity[-1]/np.linalg.norm(self.velocity[-1])]
                self.gravity=-g*self.position[-1]/np.linalg.norm(self.position[-1])

                if self.time[-1]>15 and self.time[-1]<vector_end_time:
                    self.betas[i+1]=beta_max
                    #print(i,self.betas[i])
                else:
                    self.betas[i+1]=0
                i+=1
                #out_str+=f"{self.time[-1]}, {self.velocity[-1]}, {self.acceleration[-1]}, {coriolis}, {centrifugal}\n"
                #print(self.stage_propellant_masses[current_stage])
            #print()
            out_str+=f"Stage {current_stage+1} burnout:\nat {round((np.linalg.norm(self.position[-1])-R_earth)/1e3,2)} km at {round(self.time[-1]/60,2)} minutes\n"
            #print(f"Stage {current_stage+1} burnout:\nat {round((np.linalg.norm(self.position[-1])-R_earth)/1e3,2)} km at {round(self.time[-1]/60,2)} minutes")
            self.burnouts+=[self.position[-1]]
            self.burnouts_inertial_frame+=[self.position_inertial_frame[-1]]
            self.burnout_times+=[self.time[-1]]
            
            theta_=w_earth*self.time[-1]
            Rz = np.array([
                                [ np.cos(theta_), -np.sin(theta_), 0],
                                [np.sin(theta_), np.cos(theta_), 0],
                                [0,             0,                1]
                            ])
            
            out_str+=f"downrange position from launch site (in earth centered earth fixed frame) = {np.round((self.position[-1]-self.position[0])/1e3,2)} km\n"
            out_str+=f"velocity vector at burnout (in earth centered earth fixed frame) = {np.round(self.velocity[-1],2)} m/s\n"
            out_str+=f"acceleration vector at burnout (in earth centered earth fixed frame) = {np.round(self.acceleration[-1],2)} m/s^2\n"
            out_str+=f"direction vector at burnout (in earth centered earth fixed frame) = {np.round(self.orientation[-1],4)}\n"
            out_str+=f"downrange position from launch site (in earth centered inertial frame) = {np.round((self.position_inertial_frame[-1]-Rz@self.position_inertial_frame[0])/1e3,2)} km\n"
            out_str+=f"position of launch site at burnout (in earth centered inertial frame) = {np.round((Rz@self.position_inertial_frame[0])/1e3,2)} km\n"
            out_str+=f"velocity vector at burnout (in earth centered inertial frame) = {np.round(self.velocity_inertial_frame[-1],2)} m/s\n"
            out_str+=f"acceleration vector at burnout (in earth centered inertial frame) = {np.round(self.acceleration_inertial_frame[-1],2)} m/s\n"
            out_str+=f"direction vector at burnout (in earth centered inertial frame) = {np.round(Rz@self.orientation[-1],4)}\n"
            out_str+=f"delta-v achieved = {delta_v_achieved} m/s\n"
            out_str+="\n"
            
            """print(f"downrange position from launch site (in earth centered earth fixed frame) = {np.round((self.position[-1]-self.position[0])/1e3,2)} km")
            print(f"velocity vector at burnout (in earth centered earth fixed frame) = {np.round(self.velocity[-1],2)} m/s")
            print(f"acceleration vector at burnout (in earth centered earth fixed frame) = {np.round(self.acceleration[-1],2)} m/s^2")
            print(f"direction vector at burnout (in earth centered earth fixed frame) = {np.round(self.orientation[-1],4)}")
            print(f"downrange position from launch site (in earth centered inertial frame) = {np.round((self.position_inertial_frame[-1]-Rz@self.position_inertial_frame[0])/1e3,2)} km")
            print(f"position of launch site at burnout (in earth centered inertial frame) = {np.round((Rz@self.position_inertial_frame[0])/1e3,2)} km")
            print(f"velocity vector at burnout (in earth centered inertial frame) = {np.round(self.velocity_inertial_frame[-1],2)} m/s")
            print(f"acceleration vector at burnout (in earth centered inertial frame) = {np.round(self.acceleration_inertial_frame[-1],2)} m/s")
            print(f"direction vector at burnout (in earth centered inertial frame) = {np.round(Rz@self.orientation[-1],4)}")
            print()"""
            #print(f"direction vector at burnout (in earth centered inertial frame) = {np.round(self.orientation_inertial_frame[-1],2)} m/s")
            #print(self.position[-1][2])
        #print(np.linalg.norm(self.velocity[-1]))
        self.burnouts=np.array(self.burnouts)
        self.burnouts_inertial_frame=np.array(self.burnouts_inertial_frame)
        #print(self.acceleration[1],self.velocity[1])
        out_str+="Note that the burnout altitudes are not very accurate as they depend on the actual ascent trajectory which is not yet fixed\n"
        #print("Note that the burnout altitudes are not very accurate as they depend on the actual ascent trajectory which is not yet fixed")

        if post_burnout:
            while self.time[-1]<total_time:
                coriolis=-2*np.cross(self.omega_earth,self.velocity[-1])
                centrifugal=-np.cross(self.omega_earth,np.cross(self.omega_earth,self.position[-1]))
                self.acceleration+=[self.gravity+coriolis+centrifugal]
                self.acceleration_inertial_frame+=[self.acceleration[-1]+2*np.cross(self.omega_earth,self.velocity[-1])+np.cross(self.omega_earth,np.cross(self.omega_earth,self.position[-1]))]
                
                self.velocity+=[self.velocity[-1]+self.acceleration[-1]*dt]
                self.position+=[self.position[-1]+self.velocity[-1]*dt+self.acceleration[-1]*dt*dt/2]#+self.rotation_velocity*dt] 
                self.time+=[self.time[-1]+dt]
                self.orientation+=[self.velocity[-1]/np.linalg.norm(self.velocity[-1])]
                self.gravity=-g*self.position[-1]/np.linalg.norm(self.position[-1])

                self.velocity_inertial_frame+=[self.velocity[-1]+np.cross(self.omega_earth,self.position[-1])]
                
                theta_=w_earth*self.time[-1]
                Rz = np.array([
                                    [ np.cos(theta_), -np.sin(theta_), 0],
                                    [np.sin(theta_), np.cos(theta_), 0],
                                    [0,             0,                1]
                                ])

                #self.position_inertial_frame+=[self.position_inertial_frame[-1]+self.velocity_inertial_frame[-1]*dt+self.acceleration_inertial_frame[-1]*dt*dt/2]
                self.position_inertial_frame+=[Rz@self.position[-1]]

                i+=1
        return out_str
    
    def plotter_inertial_frame(self, position=True, show_earth=True):
        if position:
            self.fig = plt.figure()
            self.ax = self.fig.add_subplot(111, projection='3d')
            x = np.array(self.position_inertial_frame).T[0] / 1e3
            y = np.array(self.position_inertial_frame).T[1] / 1e3
            z = np.array(self.position_inertial_frame).T[2] / 1e3
            self.ax.plot3D(x, y, z)

            x_burnout = self.burnouts_inertial_frame[:, 0] / 1e3  # in km
            y_burnout = self.burnouts_inertial_frame[:, 1] / 1e3
            z_burnout = self.burnouts_inertial_frame[:, 2] / 1e3
            self.ax.scatter(x_burnout, y_burnout, z_burnout, color='red', s=20, label='Burnout Points')


            if show_earth:
            # --- Add Earth surface ---
                R_earth = 6371  # Earth's radius in km
                u = np.linspace(0, 2 * np.pi, 100)
                v = np.linspace(0, np.pi, 100)
                x_earth = R_earth * np.outer(np.cos(u), np.sin(v))
                y_earth = R_earth * np.outer(np.sin(u), np.sin(v))
                z_earth = R_earth * np.outer(np.ones(np.size(u)), np.cos(v))

                self.ax.plot_surface(x_earth, y_earth, z_earth, color='lightblue', alpha=0.5, edgecolor='none')

            # --- Equal aspect ratio ---
            max_range = np.array([x.max()-x.min(), y.max()-y.min(), z.max()-z.min()]).max() / 2.0
            mid_x = (x.max() + x.min()) * 0.5
            mid_y = (y.max() + y.min()) * 0.5
            mid_z = (z.max() + z.min()) * 0.5
            self.ax.set_xlim(mid_x - max_range, mid_x + max_range)
            self.ax.set_ylim(mid_y - max_range, mid_y + max_range)
            self.ax.set_zlim(mid_z - max_range, mid_z + max_range)

            self.ax.set_xlabel('X (km)')
            self.ax.set_ylabel('Y (km)')
            self.ax.set_zlabel('Z (km)')            

            
            plt.show()

    def plotter(self,show_earth=True):
            self.fig = plt.figure()
            self.ax = self.fig.add_subplot(111, projection='3d')
            x = np.array(self.position).T[0] / 1e3
            y = np.array(self.position).T[1] / 1e3
            z = np.array(self.position).T[2] / 1e3
            self.ax.plot3D(x, y, z)

            x_burnout = self.burnouts[:, 0] / 1e3  # in km
            y_burnout = self.burnouts[:, 1] / 1e3
            z_burnout = self.burnouts[:, 2] / 1e3
            self.ax.scatter(x_burnout, y_burnout, z_burnout, color='red', s=20, label='Burnout Points')


            if show_earth:
            # --- Add Earth surface ---
                R_earth = 6371  # Earth's radius in km
                u = np.linspace(0, 2 * np.pi, 100)
                v = np.linspace(0, np.pi, 100)
                x_earth = R_earth * np.outer(np.cos(u), np.sin(v))
                y_earth = R_earth * np.outer(np.sin(u), np.sin(v))
                z_earth = R_earth * np.outer(np.ones(np.size(u)), np.cos(v))

                self.ax.plot_surface(x_earth, y_earth, z_earth, color='lightblue', alpha=0.5, edgecolor='none')

            # --- Equal aspect ratio ---
            max_range = np.array([x.max()-x.min(), y.max()-y.min(), z.max()-z.min()]).max() / 2.0
            mid_x = (x.max() + x.min()) * 0.5
            mid_y = (y.max() + y.min()) * 0.5
            mid_z = (z.max() + z.min()) * 0.5
            self.ax.set_xlim(mid_x - max_range, mid_x + max_range)
            self.ax.set_ylim(mid_y - max_range, mid_y + max_range)
            self.ax.set_zlim(mid_z - max_range, mid_z + max_range)

            self.ax.set_xlabel('X (km)')
            self.ax.set_ylabel('Y (km)')
            self.ax.set_zlabel('Z (km)')            

            output_file = f"output files/{self.data["rocket name"]} files/3d trajectory.png"
            self.fig.savefig(output_file)
            plt.show()

    def plotter_both(self, position=True, show_earth=True):

        if position:
            self.fig = plt.figure()
            self.ax = self.fig.add_subplot(111, projection='3d')
            x = np.array(self.position_inertial_frame).T[0] / 1e3
            y = np.array(self.position_inertial_frame).T[1] / 1e3
            z = np.array(self.position_inertial_frame).T[2] / 1e3
            self.ax.plot3D(x, y, z, label='Inertial')

            x_burnout = self.burnouts_inertial_frame[:, 0] / 1e3  # in km
            y_burnout = self.burnouts_inertial_frame[:, 1] / 1e3
            z_burnout = self.burnouts_inertial_frame[:, 2] / 1e3
            self.ax.scatter(x_burnout, y_burnout, z_burnout, color='red', s=20, label='Burnout Points')


            if show_earth:
            # --- Add Earth surface ---
                R_earth = 6371  # Earth's radius in km
                u = np.linspace(0, 2 * np.pi, 100)
                v = np.linspace(0, np.pi, 100)
                x_earth = R_earth * np.outer(np.cos(u), np.sin(v))
                y_earth = R_earth * np.outer(np.sin(u), np.sin(v))
                z_earth = R_earth * np.outer(np.ones(np.size(u)), np.cos(v))

                self.ax.plot_surface(x_earth, y_earth, z_earth, color='lightblue', alpha=0.5, edgecolor='none')

            # --- Equal aspect ratio ---
            max_range = np.array([x.max()-x.min(), y.max()-y.min(), z.max()-z.min()]).max() / 2.0
            mid_x = (x.max() + x.min()) * 0.5
            mid_y = (y.max() + y.min()) * 0.5
            mid_z = (z.max() + z.min()) * 0.5
            self.ax.set_xlim(mid_x - max_range, mid_x + max_range)
            self.ax.set_ylim(mid_y - max_range, mid_y + max_range)
            self.ax.set_zlim(mid_z - max_range, mid_z + max_range)

            self.ax.set_xlabel('X (km)')
            self.ax.set_ylabel('Y (km)')
            self.ax.set_zlabel('Z (km)')            

            
            x = np.array(self.position).T[0] / 1e3
            y = np.array(self.position).T[1] / 1e3
            z = np.array(self.position).T[2] / 1e3
            self.ax.plot3D(x, y, z,label='Non-Inertial')

            x_burnout = self.burnouts[:, 0] / 1e3  # in km
            y_burnout = self.burnouts[:, 1] / 1e3
            z_burnout = self.burnouts[:, 2] / 1e3
            self.ax.scatter(x_burnout, y_burnout, z_burnout, color='red', s=20, label='Burnout Points')
            self.ax.legend(loc="upper right")

            
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

    def plot_altitudes_inertial_frame(self,show_plots=False):
        plt.figure()
        self.altitudes=[(np.linalg.norm(self.position_inertial_frame[i])-R_earth)/1e3 for i in range(len(self.position_inertial_frame))]
        self.altitudes_burnout_inertial_frame=[(np.linalg.norm(self.burnouts_inertial_frame[i])-R_earth)/1e3 for i in range(len(self.burnouts_inertial_frame))]
        plt.plot(self.time,self.altitudes)
        plt.scatter(self.burnout_times,self.altitudes_burnout_inertial_frame)
        plt.title("Altitude vs Time")
        output_file = f"output files/{self.data["rocket name"]} files/altitude vs time in inertial frame.png"
        plt.savefig(output_file)
        if show_plots:
            plt.show()
        plt.figure()
        self.velocities=[(np.linalg.norm(self.velocity_inertial_frame[i])) for i in range(len(self.velocity_inertial_frame))]
        plt.plot(self.time,self.velocities)
        plt.title("Velocity vs Time")
        output_file = f"output files/{self.data["rocket name"]} files/velocity vs time in inertial frame.png"
        plt.savefig(output_file)
        if show_plots:
            plt.show()

    def plot_altitudes(self,show_plots=False):
        plt.figure()
        self.altitudes=[(np.linalg.norm(self.position[i])-R_earth)/1e3 for i in range(len(self.position))]
        plt.plot(self.time,self.altitudes)
        plt.title("Altitude vs Time")
        if show_plots:
            plt.show()
        plt.figure()
        self.velocities=[(np.linalg.norm(self.velocity[i])) for i in range(len(self.velocity))]
        plt.plot(self.time,self.velocities)
        plt.title("Velocity vs Time")
        if show_plots:
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

    
def optimize_betas_2(data,t1=20,t2=100,n_ts=20,beta_max=0.2,n_betas=20):
        print("Optimising vector end time")
        plt.figure()
        plt.title("Altitude vs Time")
        min_loss=1e100
        ans1=0
        ans2=0
        count=1
        for beta in np.linspace(0.05,beta_max,n_betas):
            for t in np.linspace(t1,t2,n_ts):
                traj=trajectory2(copy.deepcopy(data))  # reset the simulation
                traj.model(vector_end_time=t,beta_max=beta,post_burnout=False)
                traj.altitudes=[(np.linalg.norm(traj.position_inertial_frame[i])-R_earth)/1e3 for i in range(len(traj.position_inertial_frame))]
                traj.altitudes_burnout_inertial_frame=[(np.linalg.norm(traj.burnouts_inertial_frame[i])-R_earth)/1e3 for i in range(len(traj.burnouts_inertial_frame))]
                
                final_altitude=np.linalg.norm(traj.position[-1]) - R_earth
                current_loss1=np.abs(final_altitude - traj.data["target orbit"]*1e3)/1e3
                current_loss2=np.abs(np.dot(traj.velocity_inertial_frame[-1],traj.position_inertial_frame[-1]))/np.linalg.norm(traj.velocity_inertial_frame[-1])/np.linalg.norm(traj.position_inertial_frame[-1])
                current_loss=0*current_loss1+100*current_loss2
                if np.any(traj.altitudes<np.zeros_like(traj.altitudes)):
                    current_loss+=99999
                #print(f"Current Loss = {current_loss1} * {current_loss2} = {current_loss}")
                if current_loss<min_loss:
                    ans1=t
                    ans2=beta
                    min_loss=current_loss
                    plt.clf()
                    plt.plot(traj.time,traj.altitudes)
                    plt.scatter(traj.burnout_times,traj.altitudes_burnout_inertial_frame)
                print(np.round(count/n_betas/n_ts*100,2),"% completed")
                count+=1
        print(f"Vector end time = {ans1}, Beta = {ans2}, min loss= {min_loss}")
        plt.show() 
        return ans1,ans2