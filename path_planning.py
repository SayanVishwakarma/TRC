from __init__ import *
class path_planning:
    def __init__(self,thrust_by_weight,targt_orbit,isp):
        self.position=[np.array([R_earth*np.cos(latitude),0,R_earth*np.sin(latitude)])]
        self.omega_earth=np.array([0,0,w_earth])
        self.velocity=[np.cross(np.array([0,0,w_earth]),self.position[-1])]
        self.acceleration=[np.array([0,0,0])]
        self.aoa=[0]
        self.gravity=-g*np.array([np.cos(latitude),0,np.sin(latitude)])
        self.thrust_position=np.array([0,0,-30])
        self.time=[0.0]
        self.betas=[0.0]
        self.mass=4e6
        self.mass_initial=4e6
        self.inertia_matrix=self.rocket_inertia_tensor(self.mass,100,8)
        self.thrust_by_weight=thrust_by_weight
        self.q=np.array([np.cos((np.pi/2-latitude)/2),0,np.sin((np.pi/2-latitude)/2),0])
        self.orientation=[self.quaternion_to_dcm(self.q)@[0,0,1]]
        self.target_orbit=targt_orbit
        self.isp=isp
        
        self.gravity_delta_v=0
        self.thrust_delta_v=0
        self.steering_loss=0
        self.delta_v=0

    def rocket_inertia_tensor(self,mass, radius, length):
        I_xx = (1/12) * mass * (3 * radius**2 + length**2)
        I_zz = 0.5 * mass * radius**2
        return np.diag([I_xx, I_xx, I_zz])
    
    def normalize_quaternion(self,q):
        return q / np.linalg.norm(q)

    def quaternion_derivative(self,q, omega):
        p, q_b, r = omega
        omega_mat = np.array([
            [ 0,   -p,   -q_b, -r],
            [ p,    0,    r,  -q_b],
            [ q_b, -r,    0,   p],
            [ r,   q_b,  -p,   0]
        ])
        return 0.5 * omega_mat @ q

    def quaternion_to_dcm(self,q):
        q0, q1, q2, q3 = q
        return np.array([
            [1 - 2*(q2**2 + q3**2),     2*(q1*q2 - q0*q3),     2*(q1*q3 + q0*q2)],
            [2*(q1*q2 + q0*q3),         1 - 2*(q1**2 + q3**2), 2*(q2*q3 - q0*q1)],
            [2*(q1*q3 - q0*q2),         2*(q2*q3 + q0*q1),     1 - 2*(q1**2 + q2**2)]
        ])
    
    def get_component(self,a,b):
        return np.dot(a,b)/np.linalg.norm(b)


    def model(self,post_burnout=False,dt=dt,vector_start_time=10,vector_end_time=120,beta_max=0.2,total_time=total_time,optimising_controls=False):
        
        self.thrust_magnitude=self.thrust_by_weight*self.mass*g

        while (np.linalg.norm(self.position[-1])-R_earth-self.target_orbit)<0:
            if self.thrust_magnitude/self.mass/g>4:
                self.thrust_by_weight=4
                self.thrust_magnitude=self.thrust_by_weight*self.mass*g
            
            if self.time[-1]>=vector_start_time and self.time[-1]<=vector_end_time:
                theta = beta_max
            else:
                theta=0
            Rx = np.array([ [1      , 0            , 0],
                            [0      , np.cos(theta), -np.sin(theta)],
                            [0      , np.sin(theta),  np.cos(theta)]
                        ])
            
            thrust_body_frame=Rx@np.array([0,0,self.thrust_magnitude])
            torque_body_frame=np.cross(self.thrust_position,thrust_body_frame)
            omega_body_frame=np.linalg.solve(self.inertia_matrix,torque_body_frame)
            
            R_b_to_i = self.quaternion_to_dcm(self.q)
            thrust_eci_frame = R_b_to_i @ thrust_body_frame

            gravity_eci_frame=-g*self.normalize_quaternion(self.position[-1])
        
            self.gravity_delta_v+=self.get_component(gravity_eci_frame,self.velocity[-1])*dt
            self.thrust_delta_v+=self.get_component(thrust_eci_frame/self.mass,self.velocity[-1])*dt
            self.steering_loss+=(np.linalg.norm(thrust_eci_frame/self.mass)**2-self.get_component(thrust_eci_frame/self.mass,self.velocity[-1])**2)**0.5*dt
            
            self.acceleration+=[thrust_eci_frame/self.mass+gravity_eci_frame]
            self.velocity+=[self.velocity[-1]+self.acceleration[-1]*dt]
            self.position+=[self.position[-1]+self.velocity[-1]*dt+0.5*self.acceleration[-1]*dt*dt]
            self.orientation+=[self.quaternion_to_dcm(self.q)@[0,0,1]]

            dqdt = self.quaternion_derivative(self.q, omega_body_frame)
            self.q += dqdt * dt
            self.q = self.normalize_quaternion(self.q)
            self.mass-=self.thrust_magnitude/self.isp/g*dt    
        
            self.time+=[self.time[-1]+dt]
            
        self.delta_v=self.isp*g*np.log(self.mass_initial/self.mass)    
            
        
        
        
        

        if post_burnout:
            while self.time[-1]<total_time:

                gravity_eci_frame=-g*self.normalize_quaternion(self.position[-1])
                
                self.acceleration+=[gravity_eci_frame]
                self.velocity+=[self.velocity[-1]+self.acceleration[-1]*dt]
                self.position+=[self.position[-1]+self.velocity[-1]*dt+0.5*self.acceleration[-1]*dt*dt]

                self.time+=[self.time[-1]+dt]

        angle=(np.dot(self.position[-1],self.velocity[-1])/np.linalg.norm(self.position[-1])/np.linalg.norm(self.velocity[-1]))/np.pi*180
        orbital_velocity_difference=((np.linalg.norm(self.velocity[-1])-(g*R_earth**2/np.linalg.norm(self.position[-1]))**0.5))

        if optimising_controls:
            return angle,orbital_velocity_difference

        print(f"GRAVITY DELTA V : {self.gravity_delta_v}, THRUST DELTA V : {self.thrust_delta_v}, STEERING LOSS : {self.steering_loss}")
        return self.delta_v
    
    def plot_altitudes(self,show_velocity_plots=False):
        plt.figure()
        self.altitudes=[(np.linalg.norm(self.position[i])-R_earth)/1e3 for i in range(len(self.position))]
        plt.plot(self.time,self.altitudes)
        #self.altitudes_burnout=[(np.linalg.norm(self.burnouts[i])-R_earth)/1e3 for i in range(len(self.burnouts))]
        plt.plot(self.time,self.altitudes)
        #plt.scatter(self.burnout_times,self.altitudes_burnout)
        plt.title("Altitude vs Time")
        plt.show()
        if show_velocity_plots:
            plt.figure()
            self.velocities=[(np.linalg.norm(self.velocity[i])) for i in range(len(self.velocity))]
            plt.plot(self.time,self.velocities)
            plt.title("Velocity vs Time")
            plt.show()

    def plotter(self,show_earth=True):
        self.fig = plt.figure()
        self.ax = self.fig.add_subplot(111, projection='3d')
        x = np.array(self.position).T[0] / 1e3
        y = np.array(self.position).T[1] / 1e3
        z = np.array(self.position).T[2] / 1e3
        self.ax.plot3D(x, y, z)


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

    
    def sweep(self,start_time_array,end_time_array,beta_array):
        thrust_by_weight=self.thrust_by_weight*1
        isp=self.isp
        target_orbit=self.target_orbit
        angle_min=1e10
        orbital_velocity_difference_min=1e10
        prod_min=1e10
        start_time_min=0
        end_time_min=0
        beta_min=0
        num=0
        loss=[]
        for i in start_time_array:
            for j in end_time_array:
                for k in beta_array:
                    x,y=self.model(vector_start_time=i,vector_end_time=j,beta_max=k,optimising_controls=True)
                    if prod_min>x*y:#x<angle_min and y<orbital_velocity_difference_min:
                        start_time_min=i*1
                        end_time_min=j*1
                        beta_min=k*1
                        angle_min=x
                        orbital_velocity_difference_min=y
                        prod_min=angle_min*orbital_velocity_difference_min
                        print(x,y)
                    loss+=[x*y]
                    self.__init__(thrust_by_weight=thrust_by_weight,isp=isp,targt_orbit=target_orbit)
                    num+=1
                    print(f"{(int)(num/len(start_time_array)/len(end_time_array)/len(beta_array)*100)}% completed")
        plt.figure()
        plt.plot(loss)
        plt.show()
        return start_time_min,end_time_min,beta_min

obj=path_planning(1.3,500e3,330)
#print(obj.model(beta_max=0,dt=1))
#obj.plot_altitudes()
start_time_array=np.linspace(5,50,10)
end_time_array=np.linspace(50,250,10)
beta_array=np.linspace(-0.2,0.2,10)
mins=obj.sweep(start_time_array=start_time_array,end_time_array=end_time_array,beta_array=beta_array)
obj.model(vector_start_time=mins[0],vector_end_time=mins[1],beta_max=mins[2])
print(mins)
obj.plot_altitudes()
obj.plotter()