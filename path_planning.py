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
        self.mass=16e6
        self.mass_initial=16e6
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
        self.orbital_velocity_difference=0

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


    def model(self,vector_start_time=10,vector_end_time=120,beta_max=0.2,thrust_cutoff=60*8,dt=dt,total_time=total_time,post_burnout=False,display_breakdown=False):
        #self.__init__(thrust_by_weight=self.thrust_by_weight,isp=self.isp,targt_orbit=self.target_orbit)
        self.thrust_magnitude=self.thrust_by_weight*self.mass*g

        while (np.linalg.norm(self.position[-1])-R_earth-self.target_orbit)<0 and np.linalg.norm(self.position[-1])-R_earth>=0:
            if self.thrust_magnitude/self.mass/g>G_force_limit:
                #self.thrust_by_weight=G_force_limit
                self.thrust_magnitude=G_force_limit*self.mass*g

            if self.time[-1]>thrust_cutoff:
                self.thrust_magnitude=0
            
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

            gravity_eci_frame=-g*(R_earth**2/np.linalg.norm(self.position[-1])**2)*self.normalize_quaternion(self.position[-1])
        
            self.gravity_delta_v+=self.get_component(gravity_eci_frame,self.velocity[-1])*dt
            self.thrust_delta_v+=self.get_component(thrust_eci_frame/self.mass,self.velocity[-1])*dt
            self.steering_loss+=(np.linalg.norm(thrust_eci_frame/self.mass)**2-self.get_component(thrust_eci_frame/self.mass,self.velocity[-1])**2)**0.5*dt
            self.drag_delta_v=50
            
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
        
        #self.angle=(np.dot(self.position[-1],self.velocity[-1])/np.linalg.norm(self.position[-1])/np.linalg.norm(self.velocity[-1]))/np.pi*180
        self.angle=90-np.arccos(np.dot(self.position[-1],self.velocity[-1])/np.linalg.norm(self.position[-1])/np.linalg.norm(self.velocity[-1]))/np.pi*180
        orbital_velocity_difference=((np.linalg.norm(self.velocity[-1])-(g*R_earth**2/np.linalg.norm(self.position[-1]))**0.5))
        self.orbital_velocity_difference=orbital_velocity_difference
        
        self.total_delta_v=self.delta_v+np.abs(self.orbital_velocity_difference)+self.drag_delta_v   
        
        if np.linalg.norm(self.position[-1])-R_earth<0:
            self.angle=1e3
            self.orbital_velocity_difference=1e6
            
        if post_burnout:
            dt=dt/100
            #print((g*R_earth**2/np.linalg.norm(self.position[-1]))**0.5)
            #print(self.angle,self.orbital_velocity_difference)
            #print(np.arccos(np.dot(self.position[-1],self.velocity[-1])/np.linalg.norm(self.velocity[-1])/np.linalg.norm(self.position[-1]))*180/np.pi)
            #print(np.arccos((self.velocity[-1])/np.linalg.norm(self.velocity[-1]))/pi*180)
            #print(np.linalg.norm(self.velocity[-1]))
            #self.velocity[-1]*=((g*R_earth**2/np.linalg.norm(self.position[-1]))**0.5+600)/np.linalg.norm(self.velocity[-1])
            self.velocity[-1]=np.cross(np.cross(self.position[-1],self.velocity[-1]),self.position[-1])*((g*R_earth**2/np.linalg.norm(self.position[-1]))**0.5)/np.linalg.norm(self.position[-1])**2/np.linalg.norm(self.velocity[-1])
            #print(np.linalg.norm(np.cross(np.cross(self.position[-1],self.velocity[-1]),self.position[-1])/np.linalg.norm(self.position[-1])**2/np.linalg.norm(self.velocity[-1])))
            #print(np.arccos(np.dot(self.position[-1],self.velocity[-1])/np.linalg.norm(self.velocity[-1])/np.linalg.norm(self.position[-1]))*180/pi)
            #print(np.arccos((self.velocity[-1])/np.linalg.norm(self.velocity[-1]))/pi*180)            
            #print(np.linalg.norm(self.velocity[-1]))
            #rint(np.dot(self.position[-1],self.velocity[-1])/np.linalg.norm(self.velocity[-1])/np.linalg.norm(self.position[-1])*180/pi)
            while self.time[-1]<total_time:

                gravity_eci_frame=-g*(R_earth**2/np.linalg.norm(self.position[-1])**2)*self.normalize_quaternion(self.position[-1])
                
                self.acceleration+=[gravity_eci_frame]
                self.velocity+=[self.velocity[-1]+self.acceleration[-1]*dt]
                self.position+=[self.position[-1]+self.velocity[-1]*dt+0.5*self.acceleration[-1]*dt*dt]

                self.time+=[self.time[-1]+dt]

        
        if display_breakdown:
            print(f"TOTAL DELTA V : {self.total_delta_v}, GRAVITY DELTA V : {self.gravity_delta_v}, STEERING LOSS : {self.steering_loss}, DRAG DELTA V : {self.drag_delta_v}, ORBITAL VELOCITY DIFFERENCE : {orbital_velocity_difference}, INSERTION ANGLE : {self.angle}")
        return self.delta_v+np.abs(self.orbital_velocity_difference)

    def model_obsolete(self,post_burnout=False,dt=dt,vector_start_time=10,vector_end_time=120,beta_max=0.2,total_time=total_time,optimising_controls=False,thrust_cutoff=60*8):
        #self.__init__(thrust_by_weight=self.thrust_by_weight,isp=self.isp,targt_orbit=self.target_orbit)
        self.thrust_magnitude=self.thrust_by_weight*self.mass*g

        while (np.linalg.norm(self.position[-1])-R_earth-self.target_orbit)<0:
            if self.thrust_magnitude/self.mass/g>4:
                self.thrust_by_weight=4
                self.thrust_magnitude=self.thrust_by_weight*self.mass*g

            if self.time[-1]>thrust_cutoff:
                self.thrust_magnitude=0
            
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

            gravity_eci_frame=-g*(R_earth**2/np.linalg.norm(self.position[-1])**2)*self.normalize_quaternion(self.position[-1])
        
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
            if False:#self.time[-1]>total_time:
                if optimising_controls:
                    return 1e3,1e6
                else:
                    return 1e6
            
        self.delta_v=self.isp*g*np.log(self.mass_initial/self.mass)    
            
        
        
        
        

        if post_burnout:
            print(np.dot(self.position[-1],self.velocity[-1])/np.linalg.norm(self.velocity[-1])/np.linalg.norm(self.position[-1])*180/pi)
            print(np.arccos((self.velocity[-1])/np.linalg.norm(self.velocity[-1]))/pi*180)
            print(np.linalg.norm(self.velocity[-1]))
            #self.velocity[-1]*=((g*R_earth**2/np.linalg.norm(self.position[-1]))**0.5+600)/np.linalg.norm(self.velocity[-1])
            self.velocity[-1]=np.cross(np.cross(self.position[-1],self.velocity[-1]),self.position[-1])*((g*R_earth**2/np.linalg.norm(self.position[-1]))**0.5)/np.linalg.norm(self.position[-1])**2/np.linalg.norm(self.velocity[-1])
            print(np.dot(self.position[-1],self.velocity[-1])/np.linalg.norm(self.velocity[-1])/np.linalg.norm(self.position[-1])*180/pi)
            print(np.arccos((self.velocity[-1])/np.linalg.norm(self.velocity[-1]))/pi*180)            
            print(np.linalg.norm(self.velocity[-1]))
            #rint(np.dot(self.position[-1],self.velocity[-1])/np.linalg.norm(self.velocity[-1])/np.linalg.norm(self.position[-1])*180/pi)
            while self.time[-1]<total_time:

                gravity_eci_frame=-g*(R_earth**2/np.linalg.norm(self.position[-1])**2)*self.normalize_quaternion(self.position[-1])
                
                self.acceleration+=[gravity_eci_frame]
                self.velocity+=[self.velocity[-1]+self.acceleration[-1]*dt]
                self.position+=[self.position[-1]+self.velocity[-1]*dt+0.5*self.acceleration[-1]*dt*dt]

                self.time+=[self.time[-1]+dt]

        angle=(np.dot(self.position[-1],self.velocity[-1])/np.linalg.norm(self.position[-1])/np.linalg.norm(self.velocity[-1]))/np.pi*180
        orbital_velocity_difference=((np.linalg.norm(self.velocity[-1])-(g*(R_earth**2/np.linalg.norm(self.position[-1])**2)*R_earth**2/np.linalg.norm(self.position[-1]))**0.5))
        self.orbital_velocity_difference=orbital_velocity_difference
        if optimising_controls:
            return angle,orbital_velocity_difference

        print(f"GRAVITY DELTA V : {self.gravity_delta_v}, THRUST DELTA V : {self.thrust_delta_v}, STEERING LOSS : {self.steering_loss}")
        return self.delta_v+np.abs(self.orbital_velocity_difference)
    
    
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

    def sweep(self,start_time_array,end_time_array,beta_array,thrust_cutoff_array):
        thrust_by_weight=self.thrust_by_weight*1
        isp=self.isp
        target_orbit=self.target_orbit
        angle_min=1e10
        orbital_velocity_difference_min=1e10
        delta_v_min=1e10
        start_time_min=0
        end_time_min=0
        beta_min=0
        thrust_cutoff_min=0
        with tqdm.tqdm(total=len(start_time_array)*len(end_time_array)*len(beta_array)*len(thrust_cutoff_array)) as pbar:
            for i in start_time_array:
                for j in end_time_array:
                    for k in beta_array:
                        for l in thrust_cutoff_array:
                            self.model(vector_start_time=i,vector_end_time=j,beta_max=k,thrust_cutoff=l)#,optimising_controls=True)
                            x=np.abs(self.angle)
                            y=self.orbital_velocity_difference
                            if self.total_delta_v<delta_v_min and x<3 and y<0:#(np.abs(x)**2*np.abs(y))<prod_min:
                                start_time_min=i*1
                                end_time_min=j*1
                                beta_min=k*1
                                thrust_cutoff_min=l*1
                                angle_min=np.abs(x*1)
                                orbital_velocity_difference_min=np.abs(y*1)
                                delta_v_min=self.total_delta_v*1
                            self.__init__(thrust_by_weight=thrust_by_weight,isp=isp,targt_orbit=target_orbit)
                            pbar.update(1)
                            '''num+=1
                            percent_complete=(num/len(start_time_array)/len(end_time_array)/len(beta_array)/len(thrust_cutoff_array))*100
                            if (int)(percent_complete)/report_percentage/lol==1 or num==10:
                                print(f"{(int)(percent_complete)}% completed. ETA {int((time.time()-startingtime)*(100/percent_complete-1)/60)} minutes {int((time.time()-startingtime)*(100/percent_complete-1)%60)} seconds. {(int)(num/(time.time()-startingtime))} iters/s")
                                lol+=1'''
        return start_time_min,end_time_min,beta_min,thrust_cutoff_min,angle_min,orbital_velocity_difference_min,delta_v_min
    
    def epoch(self,start_time=25,end_time=200,beta=0.1,thrust_cutoff=5*60,n=5,delta_start_time=20,delta_end_time=150,delta_beta=0.09,delta_thrust_cutoff=3*60,n_epochs=10,dynamic_n=False):
        n_epoch=n_epochs
        n1=n*1
        #with tqdm.tqdm(total=n_epochs) as pbar:
        for i in range(0,n_epoch):
                if dynamic_n and n_epoch>1:
                    n1=n+(3-n)/(n_epoch-1)*i
                n1=int(n1)
                #start_time_array=np.linspace(start_time-delta_start_time,start_time+delta_start_time,n1)
                start_time_array=[5]
                end_time_array=np.linspace(end_time-delta_end_time,end_time+delta_end_time,n1)
                beta_array=np.linspace(beta-delta_beta,beta+delta_beta,n1)
                thrust_cutoff_array=np.linspace(thrust_cutoff-delta_thrust_cutoff,thrust_cutoff+delta_thrust_cutoff,n1)
                print(f"Commencing epoch {i+1}/{n_epoch}")
                time.sleep(1)
                start_time,end_time,beta,thrust_cutoff,angle_min,orbital_v_diff_min,delta_v_min=self.sweep(start_time_array=start_time_array,end_time_array=end_time_array,beta_array=beta_array,thrust_cutoff_array=thrust_cutoff_array)
                a=0.95
                delta_start_time=a*delta_start_time/n
                delta_end_time=a*delta_end_time/n
                delta_beta=a*delta_beta/n
                delta_thrust_cutoff=a*delta_thrust_cutoff/n 
                #print(start_time,delta_start_time,end_time,delta_end_time,beta,delta_beta,thrust_cutoff,delta_thrust_cutoff)   
                print(f"Current Angle Difference= {angle_min}, Current Orbital Velocity Difference = {orbital_v_diff_min}, Current Delta V = {delta_v_min}")   
                #pbar.update(1)
        return start_time,end_time,beta,thrust_cutoff

    def path_planner(self,n=3,n_epochs=8,dynamic_n=False,check_history=True):
        with open('data files/delta v tables.txt','r') as f:
            lines=f.readlines()
            min_d_v=1e10
            if check_history:
                for line in lines[1:]:
                    if float(line.split()[0])==self.thrust_by_weight and float(line.split()[1])==self.target_orbit and float(line.split()[2])==self.isp:
                        if float(line.split()[7])<min_d_v:
                            self.mins=[float(i) for i in line.split()[3:]]
                            mins=self.mins
                            min_d_v=float(line.split()[7])
                            print("Optimal Parameters Found In Database With Delta V = ",min_d_v,".\nOptimal Parameters= ",self.mins,"\n")
                if min_d_v<1e10:
                    self.model(vector_start_time=mins[0],vector_end_time=mins[1],beta_max=mins[2],thrust_cutoff=mins[3],post_burnout=True,display_breakdown=True)
                    return
        mins=self.epoch(n=n,n_epochs=n_epochs,dynamic_n=dynamic_n)
        self.mins=mins
        if np.all(self.mins)!=0:
            print("Optimal Parameters = ",mins,"\n")
            self.model(vector_start_time=mins[0],vector_end_time=mins[1],beta_max=mins[2],thrust_cutoff=mins[3],post_burnout=False,display_breakdown=True)
            with open('data files/delta v tables.txt','a') as f:
                f.write(f"{self.thrust_by_weight} {self.target_orbit} {self.isp} {' '.join([str(i) for i in mins])} {self.total_delta_v}\n")
            #self.plot_altitudes()

    def plot_variations(self,start_time_array,end_time_array,beta_array,thrust_cutoff_array):
        plt.figure()
        thrust_by_weight=self.thrust_by_weight*1
        isp=self.isp
        target_orbit=self.target_orbit
        delta_vs=[]
        with tqdm.tqdm(total=len(start_time_array)*len(end_time_array)*len(beta_array)*len(thrust_cutoff_array)) as pbar:
            for i in start_time_array:
                for j in end_time_array:
                    for k in beta_array:
                        for l in thrust_cutoff_array:
                            self.model(vector_start_time=i,vector_end_time=j,beta_max=k,thrust_cutoff=l)#,optimising_controls=True)
                            x=np.abs(self.angle)
                            y=self.orbital_velocity_difference
                            if x<3 and y<0:#(np.abs(x)**2*np.abs(y))<prod_min:
                                delta_vs+=[self.total_delta_v*1]
                            else:
                                delta_vs+=[0]
                            self.__init__(thrust_by_weight=thrust_by_weight,isp=isp,targt_orbit=target_orbit)
                            pbar.update(1)
        for arr in [start_time_array,end_time_array,beta_array,thrust_cutoff_array]:
            if len(arr)>1:
                plt.plot(arr,delta_vs)
                for name, value in locals().items():
                    if value is arr:
                        output_file = f"data files/delta_v variation plots/{self.thrust_by_weight,self.target_orbit,self.isp,name}.png"
                        break
                break
        plt.savefig(output_file)
        print(f"Plot saved to {output_file}")
        plt.show()

class step_one_boostback_path_planning:

    def __init__(self,data):
        self.position=[data['initial position']]
        self.omega_earth=[np.array([0,0,w_earth])]
        self.velocity=[data['initial velocity']]
        #self.velocity=[np.array([0.001,0,0])]
        self.acceleration=[data['initial acceleration']]
        #self.aoa=[0]
        self.gravity=-g*np.array(self.position)/np.linalg.norm(self.position)
        self.thrust_position=np.array([0,0,-30])
        self.time=[0.0]
        self.data=data
        self.data_copy=copy.deepcopy(data)
        self.betas=[0.0]
        self.step_mass=data['step masses']
        self.step_propellant_mass=data['propellant mass']
        self.inertia_matrix=self.rocket_inertia_tensor(self.step_mass,70,8)
        #self.thrust_by_weight=data['thrust by weight']
        #self.target_orbit=500e3
        #self.number_of_stages=len(self.stage_masses)
        self.burnouts=[]
        self.burnout_times=[]
        self.q=self.data["q"]
        self.orientation=[self.quaternion_to_dcm(self.q)@[0,0,1]]
        self.single_engine_thrust=data["single engine thrust"]
        self.aoa=[np.arccos(np.dot((self.velocity[-1]-np.cross(np.array([0,0,w_earth]),self.position[-1])),self.quaternion_to_dcm(self.q)@[0,0,1])/np.linalg.norm((self.velocity[-1]-np.cross(np.array([0,0,w_earth]),self.position[-1]))))]

        #self.thrust_by_weight=self.data['thrust by weight']
        #self.step_thrust_magnitude=self.thrust_by_weight*self.step_mass*g
        self.step_thrust_magnitude=3*self.single_engine_thrust
        self.thrust_magnitude=self.step_thrust_magnitude

        self.velocity_reduction_boostback_duration=0
        self.orientation_flipped=False

        self.gravity_delta_v=0
        self.thrust_delta_v=0
        self.steering_loss=0
        self.drag_delta_v=0
        self.total_delta_v=0

        self.total_rotation=np.array([0.0,0.0,0.0])

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
    
    def quat_mul(self,a,b):
        a0,a1,a2,a3 = a
        b0,b1,b2,b3 = b
        return np.array([
            a0*b0 - a1*b1 - a2*b2 - a3*b3,
            a0*b1 + a1*b0 + a2*b3 - a3*b2,
            a0*b2 - a1*b3 + a2*b0 + a3*b1,
            a0*b3 + a1*b2 - a2*b1 + a3*b0
        ])

    def dcm_to_quat(R):
        """Convert 3x3 DCM to quaternion [w,x,y,z]."""
        R = np.array(R, dtype=float)
        tr = np.trace(R)
        if tr > 0:
            S = np.sqrt(tr+1.0) * 2
            w = 0.25 * S
            x = (R[2,1] - R[1,2]) / S
            y = (R[0,2] - R[2,0]) / S
            z = (R[1,0] - R[0,1]) / S
        elif (R[0,0] > R[1,1]) and (R[0,0] > R[2,2]):
            S = np.sqrt(1.0 + R[0,0] - R[1,1] - R[2,2]) * 2
            w = (R[2,1] - R[1,2]) / S
            x = 0.25 * S
            y = (R[0,1] + R[1,0]) / S
            z = (R[0,2] + R[2,0]) / S
        elif R[1,1] > R[2,2]:
            S = np.sqrt(1.0 + R[1,1] - R[0,0] - R[2,2]) * 2
            w = (R[0,2] - R[2,0]) / S
            x = (R[0,1] + R[1,0]) / S
            y = 0.25 * S
            z = (R[1,2] + R[2,1]) / S
        else:
            S = np.sqrt(1.0 + R[2,2] - R[0,0] - R[1,1]) * 2
            w = (R[1,0] - R[0,1]) / S
            x = (R[0,2] + R[2,0]) / S
            y = (R[1,2] + R[2,1]) / S
            z = 0.25 * S
        q = np.array([w, x, y, z])
        return q / np.linalg.norm(q)
    
    def orientation_to_dcm(forward, up=np.array([0,0,1])):
        f = np.array(forward, dtype=float)
        f /= np.linalg.norm(f)
        # if forward is (anti)parallel to up, choose a default right axis
        if np.linalg.norm(np.cross(up, f)) < 1e-8:
            # forward parallel to up -> no yaw info; choose right = +x
            r = np.array([1.0, 0.0, 0.0])
            # if forward is -up (i.e. [0,0,-1]) flip right to keep right-handedness
            if np.dot(f, up) < 0:
                r = -r
        else:
            r = np.cross(up, f)
            r /= np.linalg.norm(r)
        u = np.cross(f, r)
        R = np.column_stack((r, u, f))
        return R

    def normalize(self,v):
        return v / np.linalg.norm(v)
    
    def get_component(self,a,b):
        return np.dot(a,b)/np.linalg.norm(b)

    def model(self,post_burnout=False,dt=dt,vector_start_time=10,vector_end_time=120,beta_max=0.2,beta_max_y=0.01,thrust_cutoff_time=30,total_time=total_time,display_breakdown=False):
        n=1
        out_str=f"\nRocket liftoff from {np.round(self.position[0]/1000,2)} km. Earth rotates along the Z-axis\n"
        while (np.linalg.norm(self.position[-1])-R_earth)>0:# and np.linalg.norm(self.position[-1])-R_earth>=0:
            if self.time[-1]>=vector_start_time:
                #if self.stage_propellant_masses[current_stage]<0:
                    #print("Thrust cutoff due to lack of fuel")
                #print("flip maneuver started at ",self.time[-1])
                self.thrust_magnitude= self.step_thrust_magnitude
                #print(f"Thrust magnitude set to 0 at {self.time[-1]}s and {(np.linalg.norm(self.position[-1])-R_earth)}")
            
            if self.time[-1]>vector_end_time and np.linalg.norm(self.position[-1])-R_earth>60e3:
                self.thrust_magnitude=0
                #print("flip maneuver ended at ",self.time[-1])
                if len(self.burnouts)==0:
                    self.burnouts+=[self.position[-1]]
                    self.burnout_times+=[self.time[-1]]
                #print(f"Thrust magnitude set to 0 at {self.time[-1]}s and {(np.linalg.norm(self.position[-1])-R_earth)}")
            
            if np.linalg.norm(self.position[-1])-R_earth<=60e3:
                if len(self.burnouts)==1:
                    self.burnouts+=[self.position[-1]]
                    self.burnout_times+=[self.time[-1]]
                if self.orientation_flipped==False:
                    #print(f"Flipping orientation at {self.time[-1]}s")
                    #print(f"Original orientation: {self.orientation[-1]}")
                    #print(f"Original quaternion: {self.q}")
                    Rxflip = np.array([ [1      , 0            , 0],
                            [0      , np.cos(pi), -np.sin(pi)],
                            [0      , np.sin(pi),  np.cos(pi)]
                        ])
                    self.orientation[-1]=Rxflip@self.orientation[-1]
                    self.q =self.quat_mul(np.array([0,1,0,0]),self.q) 
                    self.q = self.normalize_quaternion(self.q)
                    self.orientation_flipped=True

                    self.total_rotation+=np.array([pi,0,0])
                    #print(f"New orientation: {self.orientation[-1]}")
                    #print(f"New quaternion: {self.q}")
                #print(f"Thrust cutoff due to altitude {np.linalg.norm(self.position[-1])-R_earth}")
                if self.velocity_reduction_boostback_duration<thrust_cutoff_time:
                    self.thrust_magnitude=self.single_engine_thrust
                    self.velocity_reduction_boostback_duration+=dt
                    #print("velocity reduction started at",self.time[-1])
                else:
                    self.thrust_magnitude=0
                    #print("velocity reduction ended at",self.time[-1])
                    if len(self.burnouts)==2:
                        self.burnouts+=[self.position[-1]]
                        self.burnout_times+=[self.time[-1]]
                    #print(f"Thrust magnitude set to 0 at {self.time[-1]}s and {(np.linalg.norm(self.position[-1])-R_earth)}")  
                #print(f"Thrust magnitude set to 0 at {self.time[-1]}s and {(np.linalg.norm(self.position[-1])-R_earth)}")
            
            if np.linalg.norm(self.position[-1])-R_earth<2e3:
                self.thrust_magnitude=self.single_engine_thrust

            if self.step_propellant_mass<=0:# and self.thrust_magnitude>0:
                #print("Fuel depleted, thrust magnitude set to 0")
                self.thrust_magnitude=0
                self.burnouts+=[self.position[-1]]
                self.burnout_times+=[self.time[-1]]

            if self.time[-1]>=vector_start_time and self.time[-1]<=vector_end_time:
                theta = -beta_max
                theta2 = -beta_max_y
            else:
                theta=0
                theta2=0
            Rx = np.array([ [1      , 0            , 0],
                            [0      , np.cos(theta), -np.sin(theta)],
                            [0      , np.sin(theta),  np.cos(theta)]
                        ])
            
            Ry = np.array([
                            [ np.cos(theta2), 0, np.sin(theta2)],
                            [ 0,             1, 0           ],
                            [-np.sin(theta2), 0, np.cos(theta2)]
                        ])
            
            thrust_body_frame=Ry@Rx@np.array([0,0,self.thrust_magnitude])
            torque_body_frame=np.cross(self.thrust_position,thrust_body_frame)*1e1
            omega_body_frame=np.linalg.solve(self.inertia_matrix,torque_body_frame)

            self.total_rotation+=omega_body_frame*dt
            
            R_b_to_i = self.quaternion_to_dcm(self.q)
            thrust_eci_frame = R_b_to_i @ thrust_body_frame

            gravity_eci_frame=-g*(R_earth**2/np.linalg.norm(self.position[-1])**2)*self.normalize_quaternion(self.position[-1])

            self.gravity_delta_v+=self.get_component(gravity_eci_frame,self.velocity[-1])*dt
            self.thrust_delta_v+=self.get_component(thrust_eci_frame/self.step_mass,self.velocity[-1])*dt
            self.steering_loss+=(np.linalg.norm(thrust_eci_frame/self.step_mass)**2-self.get_component(thrust_eci_frame/self.step_mass,self.velocity[-1])**2)**0.5*dt
            self.drag_delta_v=50
            self.total_delta_v+=self.thrust_magnitude/self.step_mass*dt

            self.acceleration+=[thrust_eci_frame/self.step_mass+gravity_eci_frame]
            self.velocity+=[self.velocity[-1]+self.acceleration[-1]*dt]
            self.position+=[self.position[-1]+self.velocity[-1]*dt+0.5*self.acceleration[-1]*dt*dt]
            self.orientation+=[self.quaternion_to_dcm(self.q)@[0,0,1]]
            self.aoa+=[np.arccos(np.dot((self.velocity[-1]-np.cross(np.array([0,0,w_earth]),self.position[-1])),self.quaternion_to_dcm(self.q)@[0,0,1])/np.linalg.norm((self.velocity[-1]-np.cross(np.array([0,0,w_earth]),self.position[-1]))))]

            dqdt = self.quaternion_derivative(self.q, omega_body_frame)
            self.q += dqdt * dt
            self.q = self.normalize_quaternion(self.q)
            self.step_mass-=self.thrust_magnitude/self.data['isp']/g*dt
            self.step_propellant_mass-=self.thrust_magnitude/self.data['isp']/g*dt
            self.inertia_matrix=self.rocket_inertia_tensor(self.step_mass,70,8)

            self.time+=[self.time[-1]+dt]

        if display_breakdown:
            print(f"GRAVITY DELTA V : {self.gravity_delta_v}, STEERING LOSS : {self.steering_loss}, DRAG DELTA V : {self.drag_delta_v}")

        self.burnouts=np.array(self.burnouts)
        theta_final=w_earth*self.time[-1]
        self.positon_difference=np.linalg.norm(self.position[-1])-R_earth*np.array([np.cos(latitude)*np.cos(theta_final),np.cos(latitude)*np.sin(theta_final),np.sin(latitude)])
        self.absolute_position_difference=np.linalg.norm(self.positon_difference)
        self.final_velocity=np.linalg.norm(self.velocity[-1])

        return out_str

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
            #self.fig.savefig(output_file)
            plt.show()

    def dyn_plotter(self,show_earth=True,animation_speed=1):
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

        dt_interval = dt*1000
        ani = FuncAnimation(self.fig, self.update, frames=len(self.position), init_func=self.init,
                    interval=dt_interval, blit=True)
        

        x_burnout = self.burnouts[:, 0] / 1e3  # in km
        y_burnout = self.burnouts[:, 1] / 1e3
        z_burnout = self.burnouts[:, 2] / 1e3
        self.ax.scatter(x_burnout, y_burnout, z_burnout, color='red', s=20, label='Burnout Points')

        self.ax.set_xlabel('X (km)')
        self.ax.set_ylabel('Y (km)')
        self.ax.set_zlabel('Z (km)')

        plt.show()

    def init(self):
        self.line.set_data([], [])
        self.line.set_3d_properties([])
        self.point.set_data([], [])
        self.point.set_3d_properties([])
        return self.line, self.point

    def update(self,frame):
        x = self.positions[:frame, 0]
        y = self.positions[:frame, 1]
        z = self.positions[:frame, 2]
        self.line.set_data(x, y)
        self.line.set_3d_properties(z)
    
        self.point.set_data(x[-1:], y[-1:])
        self.point.set_3d_properties(z[-1:])
        return self.line, self.point

    def plot_magnitudes(self,altitude=True,velocity=False,aoa=False):
        if altitude:
            plt.figure()
            self.altitudes=[(np.linalg.norm(self.position[i])-R_earth)/1e3 for i in range(len(self.position))]
            plt.plot(self.time,self.altitudes)
            self.altitudes_burnout=[(np.linalg.norm(self.burnouts[i])-R_earth)/1e3 for i in range(len(self.burnouts))]
            plt.plot(self.time,self.altitudes)
            plt.scatter(self.burnout_times,self.altitudes_burnout)
            plt.title("Altitude vs Time")
            plt.show()
        if velocity:
            plt.figure()
            self.velocities=[(np.linalg.norm(self.velocity[i])) for i in range(len(self.velocity))]
            plt.plot(self.time,self.velocities)
            for times in self.burnout_times:
                plt.plot([times,times],[np.min(self.velocities),np.max(self.velocities)])
            plt.title("Velocity vs Time")
            plt.show()
        if aoa:
            plt.figure()
            plt.plot(self.time,np.array(self.aoa)/pi*180)
            for times in self.burnout_times:
                plt.plot([times,times],[np.min(self.aoa)/pi*180,np.max(self.aoa)/pi*180])
            plt.title("AOA (degrees) vs Time")
            plt.show()

    def sweep(self,start_time_array,end_time_array,beta_array,beta_y_array,thrust_cutoff_array):

        angle_min=1e10
        orbital_velocity_difference_min=1e10
        total_delta_v=1e10
        orbit_difference_min=0
        propellant_left=0
        start_time_min=0
        end_time_min=0
        beta_min=0
        thrust_cutoff_min=0

        with tqdm.tqdm(total=len(start_time_array)*len(end_time_array)*len(beta_array)*len(beta_y_array),*len(thrust_cutoff_array)) as pbar:
            for i in start_time_array:
                for j in end_time_array:
                    for k in beta_array:
                        for m in beta_y_array:
                            for l in thrust_cutoff_array:
                                self.model(vector_start_time=i,vector_end_time=j,beta_max=k,beta_max_y=m,thrust_cutoff_time=l)
                                theta_final=w_earth*self.time[-1]
                                positon_difference=np.linalg.norm(self.position[-1])-R_earth*np.array([np.cos(latitude)*np.cos(theta_final),np.cos(latitude)*np.sin(theta_final),np.sin(latitude)])
                                absolute_position_difference=np.linalg.norm(positon_difference)
                                final_velocity=np.linalg.norm(self.velocity[-1])
                                x=absolute_position_difference
                                y=self.orbital_velocity_difference
                                if self.total_delta_v<total_delta_v and x<3 and y<0 and (self.orbit_difference)>0:#(np.abs(x)**2*np.abs(y))<prod_min:
                                    start_time_min=i*1
                                    end_time_min=j*1
                                    beta_min=k*1
                                    thrust_cutoff_min=l*1
                                    angle_min=np.abs(x*1)
                                    orbital_velocity_difference_min=np.abs(y*1)
                                    total_delta_v=self.total_delta_v*1
                                    orbit_difference_min=self.orbit_difference
                                    propellant_left=self.propellant_left
                                self.__init__(self.data_copy)
                                pbar.update(1)
                            
        return start_time_min,end_time_min,beta_min,thrust_cutoff_min,angle_min,orbital_velocity_difference_min,orbit_difference_min,propellant_left,total_delta_v
    
    def epoch(self,start_time=25,end_time=200,beta=0.1,thrust_cutoff=5*60,n=5,delta_start_time=20,delta_end_time=150,delta_beta=0.09,delta_thrust_cutoff=3*60,n_epochs=10,dynamic_n=False):
        n_epoch=n_epochs
        n1=n*1
        #with tqdm.tqdm(total=n_epochs) as pbar:
        for i in range(0,n_epoch):
                if dynamic_n and n_epoch>1:
                    n1=n+(3-n)/(n_epoch-1)*i
                n1=int(n1)
                #start_time_array=np.linspace(start_time-delta_start_time,start_time+delta_start_time,n1)
                start_time_array=[5]
                end_time_array=np.linspace(end_time-delta_end_time,end_time+delta_end_time,n1)
                beta_array=np.linspace(beta-delta_beta,beta+delta_beta,n1)
                thrust_cutoff_array=np.linspace(thrust_cutoff-delta_thrust_cutoff,thrust_cutoff+delta_thrust_cutoff,n1)
                print(f"Commencing epoch {i+1}/{n_epoch}")
                time.sleep(1)
                start_time,end_time,beta,thrust_cutoff,angle_min,orbital_v_diff_min,orbit_difference_min,propellant_left,total_delta_v=self.sweep(start_time_array=start_time_array,end_time_array=end_time_array,beta_array=beta_array,thrust_cutoff_array=thrust_cutoff_array)
                a=0.9
                delta_start_time=a*delta_start_time/n
                delta_end_time=a*delta_end_time/n
                delta_beta=a*delta_beta/n
                delta_thrust_cutoff=a*delta_thrust_cutoff/n 
                #print(start_time,delta_start_time,end_time,delta_end_time,beta,delta_beta,thrust_cutoff,delta_thrust_cutoff)   
                print(f"Total delta v= {total_delta_v}, Current Angle Difference= {angle_min}, Current Orbital Velocity Difference = {orbital_v_diff_min}, Current Orbit Difference = {orbit_difference_min/1e3}, Propellant left = ={propellant_left}")   
                    
        return start_time,end_time,beta,thrust_cutoff

    def simulate_trajectory(self,n=3,n_epochs=8,start_time=25,end_time=200,beta=0.1,thrust_cutoff=5*60,delta_start_time=20,delta_end_time=150,delta_beta=0.09,delta_thrust_cutoff=3*60,dynamic_n=False):
        print("Finding optimal parameters for the rocket")
        time.sleep(0.5)
        #mins=self.epoch(n=n,n_epochs=n_epochs,start_time=dynamic_n=dynamic_n)
        mins=self.epoch(start_time=start_time,end_time=end_time,beta=beta,thrust_cutoff=thrust_cutoff,n=n,delta_start_time=delta_start_time,delta_end_time=delta_end_time,delta_beta=delta_beta,delta_thrust_cutoff=delta_thrust_cutoff,n_epochs=n_epochs,dynamic_n=dynamic_n)
        self.mins=mins
        print("Optimal Parameters = ",mins,"\n")
        #self.model(vector_start_time=mins[0],vector_end_time=mins[1],beta_max=mins[2],thrust_cutoff_time=mins[3],post_burnout=True,display_breakdown=True)

obj= step_one_boostback_path_planning(
    data={
        "rocket name": "Test Rocket",
        "initial position": np.array([6039478.38247853,  580788.51691171, 2138691.4636333]),  
        "initial velocity": np.array([356.8940986,5314.96084655,126.38282874]),  
        "initial acceleration": np.array([5.04961651, 35.44339106,  1.78816299]),  
        "isp": 350,  
        "single engine thrust": 2400e3,  
        "step masses": 128130+70286,  
        "propellant mass": 190000,  
        "q": np.array([0.67876754, -0.45406834,  0.47970547,  0.32090378]),  
    })  
"""
obj.model(vector_start_time=0,vector_end_time=60,beta_max=0.01,beta_max_y=0,thrust_cutoff_time=0,dt=0.1)
#obj.plotter()
obj.plot_magnitudes(altitude=True,aoa=True)
print(obj.absolute_position_difference,obj.total_rotation)
#print(obj.total_delta_v)
"""

#START TIME HAS BEEN FIXED REMEMBR TO CHANGE IT

"""
o=path_planning(thrust_by_weight=1.1,targt_orbit=500e3,isp=350)
o.path_planner(n=60,n_epochs=1,dynamic_n=True,check_history=False)
o=path_planning(thrust_by_weight=1.2,targt_orbit=500e3,isp=350)
o.path_planner(n=60,n_epochs=1,dynamic_n=True,check_history=False)
o=path_planning(thrust_by_weight=1.3,targt_orbit=500e3,isp=350)
o.path_planner(n=60,n_epochs=1,dynamic_n=True,check_history=False)
o=path_planning(thrust_by_weight=1.4,targt_orbit=500e3,isp=350)
o.path_planner(n=60,n_epochs=1,dynamic_n=True,check_history=False)
o=path_planning(thrust_by_weight=1.5,targt_orbit=500e3,isp=350)
o.path_planner(n=60,n_epochs=1,dynamic_n=True,check_history=False)
o=path_planning(thrust_by_weight=1.6,targt_orbit=500e3,isp=350)
o.path_planner(n=60,n_epochs=1,dynamic_n=True,check_history=False)
o=path_planning(thrust_by_weight=1.7,targt_orbit=500e3,isp=350)
o.path_planner(n=60,n_epochs=1,dynamic_n=True,check_history=False)
o=path_planning(thrust_by_weight=1.8,targt_orbit=500e3,isp=350)
o.path_planner(n=60,n_epochs=1,dynamic_n=True,check_history=False)
"""
#o.plot_variations(np.linspace(2,20,1000), [60.02816091954023], [0.10903009482758622], [240.53743103448278])
#o.plot_variations([2], np.linspace(10,250,1000), [0.10903009482758622], [240.53743103448278])
#o.plot_variations([2], [60.02816091954023], np.linspace(0.01,0.2,1000), [240.53743103448278])
#o.plot_variations([2], [60.02816091954023], [0.10903009482758622], np.linspace(180,480,1000))