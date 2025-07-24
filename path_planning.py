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
    

        x0=np.array((20.000064, 50.0, 0.14553971200000002, 210.0))
        bounds=((5,50),(50,300),(0.05,0.2),(3*60,8*60))
        constraints=[ {'type':'ineq','fun':self.constraint2}]
        solution=minimize(self.objective_function,x0,bounds=bounds,constraints=constraints)
        print(solution.success)
        print(solution.x)
        return solution.x

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
                #print(np.linalg.norm(self.position))
                self.model(vector_start_time=mins[0],vector_end_time=mins[1],beta_max=mins[2],thrust_cutoff=mins[3],post_burnout=True,display_breakdown=True)
                #self.plot_altitudes()
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
        
#START TIME HAS BEEN FIXED REMEMBR TO CHANGE IT


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
#o.plot_variations(np.linspace(2,20,1000), [60.02816091954023], [0.10903009482758622], [240.53743103448278])
#o.plot_variations([2], np.linspace(10,250,1000), [0.10903009482758622], [240.53743103448278])
#o.plot_variations([2], [60.02816091954023], np.linspace(0.01,0.2,1000), [240.53743103448278])
#o.plot_variations([2], [60.02816091954023], [0.10903009482758622], np.linspace(180,480,1000))