# class trajectory is being used for the simulation

from __init__ import *

class traj:

    def __init__(self,data):
        self.position=[np.array([R_earth*np.cos(latitude),0,R_earth*np.sin(latitude)])]
        self.omega_earth=np.array([0,0,w_earth])
        self.velocity=[np.cross(np.array([0,0,w_earth]),self.position[-1])]
        #self.velocity=[np.array([0.001,0,0])]
        self.acceleration=[np.array([0,0,0])]
        self.aoa=[0]
        self.gravity=-g*np.array([np.cos(latitude),0,np.sin(latitude)])
        self.thrust_position=np.array([0,0,-30])
        self.time=[0.0]
        self.data=data
        self.data_copy=copy.deepcopy(data)
        self.betas=[0.0]
        self.stage_masses=data['stage masses']
        self.stage_propellant_masses=data['propellant masses']
        self.inertia_matrix=[self.rocket_inertia_tensor(self.stage_masses[i],100/(i+1)**2,8) for i in range(len(self.stage_masses))]
        self.thrust_by_weights=data['thrust by weight']
        self.target_orbit=500e3
        self.number_of_stages=len(self.stage_masses)
        self.burnouts=[]
        self.burnout_times=[]
        self.q=np.array([np.cos((np.pi/2-latitude)/2),0,np.sin((np.pi/2-latitude)/2),0])
        self.orientation=[self.quaternion_to_dcm(self.q)@[0,0,1]]
        
        self.gravity_delta_v=0
        self.thrust_delta_v=0
        self.steering_loss=0
        self.drag_delta_v=0
        self.total_delta_v=0

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

    def model(self,post_burnout=False,dt=dt,vector_start_time=10,vector_end_time=120,beta_max=0.2,thrust_cutoff_time=300,total_time=total_time,display_breakdown=False):
        n=1
        out_str=f"\nRocket liftoff from {np.round(self.position[0]/1000,2)} km. Earth rotates along the Z-axis\n"
        current_stage=0
        self.thrust_by_weight=self.data['thrust by weight'][current_stage]
        self.thrust_magnitude=self.thrust_by_weight*self.stage_masses[current_stage]*g

        while (np.linalg.norm(self.position[-1])-R_earth-self.target_orbit)<0 and np.linalg.norm(self.position[-1])-R_earth>=0:
            
            if self.thrust_magnitude/self.stage_masses[current_stage]/g>G_force_limit:
                self.thrust_by_weight=G_force_limit
                self.thrust_magnitude=self.thrust_by_weight*self.stage_masses[current_stage]*g
            
            if self.time[-1]>thrust_cutoff_time or self.stage_propellant_masses[current_stage]<0:
                #if self.stage_propellant_masses[current_stage]<0:
                    #print("Thrust cutoff due to lack of fuel")
                self.thrust_magnitude=0
                #print(f"Thrust magnitude set to 0 at {self.time[-1]}s and {(np.linalg.norm(self.position[-1])-R_earth)}")

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
            omega_body_frame=np.linalg.solve(self.inertia_matrix[current_stage],torque_body_frame)
            
            R_b_to_i = self.quaternion_to_dcm(self.q)
            thrust_eci_frame = R_b_to_i @ thrust_body_frame

            gravity_eci_frame=-g*(R_earth**2/np.linalg.norm(self.position[-1])**2)*self.normalize_quaternion(self.position[-1])

            self.gravity_delta_v+=self.get_component(gravity_eci_frame,self.velocity[-1])*dt
            self.thrust_delta_v+=self.get_component(thrust_eci_frame/self.stage_masses[current_stage],self.velocity[-1])*dt
            self.steering_loss+=(np.linalg.norm(thrust_eci_frame/self.stage_masses[current_stage])**2-self.get_component(thrust_eci_frame/self.stage_masses[current_stage],self.velocity[-1])**2)**0.5*dt
            self.drag_delta_v=50
            self.total_delta_v+=self.thrust_magnitude/self.stage_masses[current_stage]*dt

            self.acceleration+=[thrust_eci_frame/self.stage_masses[current_stage]+gravity_eci_frame]
            self.velocity+=[self.velocity[-1]+self.acceleration[-1]*dt]
            self.position+=[self.position[-1]+self.velocity[-1]*dt+0.5*self.acceleration[-1]*dt*dt]
            self.orientation+=[self.quaternion_to_dcm(self.q)@[0,0,1]]

            dqdt = self.quaternion_derivative(self.q, omega_body_frame)
            self.q += dqdt * dt
            self.q = self.normalize_quaternion(self.q)
            self.stage_masses[current_stage]-=self.thrust_magnitude/self.data['isp'][current_stage]/g*dt
            self.stage_propellant_masses[current_stage]-=self.thrust_magnitude/self.data['isp'][current_stage]/g*dt

            self.time+=[self.time[-1]+dt]
            
            if self.stage_propellant_masses[current_stage]<0 and n<=self.number_of_stages:
                print(f"Burnout positon (in earth centered earth fixed frame) = {self.position[-1]} km")
                print(f"Orientation vector at burnout (in earth centered earth fixed frame) = {np.round(self.orientation[-1],9)}")
                print(f"Downrange position from launch site (in earth centered earth fixed frame) = {np.round((self.position[-1]-self.position[0])/1e3,9)} km")
                print(f"Downrange distance from launch site (in earth centered earth fixed frame) = {np.linalg.norm(self.position[-1]-self.position[0])/1e3} km")
                print(f"Velocity vector at burnout (in earth centered earth fixed frame) = {np.round(self.velocity[-1],9)} m/s")
                print(f"quaternion at burnout = {np.round(self.q,9)}")
                print(f"Acceleration vector at burnout (in earth centered earth fixed frame) = {np.round(self.acceleration[-1],9)} m/s^2")
                print(f"Velocity magnitude at burnout (in earth centered earth fixed frame) = {np.linalg.norm(self.velocity[-1])} m/s")
                print(f"Altitude at burnout = {(np.linalg.norm(self.position[-1])-R_earth)/1000} km")
                #print(f"Mass of stage {current_stage+1} at burnout = {self.stage_masses[current_stage]} kg")
                out_str+=f"Stage {current_stage+1} burnout:\nat {round((np.linalg.norm(self.position[-1])-R_earth)/1e3,2)} km at {round(self.time[-1]/60,2)} minutes\n"
                self.burnouts+=[self.position[-1]]
                self.burnout_times+=[self.time[-1]]
                
                out_str+="\n"
                #current_stage+=1
                n+=1
                if current_stage+1<self.number_of_stages:
                    current_stage+=1    
                    self.thrust_by_weight=self.data['thrust by weight'][current_stage]
                    self.thrust_magnitude=self.thrust_by_weight*self.stage_masses[current_stage]*g
        
        #print(f"Out of loop, current altitude {(np.linalg.norm(self.position[-1])-R_earth)}")

        #self.burnouts=np.array(self.burnouts)
        #self.burnouts+=[self.position[-1]]
        #self.burnout_times+=[self.time[-1]]
        #out_str+="Note that the burnout altitudes are not very accurate as they depend on the actual ascent trajectory which is not yet fixed\n"
        
        
        self.angle=90-np.arccos(np.dot(self.position[-1],self.velocity[-1])/np.linalg.norm(self.position[-1])/np.linalg.norm(self.velocity[-1]))/np.pi*180
        orbital_velocity_difference=((np.linalg.norm(self.velocity[-1])-(g*(R_earth**2/np.linalg.norm(self.position[-1])**2)*R_earth**2/np.linalg.norm(self.position[-1]))**0.5))
        self.orbital_velocity_difference=orbital_velocity_difference
        self.orbit_difference=(np.linalg.norm(self.position[-1])-R_earth-self.target_orbit)
        self.propellant_left=self.stage_propellant_masses[-1]

        self.total_delta_v+=np.abs(orbital_velocity_difference)

        #print(F"Required delta v {orbital_velocity_difference}")
        available_delta_v=self.data["isp"][current_stage]*g*np.log(self.stage_masses[current_stage]/(self.stage_masses[current_stage]-self.stage_propellant_masses[current_stage]))
        
        if available_delta_v>np.abs(orbital_velocity_difference):
            self.velocity[-1]=((g*R_earth**2/np.linalg.norm(self.position[-1]))**0.5)*(np.cross(np.cross(self.position[-1],self.velocity[-1]),self.position[-1])/np.linalg.norm(self.position[-1])**2/np.linalg.norm(self.velocity[-1]))
            self.velocity[-1]=np.cross(np.cross(self.position[-1],self.velocity[-1]),self.position[-1])/np.linalg.norm(self.position[-1])**2/np.linalg.norm(self.velocity[-1])
            self.velocity[-1]*=((g*R_earth**2/np.linalg.norm(self.position[-1]))**0.5)
            self.burnouts+=[self.position[-1]]
            self.burnout_times+=[self.time[-1]]
            propellant_used=self.stage_masses[current_stage]*(1-np.exp(orbital_velocity_difference/self.data["isp"][current_stage]/g))
            if display_breakdown:
                print(f"Final orbit {(np.linalg.norm(self.position[-1])-R_earth)/1000} km")
                print(f"Final velocity magnitude {np.linalg.norm(self.velocity[-1])} m/s")
                print(f"Required velocity magnitude {(g*R_earth**2/np.linalg.norm(self.position[-1]))**0.5} m/s")
                print(f"Angle difference {90-np.arccos(np.dot(self.position[-1],self.velocity[-1])/np.linalg.norm(self.position[-1])/np.linalg.norm(self.velocity[-1]))/np.pi*180} deg")
                print(f"Ascent propellant left {self.stage_propellant_masses[current_stage]-propellant_used} kg")
                print(f"Final insertion engine fire duration {np.abs(orbital_velocity_difference)/G_force_limit/g} s")
                print(f"Total delta v {self.total_delta_v} m/s")

                out_str+="********************************************************************************************************************************"
                out_str+=f"Final orbit {(np.linalg.norm(self.position[-1])-R_earth)/1000} km\n"
                out_str+=f"Final velocity magnitude {np.linalg.norm(self.velocity[-1])} m/s\n"
                out_str+=f"Required velocity magnitude {(g*R_earth**2/np.linalg.norm(self.position[-1]))**0.5} m/s\n"
                #out_str+=f"Angle difference {90-np.arccos(np.dot(self.position[-1],self.velocity[-1])/np.linalg.norm(self.position[-1])/np.linalg.norm(self.velocity[-1]))/np.pi*180} deg\n"
                out_str+=f"Ascent propellant left {self.stage_propellant_masses[current_stage]-propellant_used} kg\n"


        self.burnouts=np.array(self.burnouts)
        if post_burnout:
            dt=dt/100
            while self.time[-1]<total_time:

                gravity_eci_frame=-g*(R_earth**2/np.linalg.norm(self.position[-1])**2)*self.normalize_quaternion(self.position[-1])
                
                self.acceleration+=[gravity_eci_frame]
                self.velocity+=[self.velocity[-1]+self.acceleration[-1]*dt]
                self.position+=[self.position[-1]+self.velocity[-1]*dt+0.5*self.acceleration[-1]*dt*dt]

                self.time+=[self.time[-1]+dt]
        
        #print(f"Final altitutde reached={np.linalg.norm(self.position[-1])-R_earth}")

        if display_breakdown:
            print(f"GRAVITY DELTA V : {self.gravity_delta_v}, STEERING LOSS : {self.steering_loss}, DRAG DELTA V : {self.drag_delta_v}, ORBITAL VELOCITY DIFFERENCE : {orbital_velocity_difference}, INSERTION ANGLE : {self.angle}")

        return out_str

    def model_2(self,post_burnout=True,dt=dt,vector_start_time=10,vector_end_time=120,beta_max=0.2,thrust_cutoff_time=300,total_time=total_time):
        
        out_str=f"\nRocket liftoff from {np.round(self.position[0]/1000,2)} km. Earth rotates along the Z-axis\n"
        for current_stage in range(0,self.number_of_stages):
            delta_v_achieved=0
            self.thrust_by_weight=self.data['thrust by weight'][current_stage]
            self.thrust_magnitude=self.thrust_by_weight*self.stage_masses[current_stage]*g

            while self.stage_propellant_masses[current_stage]>0 and self.time[-1]<thrust_cutoff_time:
                self.betas+=[0]
                if self.thrust_magnitude/self.stage_masses[current_stage]/g>4:
                    self.thrust_by_weight=4
                    self.thrust_magnitude=self.thrust_by_weight*self.stage_masses[current_stage]*g
                
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
                omega_body_frame=np.linalg.solve(self.inertia_matrix[current_stage],torque_body_frame)
                
                R_b_to_i = self.quaternion_to_dcm(self.q)
                thrust_eci_frame = R_b_to_i @ thrust_body_frame

                gravity_eci_frame=-g*self.normalize_quaternion(self.position[-1])
                #gravity_body_frame=np.linalg.solve(self.quaternion_to_dcm(self.q),gravity_eci_frame)
                #velocity_body_frame=np.linalg.solve(self.quaternion_to_dcm(self.q),self.velocity[-1])

                #self.gravity_delta_v+=self.get_component(gravity_body_frame,velocity_body_frame)*dt
                #self.thrust_delta_v+=self.get_component(thrust_body_frame/self.stage_masses[current_stage],velocity_body_frame)*dt

                self.gravity_delta_v+=self.get_component(gravity_eci_frame,self.velocity[-1])*dt
                self.thrust_delta_v+=self.get_component(thrust_eci_frame/self.stage_masses[current_stage],self.velocity[-1])*dt
                self.steering_loss+=(np.linalg.norm(thrust_eci_frame/self.stage_masses[current_stage])**2-self.get_component(thrust_eci_frame/self.stage_masses[current_stage],self.velocity[-1])**2)**0.5*dt
                
                self.acceleration+=[thrust_eci_frame/self.stage_masses[current_stage]+gravity_eci_frame]
                self.velocity+=[self.velocity[-1]+self.acceleration[-1]*dt]
                self.position+=[self.position[-1]+self.velocity[-1]*dt+0.5*self.acceleration[-1]*dt*dt]
                self.orientation+=[self.quaternion_to_dcm(self.q)@[0,0,1]]

                dqdt = self.quaternion_derivative(self.q, omega_body_frame)
                self.q += dqdt * dt
                self.q = self.normalize_quaternion(self.q)
                self.stage_masses[current_stage]-=self.thrust_magnitude/self.data['isp'][current_stage]/g*dt
                self.stage_propellant_masses[current_stage]-=self.thrust_magnitude/self.data['isp'][current_stage]/g*dt

                self.time+=[self.time[-1]+dt]
                
                
            out_str+=f"Stage {current_stage+1} burnout:\nat {round((np.linalg.norm(self.position[-1])-R_earth)/1e3,2)} km at {round(self.time[-1]/60,2)} minutes\n"
            self.burnouts+=[self.position[-1]]
            self.burnout_times+=[self.time[-1]]
            
            theta_=w_earth*self.time[-1]
            Rz = np.array([
                                [ np.cos(theta_), -np.sin(theta_), 0],
                                [np.sin(theta_), np.cos(theta_), 0],
                                [0,             0,                1]
                            ])
            
            '''out_str+=f"downrange position from launch site (in earth centered earth fixed frame) = {np.round((self.position[-1]-self.position[0])/1e3,2)} km\n"
            out_str+=f"velocity vector at burnout (in earth centered earth fixed frame) = {np.round(self.velocity[-1],2)} m/s\n"
            out_str+=f"acceleration vector at burnout (in earth centered earth fixed frame) = {np.round(self.acceleration[-1],2)} m/s^2\n"
            out_str+=f"direction vector at burnout (in earth centered earth fixed frame) = {np.round(self.orientation[-1],4)}\n"
            out_str+=f"delta-v achieved = {delta_v_achieved} m/s\n"
            '''
            out_str+="\n"
            
            
        self.burnouts=np.array(self.burnouts)
        out_str+="Note that the burnout altitudes are not very accurate as they depend on the actual ascent trajectory which is not yet fixed\n"
        print(f"GRAVITY DELTA V : {self.gravity_delta_v}, THRUST DELTA V : {self.thrust_delta_v}, STEERING LOSS : {self.steering_loss}")

        print(f"Thrust cutoff. Altitutde reached={np.linalg.norm(self.position[-1])-R_earth}.\nFuel remaining={self.stage_propellant_masses[-1]}")

        if post_burnout:
            while self.time[-1]<total_time:

                gravity_eci_frame=-g*self.normalize_quaternion(self.position[-1])
                
                self.acceleration+=[gravity_eci_frame]
                self.velocity+=[self.velocity[-1]+self.acceleration[-1]*dt]
                self.position+=[self.position[-1]+self.velocity[-1]*dt+0.5*self.acceleration[-1]*dt*dt]

                self.time+=[self.time[-1]+dt]
        
        print(f"Final altitutde reached={np.linalg.norm(self.position[-1])-R_earth}")

        return out_str

    def model_obsolete(self,post_burnout=True,dt=dt,vector_start_time=10,vector_end_time=120,beta_max=0.2,thrust_cutoff_time=300,total_time=total_time):
        
        out_str=f"\nRocket liftoff from {np.round(self.position[0]/1000,2)} km. Earth rotates along the Z-axis\n"
        for current_stage in range(0,self.number_of_stages):
            delta_v_achieved=0
            self.thrust_by_weight=self.data['thrust by weight'][current_stage]
            self.thrust_magnitude=self.thrust_by_weight*self.stage_masses[current_stage]*g

            while self.stage_propellant_masses[current_stage]>0 and self.time[-1]<thrust_cutoff_time:
                self.betas+=[0]
                if self.thrust_magnitude/self.stage_masses[current_stage]/g>4:
                    self.thrust_by_weight=4
                    self.thrust_magnitude=self.thrust_by_weight*self.stage_masses[current_stage]*g
                
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
                omega_body_frame=np.linalg.solve(self.inertia_matrix[current_stage],torque_body_frame)
                
                R_b_to_i = self.quaternion_to_dcm(self.q)
                thrust_eci_frame = R_b_to_i @ thrust_body_frame

                gravity_eci_frame=-g*self.normalize_quaternion(self.position[-1])
                #gravity_body_frame=np.linalg.solve(self.quaternion_to_dcm(self.q),gravity_eci_frame)
                #velocity_body_frame=np.linalg.solve(self.quaternion_to_dcm(self.q),self.velocity[-1])

                #self.gravity_delta_v+=self.get_component(gravity_body_frame,velocity_body_frame)*dt
                #self.thrust_delta_v+=self.get_component(thrust_body_frame/self.stage_masses[current_stage],velocity_body_frame)*dt

                self.gravity_delta_v+=self.get_component(gravity_eci_frame,self.velocity[-1])*dt
                self.thrust_delta_v+=self.get_component(thrust_eci_frame/self.stage_masses[current_stage],self.velocity[-1])*dt
                self.steering_loss+=(np.linalg.norm(thrust_eci_frame/self.stage_masses[current_stage])**2-self.get_component(thrust_eci_frame/self.stage_masses[current_stage],self.velocity[-1])**2)**0.5*dt
                
                self.acceleration+=[thrust_eci_frame/self.stage_masses[current_stage]+gravity_eci_frame]
                self.velocity+=[self.velocity[-1]+self.acceleration[-1]*dt]
                self.position+=[self.position[-1]+self.velocity[-1]*dt+0.5*self.acceleration[-1]*dt*dt]
                self.orientation+=[self.quaternion_to_dcm(self.q)@[0,0,1]]

                dqdt = self.quaternion_derivative(self.q, omega_body_frame)
                self.q += dqdt * dt
                self.q = self.normalize_quaternion(self.q)
                self.stage_masses[current_stage]-=self.thrust_magnitude/self.data['isp'][current_stage]/g*dt
                self.stage_propellant_masses[current_stage]-=self.thrust_magnitude/self.data['isp'][current_stage]/g*dt

                self.time+=[self.time[-1]+dt]
                
                #print(self.stage_propellant_masses[current_stage])
            
            print(f"Thrust cutoff. Altitutde reached={np.linalg.norm(self.position[-1])-R_earth}")
                
                
            out_str+=f"Stage {current_stage+1} burnout:\nat {round((np.linalg.norm(self.position[-1])-R_earth)/1e3,2)} km at {round(self.time[-1]/60,2)} minutes\n"
            self.burnouts+=[self.position[-1]]
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
            #out_str+=f"downrange position from launch site (in earth centered inertial frame) = {np.round((self.position_inertial_frame[-1]-Rz@self.position_inertial_frame[0])/1e3,2)} km\n"
            #out_str+=f"position of launch site at burnout (in earth centered inertial frame) = {np.round((Rz@self.position_inertial_frame[0])/1e3,2)} km\n"
            #out_str+=f"velocity vector at burnout (in earth centered inertial frame) = {np.round(self.velocity_inertial_frame[-1],2)} m/s\n"
            #out_str+=f"acceleration vector at burnout (in earth centered inertial frame) = {np.round(self.acceleration_inertial_frame[-1],2)} m/s\n"
            #out_str+=f"direction vector at burnout (in earth centered inertial frame) = {np.round(Rz@self.orientation[-1],4)}\n"
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
        #print(self.acceleration[1],self.velocity[1])
        out_str+="Note that the burnout altitudes are not very accurate as they depend on the actual ascent trajectory which is not yet fixed\n"
        #print("Note that the burnout altitudes are not very accurate as they depend on the actual ascent trajectory which is not yet fixed")

        print(f"GRAVITY DELTA V : {self.gravity_delta_v}, THRUST DELTA V : {self.thrust_delta_v}, STEERING LOSS : {self.steering_loss}")

        if post_burnout:
            while self.time[-1]<total_time:

                gravity_eci_frame=-g*self.normalize_quaternion(self.position[-1])
                
                self.acceleration+=[gravity_eci_frame]
                self.velocity+=[self.velocity[-1]+self.acceleration[-1]*dt]
                self.position+=[self.position[-1]+self.velocity[-1]*dt+0.5*self.acceleration[-1]*dt*dt]

                self.time+=[self.time[-1]+dt]

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
            self.fig.savefig(output_file)
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

    def plot_altitudes(self,show_velocity_plots=False):
        plt.figure()
        self.altitudes=[(np.linalg.norm(self.position[i])-R_earth)/1e3 for i in range(len(self.position))]
        plt.plot(self.time,self.altitudes)
        self.altitudes_burnout=[(np.linalg.norm(self.burnouts[i])-R_earth)/1e3 for i in range(len(self.burnouts))]
        plt.plot(self.time,self.altitudes)
        plt.scatter(self.burnout_times,self.altitudes_burnout)
        plt.title("Altitude vs Time")
        plt.show()
        plt.figure()
        self.velocities=[(np.linalg.norm(self.velocity[i])) for i in range(len(self.velocity))]
        plt.plot(self.time,self.velocities)
        plt.title("Velocity vs Time")
        if show_velocity_plots:
            plt.show()

    def optimize_controls(self,x,print_stats=False):
        vector_start_time=x[0]
        vector_end_time=x[1]
        beta_max=x[2]
        for current_stage in range(0,self.number_of_stages):
            delta_v_achieved=0
            self.thrust_by_weight=self.data['thrust by weight'][current_stage]
            self.thrust_magnitude=self.thrust_by_weight*self.stage_masses[current_stage]*g

            while self.stage_propellant_masses[current_stage]>0:
                self.betas+=[0]
                if self.thrust_magnitude/self.stage_masses[current_stage]/g>4:
                    self.thrust_by_weight=4
                    self.thrust_magnitude=self.thrust_by_weight*self.stage_masses[current_stage]*g
                
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
                omega_body_frame=np.linalg.solve(self.inertia_matrix[current_stage],torque_body_frame)
                
                R_b_to_i = self.quaternion_to_dcm(self.q)
                thrust_eci_frame = R_b_to_i @ thrust_body_frame

                gravity_eci_frame=-g*self.normalize_quaternion(self.position[-1])
                
                self.acceleration+=[thrust_eci_frame/self.stage_masses[current_stage]+gravity_eci_frame]
                self.velocity+=[self.velocity[-1]+self.acceleration[-1]*dt]
                self.position+=[self.position[-1]+self.velocity[-1]*dt+0.5*self.acceleration[-1]*dt*dt]

                dqdt = self.quaternion_derivative(self.q, omega_body_frame)
                self.q += dqdt * dt
                self.q = self.normalize_quaternion(self.q)
                self.stage_masses[current_stage]-=self.thrust_magnitude/self.data['isp'][current_stage]/g*dt
                self.stage_propellant_masses[current_stage]-=self.thrust_magnitude/self.data['isp'][current_stage]/g*dt
                self.time+=[self.time[-1]+dt]
        
        altitudes=[np.linalg.norm(self.position[i])-R_earth for i in range(0,len(self.position))]
        average_orbit=np.average(altitudes[-30:-1])
        final_orbit=altitudes[-1]

        if print_stats:
            print(f"\nFinal Orbit : {final_orbit}\t\tTarget Orbit : {self.data["target orbit"]*1e3}\nFinal Velocity : {np.linalg.norm(self.velocity[-1])}\tRequired Velocity : {(g*R_earth**2/np.linalg.norm(self.position[-1]))**0.5}\nOffset Angle : {np.arccos(np.dot(self.position[-1],self.velocity[-1])/np.linalg.norm(self.velocity[-1])/np.linalg.norm(self.position[-1]))/np.pi*180}")

        self.__init__(self.data_copy)
        if beta_max<0 or vector_end_time>250 or vector_start_time<15:
            return 1e10,1e10,1e10

        #return (np.linalg.norm(self.velocity[-1])-(g*R_earth**2/np.linalg.norm(self.position[-1]))**0.5)*(np.arccos(np.dot(self.position[-1],self.velocity[-1])/np.linalg.norm(self.velocity[-1])/np.linalg.norm(self.position[-1]))/np.pi*180),final_orbit-self.data["target orbit"]*1e3#np.dot(self.position[-1],self.velocity[-1])
        return (np.linalg.norm(self.velocity[-1])-(g*R_earth**2/np.linalg.norm(self.position[-1]))**0.5),np.dot(self.position[-1],self.velocity[-1])*1000,final_orbit-self.data["target orbit"]*1e3

    def sweep_controls(self,start_time_array,end_time_array,beta_array):
        print(f"\nSweeping through control values")
        a=0
        plt.figure()
        for i in start_time_array:
            for j in end_time_array:
                for k in beta_array:
                    self.model(vector_start_time=i,vector_end_time=j,beta_max=k)
                    self.altitudes=[(np.linalg.norm(self.position[i])-R_earth)/1e3 for i in range(len(self.position))]
                    plt.plot(self.time,self.altitudes)
                    a+=1
                    print(f"{int(a/len(start_time_array)/len(end_time_array)/len(beta_array)*100)}% completed")
                    self.__init__(self.data_copy)
        plt.show()

    def sweep(self,start_time_array,end_time_array,beta_array,thrust_cutoff_array):

        angle_min=1e10
        orbital_velocity_difference_min=1e10
        total_delta_v=1e10
        orbit_difference_min=0
        propellant_left=0
        start_time_min=0
        end_time_min=0
        beta_min=0
        thrust_cutoff_min=0

        with tqdm.tqdm(total=len(start_time_array)*len(end_time_array)*len(beta_array)*len(thrust_cutoff_array)) as pbar:
            for i in start_time_array:
                for j in end_time_array:
                    for k in beta_array:
                        for l in thrust_cutoff_array:
                            self.model(vector_start_time=i,vector_end_time=j,beta_max=k,thrust_cutoff_time=l)
                            x=np.abs(self.angle)
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

class step_one_boostback_traj:

    def __init__(self,data):
        self.position=[data['initial position']]
        self.omega_earth=[np.array([0,0,w_earth])]
        self.velocity=[data['initial velocity']]
        #self.velocity=[np.array([0.001,0,0])]
        self.acceleration=[data['initial acceleration']]
        self.aoa=[0]
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

    def normalize(self,v):
        return v / np.linalg.norm(v)
    
    def get_component(self,a,b):
        return np.dot(a,b)/np.linalg.norm(b)

    def model(self,post_burnout=False,dt=dt,vector_start_time=10,vector_end_time=120,beta_max=0.2,beta_max_y=0.01,thrust_cutoff_time=300,total_time=total_time,display_breakdown=False):
        n=1
        out_str=f"\nRocket liftoff from {np.round(self.position[0]/1000,2)} km. Earth rotates along the Z-axis\n"
        #current_stage=0
        #self.thrust_by_weight=self.data['thrust by weight']
        #self.thrust_magnitude=self.thrust_by_weight*self.step_mass*g
        print((np.linalg.norm(self.position[-1])-R_earth)>0)
        while (np.linalg.norm(self.position[-1])-R_earth)>0:# and np.linalg.norm(self.position[-1])-R_earth>=0:
            if self.time[-1]>=vector_start_time:
                #if self.stage_propellant_masses[current_stage]<0:
                    #print("Thrust cutoff due to lack of fuel")
                print("flip maneuver startd at ",self.time[-1])
                self.thrust_magnitude= self.step_thrust_magnitude
                #print(f"Thrust magnitude set to 0 at {self.time[-1]}s and {(np.linalg.norm(self.position[-1])-R_earth)}")
            
            if self.time[-1]>vector_end_time and np.linalg.norm(self.position[-1])-R_earth>60e3:
                self.thrust_magnitude=0
                print("flip maneuver ended at ",self.time[-1])
                #print(f"Thrust magnitude set to 0 at {self.time[-1]}s and {(np.linalg.norm(self.position[-1])-R_earth)}")
            
            if np.linalg.norm(self.position[-1])-R_earth<=60e3:
                if self.orientation_flipped==False:
                    print(f"Flipping orientation at {self.time[-1]}s")
                    print(f"Original orientation: {self.orientation[-1]}")
                    print(f"Original quaternion: {self.q}")
                    Rxflip = np.array([ [1      , 0            , 0],
                            [0      , np.cos(pi), -np.sin(pi)],
                            [0      , np.sin(pi),  np.cos(pi)]
                        ])
                    self.orientation[-1]=Rxflip@self.orientation[-1]
                    self.q =self.quat_mul(np.array([0,1,0,0]),self.q) 
                    self.q = self.normalize_quaternion(self.q)
                    self.orientation_flipped=True
                    print(f"New orientation: {self.orientation[-1]}")
                    print(f"New quaternion: {self.q}")
                #print(f"Thrust cutoff due to altitude {np.linalg.norm(self.position[-1])-R_earth}")
                if self.velocity_reduction_boostback_duration<thrust_cutoff_time:
                    self.thrust_magnitude=self.single_engine_thrust
                    self.velocity_reduction_boostback_duration+=dt
                    print("velocity reduction started at",self.time[-1])
                else:
                    self.thrust_magnitude=0
                    print("velocity reduction ended at",self.time[-1])
                    #print(f"Thrust magnitude set to 0 at {self.time[-1]}s and {(np.linalg.norm(self.position[-1])-R_earth)}")  
                #print(f"Thrust magnitude set to 0 at {self.time[-1]}s and {(np.linalg.norm(self.position[-1])-R_earth)}")
            
            if np.linalg.norm(self.position[-1])-R_earth<2e3:
                self.thrust_magnitude=self.single_engine_thrust

            if self.step_propellant_mass<=0 and self.thrust_magnitude>0:
                print("Fuel depleted, thrust magnitude set to 0")
                self.thrust_magnitude=0

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
            torque_body_frame=np.cross(self.thrust_position,thrust_body_frame)
            omega_body_frame=np.linalg.solve(self.inertia_matrix,torque_body_frame)
            
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

            dqdt = self.quaternion_derivative(self.q, omega_body_frame)
            self.q += dqdt * dt
            self.q = self.normalize_quaternion(self.q)
            self.step_mass-=self.thrust_magnitude/self.data['isp']/g*dt
            self.step_propellant_mass-=self.thrust_magnitude/self.data['isp']/g*dt
            self.inertia_matrix=self.rocket_inertia_tensor(self.step_mass,70,8)

            self.time+=[self.time[-1]+dt]

        if display_breakdown:
            print(f"GRAVITY DELTA V : {self.gravity_delta_v}, STEERING LOSS : {self.steering_loss}, DRAG DELTA V : {self.drag_delta_v}, ORBITAL VELOCITY DIFFERENCE : {orbital_velocity_difference}, INSERTION ANGLE : {self.angle}")

        self.burnouts+=[self.position[-1]]
        self.burnout_times+=[self.time[-1]]
        self.burnouts=np.array(self.burnouts)

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

    def plot_altitudes(self,show_velocity_plots=False):
        plt.figure()
        self.altitudes=[(np.linalg.norm(self.position[i])-R_earth)/1e3 for i in range(len(self.position))]
        plt.plot(self.time,self.altitudes)
        self.altitudes_burnout=[(np.linalg.norm(self.burnouts[i])-R_earth)/1e3 for i in range(len(self.burnouts))]
        plt.plot(self.time,self.altitudes)
        plt.scatter(self.burnout_times,self.altitudes_burnout)
        plt.title("Altitude vs Time")
        plt.show()
        plt.figure()
        self.velocities=[(np.linalg.norm(self.velocity[i])) for i in range(len(self.velocity))]
        plt.plot(self.time,self.velocities)
        plt.title("Velocity vs Time")
        if show_velocity_plots:
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

obj= step_one_boostback_traj(
    data={
        "rocket name": "Test Rocket",
        "initial position": np.array([5671585.06898223,  401116.70863152, 3143810.93982123]),  # 100 km above Earth's surface
        "initial velocity": np.array([744.49624513, 4914.2327531,   412.68100745]),  # Starting at rest
        "initial acceleration": np.array([2.03135836, 36.88199574,  1.12600033]),  # No initial acceleration
        "isp": 350,  # Specific impulse in seconds
        "single engine thrust": 2400e3,  # Thrust of a single engine
        "step masses": 128130+70286,  # Mass of the step in kg
        "propellant mass": 150000,  # Mass of the propellant in kg
        "q": np.array([0.69470791, -0.5096918,   0.40921423,  0.30023142]),  # Initial quaternion (no rotation)
    })  
obj.model(vector_end_time=1000,beta_max=0.5)
obj.plotter()
obj.plot_altitudes()
