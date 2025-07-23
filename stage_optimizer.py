"""
epsilon (Structural factor) is the ratio of structural mass to total mass (inert mass+fuel mass) of a step. Total mass does not include payload mass
"""

from __init__ import *

class stage_optimizer():

    def __init__(self,delta_v,delta_v_boostback,payload,number_of_stages,Isps,structural_factors,boostback=["False","False"]):
        self.boostback=boostback
        self.payload=payload
        self.number_of_stages=number_of_stages
        self.Isps=np.zeros(number_of_stages)
        self.epsilons=np.zeros(number_of_stages)
        self.epsilons_ascent=np.zeros(number_of_stages)
        self.epsilon_boostback=np.zeros(number_of_stages)
        self.epsilons_boostback=np.zeros(number_of_stages)
        self.delta_v_boostback=np.zeros(number_of_stages)
        self.stage_masses=np.zeros(number_of_stages)
        self.step_masses=np.zeros(number_of_stages)
        self.alpha_min=0
        self.alpha_max=0
        self.alpha=(self.alpha_min+self.alpha_max)/2
        self.Isps=Isps
        self.epsilons=structural_factors*1
        for i in range(0,self.number_of_stages):
            if not self.boostback[i]:
                self.delta_v=delta_v
                self.epsilons_ascent[i]=self.epsilons[i]
                self.epsilons_boostback[i]=1
            else:
                #print(self.boostback[i])
                self.delta_v=delta_v
                self.delta_v_boostback[i]=delta_v_boostback[i]
                self.epsilon_boostback[i]=np.exp(-self.delta_v_boostback[i]/self.Isps[i]/g)
                #self.epsilons_boostback=np.ones_like(self.epsilons)
                self.epsilons_boostback[i]=self.epsilon_boostback[i]
                #self.epsilons_ascent=self.epsilons*1
                self.epsilons_ascent[i]=self.epsilons[i]/self.epsilon_boostback[i]
                self.epsilons[i]=self.epsilons[i]/self.epsilon_boostback[i]
                #print(self.epsilons)
        #print(self.epsilons_boostback,self.epsilons_ascent,self.epsilons)
        self.calculate_betas(self.alpha)
        alpha=fsolve(self.solve_for_betas_2,0)
        self.delta_v_achieved=self.calculate_delta_v()
        self.epsilons=structural_factors*1
        #self.calculate_betas(alpha)
        #self.epsilons=structural_factors
        #print(self.betas)
        self.calculate_rocket_mass()
        #self.display_results()
        self.tabulate_mass_breakdown()
        self.rocket=rocket.rocket(self.payload)
        self.create_rocket_object()
    
    def calculate_betas(self,alpha):
        self.betas=self.epsilons/((1-self.epsilons)*(1+alpha*self.Isps*g))
        return

    def calculate_delta_v(self):
        return np.sum(-self.Isps*g*np.log(self.betas+(1-self.betas)*self.epsilons))

    def solve_for_betas_2(self,alpha):
        self.calculate_betas(alpha)
        return self.delta_v-self.calculate_delta_v()
            
    def calculate_rocket_mass(self):
        #print(self.betas)
        self.step_masses[-1]=self.payload/self.betas[-1]-self.payload
        self.stage_masses[-1]=self.payload/self.betas[-1]-self.payload
        #print(self.stage_masses[-1])
        for i in range(1,self.number_of_stages):
            payload=self.payload+np.sum(self.step_masses[self.number_of_stages-i:])
            self.step_masses[self.number_of_stages-i-1]=payload/self.betas[self.number_of_stages-i-1]-payload
            self.stage_masses[self.number_of_stages-i-1]=payload/self.betas[self.number_of_stages-i-1]-self.payload
            #print(self.payload,payload,self.step_masses,self.stage_masses)
            #print(self.stage_masses[self.number_of_stages-i-1])
        return

    def display_results(self):
        #print()
        print()
        print("###################################   ROCKET MASS BREAKDOWN   ##############################")
        print()
        #print()
        #print(self.epsilons_ascent,self.epsilons_boostback)
        print("====================================================================================================")
        print("Step/Stage \t\t\t\t\t",end="")
        for i in range(self.number_of_stages):
            print(i+1,"\t\t\t",end="")
        print()

        print("----------------------------------------------------------------------------------------------------")
        
        print("Step Total Mass \t\t\t\t", end="")
        for i in range(self.number_of_stages):
            print(round(self.step_masses[i],0),"\t\t",end="")

        print()

        print("Step Structural Mass \t\t\t\t",end="")
        for i in range(self.number_of_stages):
            print(round(self.step_masses[i]*(self.epsilons[i]),0),"\t\t",end="")
        print()
        
        print("Step Propellant Mass \t\t\t\t", end="")
        for i in range(self.number_of_stages):
            print(round(self.step_masses[i]*(1-self.epsilons[i]),0),"\t\t",end="")
        
        print()

        if self.boostback:

            print("Fuel Mass Used For Ascent \t\t\t", end="")
            for i in range(self.number_of_stages):
                print(round(self.step_masses[i]*(1-self.epsilons_ascent[i]),0),"\t\t",end="")

            print()

            print("Fuel Mass Used For Boostback \t\t\t", end="")
            for i in range(self.number_of_stages):
                print(round(self.step_masses[i]*self.epsilons_ascent[i]*(1-self.epsilons_boostback[i]),0),"\t\t",end="")

            print()

        if not self.boostback:
            print("Step Delta-v  \t\t\t\t\t", end="")
            for i in range(self.number_of_stages):
                print(round(-self.Isps[i]*g*np.log(self.betas[i]+(1-self.betas[i])*self.epsilons_ascent[i]),0),"\t\t\t",end="")

            print()

        if self.boostback:

            print("Step Delta-v During Ascent  \t\t\t", end="")
            for i in range(self.number_of_stages):
                print(round(-self.Isps[i]*g*np.log(self.betas[i]+(1-self.betas[i])*self.epsilons_ascent[i]),0),"\t\t\t",end="")

            print()
            print("Step Delta-v for Boostback \t\t\t", end="")
            for i in range(self.number_of_stages):
                print(round(-self.Isps[0]*g*np.log(self.epsilons_boostback[i]),2),"\t\t\t",end="")

            print()

        print("----------------------------------------------------------------------------------------------------")

        print("Stage Total Mass \t\t\t\t",end="")
        for i in range(self.number_of_stages):
            print(round(self.stage_masses[i],0),"\t\t",end="")
        print()
        print("====================================================================================================")
        if not self.boostback:
            print("Delta-v achieved =",round(self.calculate_delta_v(),2))
        else:
            print("Delta-v achieved =",round(self.delta_v_achieved-self.Isps[0]*g*np.log(self.epsilon_boostback),2))
            print("Boostback Delta-v achieved =",round(-self.Isps[0]*g*np.log(self.epsilon_boostback),2))
        print("====================================================================================================")


    def tabulate_mass_breakdown(self):
        out_str=""
        columns = [f"{i+1}" for i in range(self.number_of_stages)]
        data = {}

        # Step masses
        data["Step Total Mass (kg)"] = [int(self.step_masses[i]) for i in range(self.number_of_stages)]
        data["Step Structural Mass (kg)"] = [int(self.step_masses[i] * self.epsilons[i]) for i in range(self.number_of_stages)]
        data["Step Propellant Mass (kg)"] = [int(self.step_masses[i] * (1 - self.epsilons[i])) for i in range(self.number_of_stages)]

        if self.boostback:
            data["Fuel Mass for Ascent (kg)"] = [int(self.step_masses[i] * (1 - self.epsilons_ascent[i])) for i in range(self.number_of_stages)]
            data["Fuel Mass for Boostback (kg)"] = [int(self.step_masses[i] * self.epsilons_ascent[i] * (1 - self.epsilons_boostback[i])) for i in range(self.number_of_stages)]
        else:
            data["Delta-v (m/s)"] = [
                int(-self.Isps[i]*g*np.log(self.betas[i] + (1 - self.betas[i]) * self.epsilons_ascent[i]))
                for i in range(self.number_of_stages)
            ]

        if self.boostback:
            data["Delta-v Ascent (m/s)"] = [
                int(-self.Isps[i]*g*np.log(self.betas[i] + (1 - self.betas[i]) * self.epsilons_ascent[i]))
                for i in range(self.number_of_stages)
            ]
            data["Delta-v Boostback (m/s)"] = [
                int(-self.Isps[0]*g*np.log(self.epsilons_boostback[i])) for i in range(self.number_of_stages)
            ]

        data["Stage Total Mass (kg)"] = [round(self.stage_masses[i], 2) for i in range(self.number_of_stages)]

        df = pd.DataFrame(data, index=columns).T
        out_str+="\n################################# ROCKET MASS & Delta-v BREAKDOWN ###################################\n\n"
        #print("\n################################# ROCKET MASS & Delta-v BREAKDOWN ###################################\n")
        out_str+="=====================================================================================\n"
        #print("=====================================================================================")
        out_str+=df.to_markdown()+"\n"
        #print(df.to_markdown())

        # Print total delta-v summary
        out_str+="\n===================================================================================\n"
        #print("\n===================================================================================")
        if self.boostback[0] or self.boostback[1]:
            del_v_b=0
            out_str+=f"Delta-v Available (excluding boostback): {round(self.delta_v_achieved, 2)} m/s\n"
            for i in range(0,self.number_of_stages):
                if self.boostback[i]:
                    out_str+=f"Boostback Delta-v Achieved for stage {i+1}: {round(-self.Isps[i]*g*np.log(self.epsilon_boostback[i]), 2)} m/s\n"
                    del_v_b+=-(self.Isps[i]*g*np.log(self.epsilon_boostback[i]))
            out_str+=f"Total Delta-v Achieved: {round(self.delta_v_achieved+del_v_b, 2)} m/s\n"
        else:
            out_str+=f"Total Delta-v Achieved: {round(self.calculate_delta_v(), 2)} m/s\n"
        out_str+="===================================================================================\n"
        return out_str

    def create_rocket_object(self):
        for i in range(0,self.number_of_stages):
            self.rocket.add_stage(self.step_masses[i]*(1-self.epsilons[i]),self.stage_masses[i],self.step_masses[i],self.Isps[i],self.step_masses[i]*(1-self.epsilons_ascent[i]))


    def create_rocket_json(self, filename):
        # Create a dictionary to hold the results
        results = {
            "rocket_mass_breakdown": {
                "step_total_mass": [round(mass, 0) for mass in self.step_masses],
                "step_structural_mass": [round(self.step_masses[i] * self.epsilons[i], 0) for i in range(self.number_of_stages)],
                "step_propellant_mass": [round(self.step_masses[i] * (1 - self.epsilons[i]), 0) for i in range(self.number_of_stages)],
                "stage_total_mass": [round(mass, 0) for mass in self.stage_masses]
            }
        }
        
        # If boostback is used, add additional info to the results dictionary
        if self.boostback:
            results["boostback"] = {
                "fuel_mass_used_for_ascent": [round(self.step_masses[i] * (1 - self.epsilons_ascent[i]), 0) for i in range(self.number_of_stages)],
                "fuel_mass_used_for_boostback": [round(self.step_masses[i] * self.epsilons_ascent[i] * (1 - self.epsilons_boostback[i]), 0) for i in range(self.number_of_stages)],
            }

        # Add delta-v calculations to the results dictionary
        if not self.boostback:
            results["delta_v"] = {
                "step_delta_v": [round(-self.Isps[i] * g * np.log(self.betas[i] + (1 - self.betas[i]) * self.epsilons_ascent[i]), 0) for i in range(self.number_of_stages)],
                "delta_v_achieved": round(self.calculate_delta_v(), 2),
            }
        else:
            results["delta_v"] = {
                "step_delta_v_during_ascent": [round(-self.Isps[i] * g * np.log(self.betas[i] + (1 - self.betas[i]) * self.epsilons_ascent[i]), 0) for i in range(self.number_of_stages)],
                "boostback_delta_v": [round(-self.Isps[0] * g * np.log(self.epsilons_boostback[i]), 2) for i in range(self.number_of_stages)],
                "delta_v_achieved": round(self.delta_v_achieved - self.Isps[0] * g * np.log(self.epsilon_boostback), 2),
                "boostback_delta_v_achieved": round(-self.Isps[0] * g * np.log(self.epsilon_boostback), 2),
            }
        

        # Save the results to a JSON file
        with open(filename, 'w') as json_file:
            json.dump(results, json_file, indent=4)



class stage_optimizer_2:

    def __init__(self,delta_v,payload,number_of_stages,Isps,structural_factors):
        self.delta_v=delta_v
        self.payload=payload
        self.number_of_stages=number_of_stages
        self.Isps=np.zeros(number_of_stages)
        self.epsilons=np.zeros(number_of_stages)
        self.stage_masses=np.zeros(number_of_stages)
        self.step_masses=np.zeros(number_of_stages)
        self.alpha_min=0
        self.alpha_max=0
        self.alpha=(self.alpha_min+self.alpha_max)/2
        self.Isps=Isps
        self.epsilons=structural_factors
        '''
        for i in range(0,number_of_stages):
            self.Isps[i]=input('ENTER ISP FOR STAGE NUMBER '+str(i+1)+':')
            self.epsilons[i]=input('ENTER STRUCTURAL FACTOR FOR STAGE NUMBER '+str(i+1)+':')
        '''
        
        #self.plotter()
        #self.calculate_betas(-0.001)
        #print(self.stage_masses)
        self.calculate_betas(self.alpha)
        #self.solve_for_betas()
        alpha=fsolve(self.solve_for_betas_2,0)
        #print(alpha)
        self.calculate_betas(alpha)
        print(self.betas)
        self.calculate_rocket_mass()
        self.display_results()
        self.rocket=rocket.rocket(self.payload)
        self.create_rocket_object()
        #self.create_roccket_object()
        #print(self.stage_masses)
        #print(self.step_masses)
        #print(self.betas)
    
    def calculate_betas(self,alpha):
        self.betas=self.epsilons/((1-self.epsilons)*(1+alpha*self.Isps*g))
        return

    def calculate_delta_v(self):
        return np.sum(-self.Isps*g*np.log(self.betas+(1-self.betas)*self.epsilons))

    def solve_for_betas(self):
        i=0
        d_alpha=0.1
        self.alpha=10
        self.calculate_betas(self.alpha)
        flag=self.delta_v<self.calculate_delta_v()
        if not flag:
            #print("WARNING: GIVEN STRUCTURAL RATIOS DO NOT PRODUCE SUFFICIENT DELTA-V")
            raise Exception("WARNING: GIVEN STRUCTURAL RATIOS DO NOT PRODUCE SUFFICIENT DELTA-V. ADDITIONAL DELTA-V REQUIRED=",self.delta_v-self.calculate_delta_v())
        self.alpha=0
        self.calculate_betas(self.alpha)
        while self.delta_v-self.calculate_delta_v()>1 and flag: 
            self.alpha+=(np.abs(self.calculate_delta_v()-self.delta_v))*1e-9
            '''diff=self.calculate_delta_v()-self.delta_v
            self.alpha=self.alpha+d_alpha
            self.calculate_betas(self.alpha)
            if (self.calculate_delta_v()-self.delta_v)>diff:
                self.alpha-=d_alpha
                d_alpha/=10
            i+=1'''
            #print(self.calculate_delta_v()-self.delta_v,self.alpha)
            self.calculate_betas(self.alpha)
        print("Delta v:",self.calculate_delta_v())

    def solve_for_betas_2(self,alpha):
        self.calculate_betas(alpha)
        return self.delta_v-self.calculate_delta_v()
            
    def calculate_rocket_mass(self):
        #print(self.betas)
        self.step_masses[-1]=self.payload/self.betas[-1]-self.payload
        self.stage_masses[-1]=self.payload/self.betas[-1]-self.payload
        #print(self.stage_masses[-1])
        for i in range(1,self.number_of_stages):
            payload=self.payload+np.sum(self.step_masses[self.number_of_stages-i:])
            self.step_masses[self.number_of_stages-i-1]=payload/self.betas[self.number_of_stages-i-1]-payload
            self.stage_masses[self.number_of_stages-i-1]=payload/self.betas[self.number_of_stages-i-1]-self.payload
            #print(self.payload,payload,self.step_masses,self.stage_masses)
            #print(self.stage_masses[self.number_of_stages-i-1])
        return

    def display_results(self):

        '''
        for i in range(self.number_of_stages):
            print("====================================================================================================")
            print("STAGE/STEP ",i+1)
            print("Stage mass :",self.stage_masses[i])
            print("Step mass :",self.step_masses[i])
            print("Propellant mass :",self.step_masses[i]*(1-self.epsilons[i]))
            print("Structural mass of step ",i+1,":",self.step_masses[i]*(self.epsilons[i]))
        '''
        print("====================================================================================================")
        print("Step/Stage \t\t\t\t\t",end="")
        for i in range(self.number_of_stages):
            print(i+1,"\t\t\t",end="")
        print()
        print("----------------------------------------------------------------------------------------------------")
        
        print("Step Propellant Mass \t\t\t\t", end="")
        for i in range(self.number_of_stages):
            print(round(self.step_masses[i]*(1-self.epsilons[i]),0),"\t\t",end="")
        
        print()
        print("Step Structural Mass \t\t\t\t",end="")
        for i in range(self.number_of_stages):
            print(round(self.step_masses[i]*(self.epsilons[i]),0),"\t\t",end="")
        print()

        print("Step Total Mass \t\t\t\t", end="")
        for i in range(self.number_of_stages):
            print(round(self.step_masses[i],0),"\t\t",end="")

        print()

        print("Step Delta-v \t\t\t\t\t", end="")
        for i in range(self.number_of_stages):
            print(round(-self.Isps[i]*g*np.log(self.betas[i]+(1-self.betas[i])*self.epsilons[i]),0),"\t\t\t",end="")

        print()

        print("----------------------------------------------------------------------------------------------------")

        print("Stage Total Mass \t\t\t\t",end="")
        for i in range(self.number_of_stages):
            print(round(self.stage_masses[i],0),"\t\t",end="")
        print()
        print("====================================================================================================")
        print("Delta-v achieved =",round(self.calculate_delta_v(),2))
        print("====================================================================================================")
        
    def plotter(self):
        alphas=np.linspace(-0.1,0.1,1000)
        vel=alphas*1
        for i in range(0,1000):
            self.alpha=alphas[i]
            self.calculate_betas(alphas[i])
            vel[i]=self.calculate_delta_v()
        plt.figure()
        plt.plot(alphas,vel)
        plt.show()

    def create_rocket_object(self):
        for i in range(0,self.number_of_stages):
            self.rocket.add_stage(self.step_masses[i]*(1-self.epsilons[i]),self.stage_masses[i],self.step_masses[i],self.Isps[i])

'''        
    def solve_for_betas_gd(self):
        d_alpha=1e-2
        i=0
        while i<20:#np.abs(self.calculate_delta_v(self.alpha)-self.delta_v)>1e3:
            gradient=(self.calculate_delta_v(self.alpha)-self.calculate_delta_v(self.alpha-d_alpha))/d_alpha
            if self.calculate_delta_v(self.alpha)-self.delta_v>0:
                self.alpha=self.alpha-(self.calculate_delta_v(self.alpha)-self.delta_v)/gradient/2
            else:
                self.alpha=self.alpha+(self.calculate_delta_v(self.alpha)-self.delta_v)/gradient/2
            self.calculate_betas(self.alpha)
            print((self.calculate_delta_v(self.alpha)-self.calculate_delta_v(self.alpha-d_alpha)),gradient,np.abs(self.calculate_delta_v(self.alpha)-self.delta_v),self.alpha)
            i+=1
'''




'''
"""
epsilon (Structural factor) is the ratio of structural mass to total mass (inert mass + fuel mass) of a step.
Total mass does not include payload mass.
"""

from __init__ import *

class stage_optimizer:
    """
    Optimizes a multi-stage rocket for a given delta-v requirement and payload,
    calculating mass breakdown for each stage.
    Can optionally include a boostback maneuver.
    """
    def __init__(self, delta_v, payload, number_of_stages, Isps, structural_factors, boostback=False):
        self.boostback = boostback  # Enable/disable boostback modeling
        self.payload = payload
        self.number_of_stages = number_of_stages
        self.Isps = Isps  # List of specific impulses per stage
        self.epsilons = structural_factors * 1  # Copy of structural factors
        self.stage_masses = np.zeros(number_of_stages)
        self.step_masses = np.zeros(number_of_stages)
        self.alpha_min = 0  # Starting guess for alpha (optimization param)
        self.alpha_max = 0
        self.alpha = (self.alpha_min + self.alpha_max) / 2

        # Split delta-v into ascent and boostback if needed
        if not self.boostback:
            self.delta_v = delta_v
        else:
            self.delta_v = delta_v[0]
            self.delta_v_boostback = delta_v[1]

            # Compute epsilon for boostback based on desired delta-v
            self.epsilon_boostback = np.exp(-self.delta_v_boostback / self.Isps[0] / g)

            # Adjust structural factors to account for boostback fuel fraction
            self.epsilons_boostback = np.ones_like(self.epsilons)
            self.epsilons_boostback[0] = self.epsilon_boostback

            self.epsilons_ascent = self.epsilons * 1
            self.epsilons_ascent[0] /= self.epsilon_boostback

            self.epsilons[0] /= self.epsilon_boostback  # Adjust for total flight

        # Solve for optimal betas (mass ratios) using fsolve
        self.calculate_betas(self.alpha)
        alpha = fsolve(self.solve_for_betas_2, 0)

        # Calculate resulting delta-v and rocket configuration
        self.delta_v_achieved = self.calculate_delta_v()
        self.epsilons = structural_factors * 1  # Reset original
        self.calculate_rocket_mass()
        self.tabulate_mass_breakdown()
        self.rocket = rocket.rocket(self.payload)
        self.create_rocket_object()

    def calculate_betas(self, alpha):
        # Computes beta values from structural efficiency and alpha
        self.betas = self.epsilons / ((1 - self.epsilons) * (1 + alpha * self.Isps * g))

    def calculate_delta_v(self):
        # Calculates total achievable delta-v with current betas
        return np.sum(-self.Isps * g * np.log(self.betas + (1 - self.betas) * self.epsilons))

    def solve_for_betas_2(self, alpha):
        # Root-finding function to match target delta-v
        self.calculate_betas(alpha)
        return self.delta_v - self.calculate_delta_v()

    def calculate_rocket_mass(self):
        # Back-calculates stage masses from final payload using beta
        self.step_masses[-1] = self.payload / self.betas[-1] - self.payload
        self.stage_masses[-1] = self.step_masses[-1]
        for i in range(1, self.number_of_stages):
            payload = self.payload + np.sum(self.step_masses[self.number_of_stages - i:])
            self.step_masses[self.number_of_stages - i - 1] = payload / self.betas[self.number_of_stages - i - 1] - payload
            self.stage_masses[self.number_of_stages - i - 1] = self.step_masses[self.number_of_stages - i - 1]

    def tabulate_mass_breakdown(self):
        # Displays a detailed breakdown of each stage
        columns = [f"{i + 1}" for i in range(self.number_of_stages)]
        data = {
            "Step Total Mass (kg)": [round(m, 2) for m in self.step_masses],
            "Step Structural Mass (kg)": [round(m * e, 2) for m, e in zip(self.step_masses, self.epsilons)],
            "Step Propellant Mass (kg)": [round(m * (1 - e), 2) for m, e in zip(self.step_masses, self.epsilons)],
            "Stage Total Mass (kg)": [round(m, 2) for m in self.stage_masses],
        }

        # Add boostback data if enabled
        if self.boostback:
            data["Fuel Mass for Ascent (kg)"] = [round(m * (1 - e), 2) for m, e in zip(self.step_masses, self.epsilons_ascent)]
            data["Fuel Mass for Boostback (kg)"] = [round(m * e_a * (1 - e_b), 2) for m, e_a, e_b in zip(self.step_masses, self.epsilons_ascent, self.epsilons_boostback)]
            data["Delta-v Ascent (m/s)"] = [round(-Isp * g * np.log(beta + (1 - beta) * eps), 2) for Isp, beta, eps in zip(self.Isps, self.betas, self.epsilons_ascent)]
            data["Delta-v Boostback (m/s)"] = [round(-self.Isps[0] * g * np.log(eps), 2) for eps in self.epsilons_boostback]
        else:
            data["Delta-v (m/s)"] = [round(-Isp * g * np.log(beta + (1 - beta) * eps), 2) for Isp, beta, eps in zip(self.Isps, self.betas, self.epsilons)]

        # Create and print table
        df = pd.DataFrame(data, index=columns).T
        print("\n################################# ROCKET MASS & Delta-v BREAKDOWN ###################################\n")
        print(df.to_markdown())
        print("===================================================================================")
        if self.boostback:
            print(f"Delta-v Achieved (excluding boostback): {round(self.delta_v_achieved - self.Isps[0] * g * np.log(self.epsilon_boostback), 2)} m/s")
            print(f"Boostback Delta-v Achieved: {round(-self.Isps[0] * g * np.log(self.epsilon_boostback), 2)} m/s")
        else:
            print(f"Total Delta-v Achieved: {round(self.calculate_delta_v(), 2)} m/s")
        print("===================================================================================")
        return df

    def create_rocket_object(self):
        # Adds stages to rocket object with current stage parameters
        for i in range(self.number_of_stages):
            self.rocket.add_stage(
                self.step_masses[i] * (1 - self.epsilons[i]),  # propellant
                self.stage_masses[i],                          # gross mass
                self.step_masses[i],                           # step-only mass
                self.Isps[i]
            )


class stage_optimizer_2:
    """
    Alternate version of stage_optimizer using a more console-style output format.
    """
    def __init__(self, delta_v, payload, number_of_stages, Isps, structural_factors):
        self.delta_v = delta_v
        self.payload = payload
        self.number_of_stages = number_of_stages
        self.Isps = Isps
        self.epsilons = structural_factors
        self.stage_masses = np.zeros(number_of_stages)
        self.step_masses = np.zeros(number_of_stages)
        self.alpha_min = 0
        self.alpha_max = 0
        self.alpha = (self.alpha_min + self.alpha_max) / 2

        self.calculate_betas(self.alpha)
        alpha = fsolve(self.solve_for_betas_2, 0)
        self.calculate_betas(alpha)
        print(self.betas)
        self.calculate_rocket_mass()
        self.display_results()
        self.rocket = rocket.rocket(self.payload)
        self.create_rocket_object()

    def calculate_betas(self, alpha):
        self.betas = self.epsilons / ((1 - self.epsilons) * (1 + alpha * self.Isps * g))

    def calculate_delta_v(self):
        return np.sum(-self.Isps * g * np.log(self.betas + (1 - self.betas) * self.epsilons))

    def solve_for_betas_2(self, alpha):
        self.calculate_betas(alpha)
        return self.delta_v - self.calculate_delta_v()

    def calculate_rocket_mass(self):
        self.step_masses[-1] = self.payload / self.betas[-1] - self.payload
        self.stage_masses[-1] = self.step_masses[-1]
        for i in range(1, self.number_of_stages):
            payload = self.payload + np.sum(self.step_masses[self.number_of_stages - i:])
            self.step_masses[self.number_of_stages - i - 1] = payload / self.betas[self.number_of_stages - i - 1] - payload
            self.stage_masses[self.number_of_stages - i - 1] = self.step_masses[self.number_of_stages - i - 1]

    def display_results(self):
        # Print mass breakdown and delta-v per stage in plain text format
        print("====================================================================================================")
        print("Step/Stage \t\t\t\t\t", end="")
        for i in range(self.number_of_stages):
            print(i + 1, "\t\t\t", end="")
        print("\n----------------------------------------------------------------------------------------------------")
        print("Step Propellant Mass \t\t\t\t", end="")
        for i in range(self.number_of_stages):
            print(round(self.step_masses[i] * (1 - self.epsilons[i]), 0), "\t\t", end="")
        print("\nStep Structural Mass \t\t\t\t", end="")
        for i in range(self.number_of_stages):
            print(round(self.step_masses[i] * self.epsilons[i], 0), "\t\t", end="")
        print("\nStep Total Mass \t\t\t\t", end="")
        for i in range(self.number_of_stages):
            print(round(self.step_masses[i], 0), "\t\t", end="")
        print("\nStep Delta-v \t\t\t\t\t", end="")
        for i in range(self.number_of_stages):
            print(round(-self.Isps[i] * g * np.log(self.betas[i] + (1 - self.betas[i]) * self.epsilons[i]), 0), "\t\t\t", end="")
        print("\n----------------------------------------------------------------------------------------------------")
        print("Stage Total Mass \t\t\t\t", end="")
        for i in range(self.number_of_stages):
            print(round(self.stage_masses[i], 0), "\t\t", end="")
        print("\n====================================================================================================")
        print("Delta-v achieved =", round(self.calculate_delta_v(), 2))
        print("====================================================================================================")

    def create_rocket_object(self):
        for i in range(self.number_of_stages):
            self.rocket.add_stage(
                self.step_masses[i] * (1 - self.epsilons[i]),
                self.stage_masses[i],
                self.step_masses[i],
                self.Isps[i]
            )


'''