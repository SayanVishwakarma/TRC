from __init__ import *
class engine:
    def __init__(self,engine_name,override_inputs=True):
        self.engine_name=engine_name
        self.input_file="input files/engines/"+engine_name+".json"
        self.output_file="output files/engines/"+engine_name+".txt"
        with open(self.input_file) as file:
            self.data = json.load(file)
        self.read_fuel_data()

        if override_inputs:
            self.override_inputs()
        self.R_pb=[0.1,10]
        self.R_pbf=1
        self.R_pbo=1
        self.R=self.data["oxidiser by fuel"]

        self.m_z=self.data["m_z"]
        self.A_t=self.data["A_t"]
        self.p_f1=self.data["fuel tank pressure"]
        self.p_o1=self.data["oxidiser tank pressure"]
        self.p_c=self.m_z*self.cstar(self.R)/self.A_t

        self.eta_tf=self.data["fuel rich turbine efficiency"]
        self.eta_to=self.data["oxidiser rich turbine efficiency"]
        self.eta_pf=self.data["fuel pump efficiency"]
        self.eta_po=self.data["oxidiser pump efficiency"]
        
        self.rho_o=self.data["oxidiser density"]
        self.rho_f=self.data["fuel density"]
        self.p_pbf=self.data["fuel rich preburner pressure"]
        
        
        self.R_0=8.314
        
        self.ho=1e-1
        self.hf=1e-1

        self.C_ff=0.14       # Pressure drop loss coefficient between the fuel pump and the fuel-rich preburner
        self.C_oo=0.14       # Pressure drop loss coefficient between the oxidizer pump and the oxidizer-rich preburner
        self.C_f1=0.63       # Pressure drop loss coefficient between the fuel pump and fuel injector of the fuel-rich preburner
        self.C_o1=0       # Pressure drop loss coefficient between the oxidizer pump and oxidizer injector of the oxidizer-rich preburner
        self.C_f2=0       # Pressure drop loss coefficient between the fuel-rich preburner and the fuel-rich turbine
        self.C_o2=0       # Pressure drop loss coefficient between the oxidizer-rich preburner and the oxidizer-rich turbine
        self.C_f3=0.14       # Pressure drop loss coefficient between the fuel-rich turbine and the combustion chamber
        self.C_o3=0.43       # Pressure drop loss coefficient between the oxidizer-rich turbine and the combustion chamber
        #self.C_f4=0        Pressure drop loss coefficient between the fuel-rich preburner and the oxidizer-rich preburner
        #self.C_o4=0        Pressure drop loss coefficient between the oxidizer-rich preburner and the fuel-rich preburner
        self.C_fo=0.14       # Pressure drop loss coefficient between the oxidizer pump and the fuel-rich preburner
        self.C_of=0.45       # Pressure drop loss coefficient between the fuel pump and the oxidizer-rich preburner

        self.p_pbo=(1+self.C_fo)*self.p_pbf/(1+self.C_oo)

        self.max_residual=0

        self.symbols = ["D_E", "Isp", "Thrust", "p_c", "m_f", "m_o", "m_tf", "m_to", "m_ff", "m_fo", "m_of", "m_oo", "p_f2", "p_o2", "p_f3", "p_o3", "p_f4", "p_o4", "pi_tf", "pi_to", "delta_p_f", "delta_p_o", "P_f", "P_o", "g_f", "g_o", "gamma_pbf", "gamma_pbo", "M_pbf", "M_pbo", "T_pbf", "T_pbo", "R_pbf", "R_pbo"]

    def override_inputs(self):
        self.data["m_z"]=self.data["Thrust"]/self.isp(self.data["oxidiser by fuel"])
        self.data["A_t"]=self.data["m_z"]*self.cstar(self.data["oxidiser by fuel"])/self.data["p_c"]

    
    def read_fuel_data(self):
        self.fuel_data=np.genfromtxt("data files/LOXMETHANE cleaned.txt")[1:].T
        self.fuel_data_heads=["O/F","temp","isp","mw","cstar","gamma"]
           
    def temperature(self,R):
        ratios=self.fuel_data[0]
        ans=self.fuel_data[1]
        if ratios[0]<=R:
            return ans[0]
        for i in range(1,len(ratios)-1):
            if ratios[i]<=R:
                return (ans[i+1]-ans[i])/(ratios[i+1]-ratios[i])*(R-ratios[i])+ans[i]
        return ans[len(ratios)-1]  

    def isp(self,R):
        ratios=self.fuel_data[0]
        ans=self.fuel_data[2]
        if ratios[0]<=R:
            return ans[0]
        for i in range(1,len(ratios)-1):
            if ratios[i]<=R:
                return (ans[i+1]-ans[i])/(ratios[i+1]-ratios[i])*(R-ratios[i])+ans[i]
        return ans[len(ratios)-1]    
      
    def mol_w(self,R):
        ratios=self.fuel_data[0]
        ans=self.fuel_data[3]/1000
        if ratios[0]<=R:
            return ans[0]
        for i in range(1,len(ratios)-1):
            if ratios[i]<=R:
                return (ans[i+1]-ans[i])/(ratios[i+1]-ratios[i])*(R-ratios[i])+ans[i]
        return ans[len(ratios)-1]   

    def cstar(self,R):
        ratios=self.fuel_data[0]
        ans=self.fuel_data[4]
        if ratios[0]<=R:
            return ans[0]
        for i in range(1,len(ratios)-1):
            if ratios[i]<=R:
                return (ans[i+1]-ans[i])/(ratios[i+1]-ratios[i])*(R-ratios[i])+ans[i]
        return ans[len(ratios)-1]    
     
    def gamma(self,R):
        ratios=self.fuel_data[0]
        ans=self.fuel_data[5]
        if ratios[0]<=R:
            return ans[0]
        for i in range(1,len(ratios)-1):
            if ratios[i]<=R:
                return (ans[i+1]-ans[i])/(ratios[i+1]-ratios[i])*(R-ratios[i])+ans[i]
        return ans[len(ratios)-1]    
    
    def gas_constant(self,R):
        ratios=self.fuel_data[0]
        ans=self.fuel_data[3]/1000
        if ratios[0]<=R:
            return self.R_0/ans[0]
        for i in range(1,len(ratios)-1):
            if ratios[i]<=R:
                return self.R_0/((ans[i+1]-ans[i])/(ratios[i+1]-ratios[i])*(R-ratios[i])+ans[i])
        return self.R_0/ans[len(ratios)-1]
    
    def exit_diameter(self):

        p_p0=1e5/self.p_c
        gamma=self.gamma(self.R)
        # Define the isentropic pressure ratio function
        def pressure_ratio(M):
            return (1 + 0.5*(gamma -1)*M**2)**(-gamma/(gamma -1))

        # Define function whose root we want: pressure_ratio(M) - p_p0 = 0
        def f(M):
            return pressure_ratio(M) - p_p0

        # Solve for subsonic Mach (M < 1)
        #sol_sub = root_scalar(f, bracket=[1e-5, 1])
        #if not sol_sub.converged:
        #    raise ValueError("Subsonic Mach solving failed")
        #M_sub = sol_sub.root

        # Solve for supersonic Mach (M >1)
        sol_sup = root_scalar(f, bracket=[1.0001, 20])
        if not sol_sup.converged:
            raise ValueError("Supersonic Mach solving failed")
        M_sup = sol_sup.root

        # Area-Mach relation
        def area_ratio(M):
            return (1/M) * ( (2/(gamma+1)*(1 + 0.5*(gamma-1)*M**2)) ) ** ( (gamma+1)/(2*(gamma-1)) )

        #A_Astar_sub = area_ratio(M_sub)
        A_Astar_sup = area_ratio(M_sup)
        A_e=A_Astar_sup*self.A_t
        d_e=(A_e/pi)**0.5*2
        return d_e


    '''
    def c_star(self):
        #self.c_star=((self.data['R_c']*self.data['T_c']/self.data['gamma_c'])**0.5)/(2/(self.data['gamma_c']+1))**((self.data['gamma_c']+1)/2/(self.data['gamma_c']-1))
        return ((self.gas_constant(self.R)*self.temperature(self.R)/self.gamma(self.R))**0.5)/(2/(self.gamma(self.R)+1))**((self.gamma(self.R)+1)/2/(self.gamma(self.R)-1)) 
    
    def delta_p_f(self,R_pbo,R_pbf):
        return (R_pbo-self.R)*(1+R_pbf)/(R_pbo-R_pbf)*(self.eta_tf*self.rho_f*self.eta_pf)*(self.gamma_pbf/(self.gamma_pbf-1))*(self.R_0/self.M_pbf)*self.T_pbf*(1-(1/self.data['pi_tf'])**(self.gamma_pbf-1/self.gamma_pbf))
    
    def delta_p_o(self,R_pbo,R_pbf):
        return (self.R-R_pbf)*(1+R_pbf)/(R_pbo-R_pbf)*(self.eta_to*self.rho_o*self.eta_po)*(self.gamma_pbo/(self.gamma_pbo-1))*(self.R_0/self.M_pbo)*self.T_pbo*(1-(1/self.data['pi_to'])**(self.gamma_pbo-1/self.gamma_pbo))
    '''

    def delta_p_f(self,p_pbf):
        return (1+self.C_ff)*p_pbf+self.C_f1*self.p_c-self.p_f1
    
    def delta_p_o(self,p_pbf):
        return (1+self.C_fo)*p_pbf-self.p_o1

    def pi_tf(self):
        return (self.p_pbf-self.C_f2*self.m_z*self.cstar(self.R)/self.A_t)/(self.m_z*self.cstar(self.R_pbf)/self.A_t*(1+self.C_f3))    

    def pi_to(self):
        return (self.p_pbo-self.C_o2*self.m_z*self.cstar(self.R)/self.A_t)/(self.m_z*self.cstar(self.R_pbf)/self.A_t*(1+self.C_o3))
    
    def Gf(self,R_pbo):
        gamma_pbf=self.gamma(self.R_pbf)
        gamma_pbo=self.gamma(R_pbo)
        M_pbf=self.mol_w(self.R_pbf)
        M_pbo=self.mol_w(R_pbo)
        T_pbf=self.temperature(self.R_pbf)
        T_pbo=self.temperature(R_pbo)
        #print(self.delta_p_f(self.p_pbf_k)/1e5,self.delta_p_o(self.p_pbf_k)/1e5,gamma_pbf,gamma_pbo,M_pbf,M_pbo,T_pbf,T_pbo)
        #c1=self.delta_p_f(R_pbo,self.R_pb[0])/self.eta_tf/self.rho_f/self.eta_pf
        c1=self.delta_p_f(self.p_pbf)/self.eta_tf/self.rho_f/self.eta_pf
        #c2=self.delta_p_o(R_pbo,self.R_pb[0])/self.eta_to/self.rho_o/self.eta_po
        c2=self.delta_p_o(self.p_pbf)/self.eta_to/self.rho_o/self.eta_po
        w_f=gamma_pbf/(gamma_pbf-1)*(self.R_0/M_pbf)*T_pbf*(1-(1/self.pi_tf())**(gamma_pbf-1/gamma_pbf))
        w_o=gamma_pbo/(gamma_pbo-1)*(self.R_0/M_pbo)*T_pbo*(1-(1/self.pi_to())**(gamma_pbo-1/gamma_pbo))
        return (self.R*c2*R_pbo/(1+R_pbo)-self.R*w_o)/(self.R*c2/(1+R_pbo)-w_o)
    
    def Go(self,R_pbf):        
        gamma_pbf=self.gamma(R_pbf)
        gamma_pbo=self.gamma(self.R_pbo)
        M_pbf=self.mol_w(R_pbf)
        M_pbo=self.mol_w(self.R_pbo)
        T_pbf=self.temperature(R_pbf)
        T_pbo=self.temperature(self.R_pbo)
        #c1=self.delta_p_f(R_pbo,self.R_pb[0])/self.eta_tf/self.rho_f/self.eta_pf
        c1=self.delta_p_f(self.p_pbf)/self.eta_tf/self.rho_f/self.eta_pf
        #c2=self.delta_p_o(R_pbo,self.R_pb[0])/self.eta_to/self.rho_o/self.eta_po
        c2=self.delta_p_o(self.p_pbf)/self.eta_to/self.rho_o/self.eta_po
        w_f=gamma_pbf/(gamma_pbf-1)*(self.R_0/M_pbf)*T_pbf*(1-(1/self.pi_tf())**(gamma_pbf-1/gamma_pbf))
        w_o=gamma_pbo/(gamma_pbo-1)*(self.R_0/M_pbo)*T_pbo*(1-(1/self.pi_to())**(gamma_pbo-1/gamma_pbo))
        return (-c1*R_pbf/(1+R_pbf)+self.R*w_f)/(-c1/(1+R_pbf)+w_f)
    
    def objective_function(self,x):
        self.R_pbf=x[0]
        self.R_pbo=x[1]
        errs=self.calculate_all_parameters(finding_eq=True,print_residuals=False)
        return errs    
    
    def inner_convergence2(self):
        x=self.R_pb
        ans=root(self.objective_function,x)
        print(f"Equilibrium Ratios are {ans.x}")
        self.R_pb=ans.x
        self.R_pbf=ans.x[0]
        self.R_pbo=ans.x[1]


    def inner_convergence(self):
        while True:
            Fu=np.array([self.R_pbf-self.Gf(self.R_pbo),
                        self.R_pbo-self.Go(self.R_pbf)])
            J=np.array([
                        [1, (-self.Gf(self.R_pbo+self.ho)+self.Gf(self.R_pbo))/self.ho],
                        [(-self.Go(self.R_pbf+self.hf)+self.Go(self.R_pbf))/self.hf, 1]
                        ])
            delta_R_pb=np.linalg.solve(J,-Fu)
            if np.linalg.norm(delta_R_pb/self.R_pb)<1e-9:
                print(f"{self.p_pbf} has converged")
                break
            self.R_pbf=self.R_pbf+delta_R_pb[0]#*self.hf
            self.R_pbo=self.R_pbo+delta_R_pb[1]#*self.ho
            self.R_pb=[self.R_pbf,self.R_pbo]
            #print(self.R_pb,np.linalg.norm(delta_R_pb/self.R_pb))
        #print(self.R_pbf,self.R_pbo)
        self.R_pbf=0.17
        self.R_pbo=58

    def calculate_all_parameters(self,check_residuals=True,finding_eq=False,print_residuals=False):
        self.max_residual=0
        self.parameters={}
        self.parameters["Exit Diameter"]=self.exit_diameter()#((self.A_t*29.5)/np.pi)**0.5*2
        self.parameters["Isp"]=Isp=self.isp(self.R)
        self.parameters["Thrust"]=Thrust=self.m_z*self.parameters["Isp"]
        self.parameters["chamber pressure"]=p_c=self.m_z*self.cstar(self.R)/self.A_t
        self.parameters["fuel flow rate"]=m_f=self.m_z/(1+self.R)
        self.parameters["oxidiser flow rate"]=m_o=self.m_z/(1+self.R)*self.R
        self.parameters["fuel-rich turbine gas flow rate"]=m_tf=m_f*((self.R_pbo-self.R)*(1+self.R_pbf))/(self.R_pbo-self.R_pbf)
        self.parameters["oxidiser-rich turbine gas flow rate"]=m_to=m_o*((self.R-self.R_pbf) * (1 + self.R_pbo)) / ((self.R_pbo - self.R_pbf) * self.R)
        self.parameters["fuel flow rate of the fuel-rich preburner"]=m_ff=m_tf / (1 + self.R_pbf)
        self.parameters["oxidiser flow rate of the fuel-rich preburner"]=m_fo=m_tf*self.R_pbf / (1 + self.R_pbf)
        self.parameters["fuel flow rate of the oxidiser-rich preburner"]=m_of=m_to / (1 + self.R_pbo)
        self.parameters["oxidiser flow rate of the oxidiser-rich preburner"]=m_oo=m_to*self.R_pbo / (1 + self.R_pbo)
        self.parameters["outlet pressure of the fuel pump"]=p_f2=(self.p_f1 + self.delta_p_f(self.p_pbf))
        self.parameters["outlet pressure of the oxidizer pump"]=p_o2=(self.p_o1 + self.delta_p_o(self.p_pbf))
        self.parameters["inlet pressure of the fuel-rich turbine"]=p_f3=(self.p_pbf - self.C_f2 * p_c)
        self.parameters["inlet pressure of the oxidizer-rich turbine"]=p_o3=(self.p_pbo - self.C_o2 * p_c)
        self.parameters["outlet pressure of the fuel-rich turbine"]=p_f4=p_c * (1 + self.C_f3)
        self.parameters["outlet pressure of the oxidizer-rich turbine"]=p_o4=p_c * (1 + self.C_o3)
        self.parameters["fuel-rich turbine pressure ratio"]=pi_tf=p_f3/p_f4
        self.parameters["oxidizer-rich turbine pressure ratio"]=pi_to=p_o3/p_o4
        self.parameters["fuel pump head"]=delta_p_f=self.delta_p_f(self.p_pbf)
        self.parameters["oxidiser pump head"]=delta_p_o=self.delta_p_o(self.p_pbf)
        self.parameters["fuel turbopump power"]=P_f=m_f*delta_p_f/self.rho_f/self.eta_pf
        self.parameters["oxidiser turbopump power"]=P_o=m_o*delta_p_o/self.rho_o/self.eta_po
        self.parameters["fuel-rich turbine power"]=g_f=m_tf*self.eta_tf*self.gamma(self.R_pbf)/(self.gamma(self.R_pbf)-1)*(self.R_0/self.mol_w(self.R_pbf))*self.temperature(self.R_pbf)*(1-(1/pi_tf)**(self.gamma(self.R_pbf)-1/self.gamma(self.R_pbf)))
        self.parameters["oxidiser-rich turbine power"]=g_o=m_to*self.eta_to*self.gamma(self.R_pbo)/(self.gamma(self.R_pbo)-1)*(self.R_0/self.mol_w(self.R_pbo))*self.temperature(self.R_pbo)*(1-(1/pi_to)**(self.gamma(self.R_pbo)-1/self.gamma(self.R_pbo)))
        self.parameters["gamma in fuel rich preburner"]=gamma_pbf=self.gamma(self.R_pbf)
        self.parameters["gamma in oxidiser rich preburner"]=gamma_pbo=self.gamma(self.R_pbo)
        self.parameters["molecular weight in fuel rich preburner"]=M_pbf=self.mol_w(self.R_pbf)*1e3
        self.parameters["molecular weight in oxidiser rich preburner"]=M_pbo=self.mol_w(self.R_pbo)*1e3
        self.parameters["temperature in fuel rich preburner"]=T_pbf=self.temperature(self.R_pbf)
        self.parameters["temperature in oxidiser rich preburner"]=T_pbo=self.temperature(self.R_pbo)
        self.parameters["fuel rich preburner O/F"]=R_pbf=self.R_pbf
        self.parameters["oxidiser rich preburner O/F"]=R_pbo=self.R_pbo

        if check_residuals:
            eqs=[]
            eqs2=[]
            eqs.append((self.m_z-Thrust/Isp)) #1
            eqs.append((self.m_z-m_f-m_o)) #2
            eqs.append((p_c-self.m_z/self.A_t*self.cstar(self.R))) #3
            #eqs.append((self.F-C_f*self.A_t*P_c)/a)
            eqs.append((m_f-self.m_z/(1+self.R))) #4
            eqs.append((m_f-m_ff-m_of)) #5
            eqs.append((m_o-self.R*m_f)) #6
            eqs.append((m_o-m_fo-m_oo)) #7
            eqs.append((m_ff - m_tf / (1 + R_pbf))) #8
            eqs.append((m_fo - R_pbf * m_tf / (1 + R_pbf))) #9
            eqs.append((m_of - m_to / (1 + R_pbo))) #10
            eqs.append((m_oo - R_pbo * m_to / (1 + R_pbo))) #11
            eqs.append(m_tf - m_f*((R_pbo - self.R) * (1 + R_pbf)) / (R_pbo - R_pbf)) #12
            eqs.append(m_to - m_o*((self.R - R_pbf) * (1 + R_pbo)) / ((R_pbo - R_pbf) * self.R)) #13
            eqs.append((p_f2 - (self.p_f1 + delta_p_f))) #14
            eqs.append((p_o2 - (self.p_o1 + delta_p_o))) #15
            eqs.append((self.p_pbf - (p_f2 - self.C_f1 * p_c) / (1 + self.C_ff))) #16
            eqs.append((self.p_pbo - (p_o2 - self.C_o1 * p_c) / (1 + self.C_oo))) #17
            eqs.append((p_f3 - (self.p_pbf - self.C_f2 * p_c))) #18
            eqs.append((p_o3 - (self.p_pbo - self.C_o2 * p_c))) #19
            eqs.append((p_f4 - p_c * (1 + self.C_f3))) #20
            eqs.append((p_o4 - p_c * (1 + self.C_o3))) #21
            eqs.append(pi_tf - p_f3 / p_f4) #22
            eqs.append(pi_to - p_o3 / p_o4) #23
            eqs.append((delta_p_f-((1+self.C_ff)*self.p_pbf+self.C_f1*p_c-self.p_f1))) #24
            eqs.append((delta_p_o-((1+self.C_oo)*self.p_pbo-self.p_o1))) #25
            eqs.append(pi_tf-(self.p_pbf-self.C_f2*p_c)/(p_c*(1+self.C_f3)))#26
            eqs.append(pi_to-(self.p_pbo-self.C_o2*p_c)/(p_c*(1+self.C_o3))) #27
            eqs.append((p_o2-self.p_pbf*(1+self.C_fo))) #28
            eqs.append((p_o2-self.p_pbo*(1+self.C_oo)))#29
            eqs.append((p_f2-self.p_pbo*(1+self.C_of))/p_f2) #30
            eqs.append((p_f2-(self.p_pbf*(1+self.C_ff)+self.C_f1*p_c))) #31
            eqs.append(P_f-g_f) #32
            #eqs.append(P_f-p_f)
            eqs.append(P_o-g_o) #33
            #eqs.append(P_o-p_o)
            eqs.append(g_f-m_tf*self.eta_tf*self.gamma(self.R_pbf)/(self.gamma(self.R_pbf)-1)*(self.R_0/self.mol_w(self.R_pbf))*self.temperature(self.R_pbf)*(1-(1/pi_tf)**(self.gamma(self.R_pbf)-1/self.gamma(self.R_pbf)))) #34
            eqs.append(g_o-m_to*self.eta_to*self.gamma(self.R_pbo)/(self.gamma(self.R_pbo)-1)*(self.R_0/self.mol_w(self.R_pbo))*self.temperature(self.R_pbo)*(1-(1/pi_to)**(self.gamma(self.R_pbo)-1/self.gamma(self.R_pbo)))) #35
            eqs.append((P_f-m_f*delta_p_f/self.rho_f/self.eta_pf)) #36
            eqs.append((P_o-m_o*delta_p_o/self.rho_o/self.eta_po)) #37
            #eqs.append((delta_p_f-(R_pbo-self.R)*(1+R_pbf)/(R_pbo-R_pbf)*(self.eta_tf*self.rho_f*self.eta_pf)*(self.gamma_pbf/(self.gamma_pbf-1))*(self.R_0/self.M_pbf)*self.T_pbf*(1-(1/pi_tf)**(self.gamma_pbf-1/self.gamma_pbf)))/a)
            #eqs.append((delta_p_o-(self.R-R_pbf)*(1+R_pbf)/(R_pbo-R_pbf)*(self.eta_to*self.rho_o*self.eta_po)*(self.gamma_pbo/(self.gamma_pbo-1))*(self.R_0/self.M_pbo)*self.T_pbo*(1-(1/pi_to)**(self.gamma_pbo-1/self.gamma_pbo)))/a)

            eqs2.append(self.m_z) #1
            eqs2.append(self.m_z) #2
            eqs2.append(self.m_z) #3
            eqs2.append(m_f) #4
            eqs2.append(m_f) #5
            eqs2.append(m_o) #6
            eqs2.append(m_o) #7
            eqs2.append(m_ff) #8
            eqs2.append(m_fo) #9
            eqs2.append(m_of) #10
            eqs2.append(m_oo) #11
            eqs2.append(m_tf) #12
            eqs2.append(m_to) #13
            eqs2.append(p_f2) #14
            eqs2.append(p_o2) #15
            eqs2.append(self.p_pbf) #16
            eqs2.append(self.p_pbo) #17
            eqs2.append(p_f3) #18
            eqs2.append(p_o3) #19
            eqs2.append(p_f4) #20
            eqs2.append(p_o4) #21
            eqs2.append(pi_tf) #22
            eqs2.append(pi_to) #23
            eqs2.append(delta_p_f) #24
            eqs2.append(delta_p_o) #25
            eqs2.append(pi_tf) #26
            eqs2.append(pi_to) #27
            eqs2.append(p_o2) #28
            eqs2.append(p_o2) #29
            eqs2.append(p_f2) #30
            eqs2.append(p_f2) #31
            eqs2.append(P_f) #32
            eqs2.append(P_o) #33
            eqs2.append(g_f) #34
            eqs2.append(g_o) #35
            eqs2.append(P_f) #36
            eqs2.append(P_o) #37

            #eqs.append((delta_p_f-(R_pbo-self.R)*(1+R_pbf)/(R_pbo-R_pbf)*(self.eta_tf*self.rho_f*self.eta_pf)*(self.gamma_pbf/(self.gamma_pbf-1))*(self.R_0/self.M_pbf)*self.T_pbf*(1-(1/pi_tf)**(self.gamma_pbf-1/self.gamma_pbf)))/a)
            #eqs.append((delta_p_o-(self.R-R_pbf)*(1+R_pbf)/(R_pbo-R_pbf)*(self.eta_to*self.rho_o*self.eta_po)*(self.gamma_pbo/(self.gamma_pbo-1))*(self.R_0/self.M_pbo)*self.T_pbo*(1-(1/pi_to)**(self.gamma_pbo-1/self.gamma_pbo)))/a)
            if print_residuals:
                print("Residuals:")

            for i in range(0,len(eqs)):
                if print_residuals:
                    print(f"{i+1} : {round(eqs[i]/eqs2[i]*100,3)}%")
                if np.abs((eqs[i]/eqs2[i]))>self.max_residual:
                    self.max_residual=np.abs((eqs[i]/eqs2[i]))*100
        
        if finding_eq:
            return eqs[31]/eqs2[31],eqs[32]/eqs2[32]

    def get_output(self,print_residual=True):
        sys.stdout = open(self.output_file, 'w')
        print(f"***********************************************  {self.engine_name} ENGINE  ***********************************************")
        print()
        print(f"--- INPUT PARAMETERS ---\n")
        """
        df=pd.DataFrame.from_dict(self.data,orient="index")
        print(df.to_markdown())
        """
        #print("="*70)
        for key in self.data.keys():
            if ("pressure" in key and not("ratio" in key)) or "head" in key:
                print(f"{key} = {round(self.data[key]/1e5,2)} bar")
            elif "m_z" in key:
                print(f"{key} = {round(self.data[key],2)} kg/s")
            elif "A_t" in key:
                print(f"{key} = {round(self.data[key],4)} m^2")
            elif not("Thrust" in key or "p_c" in key):
                print(f"{key} = {self.data[key]}")
            #if not("Thrust" in key or "p_c" in key):
                #print("="*70)
        print()
        print("--- OUTPUT PARAMETERS ---\n")
        i=0
        #print("="*70)
        for key in self.parameters.keys():
            if "Exit Diameter" in key:
                print(f"{key} ({self.symbols[i]}) = {round(self.parameters[key],3)} m")
            elif ("pressure" in key and not("ratio" in key)) or "head" in key:
                print(f"{key} ({self.symbols[i]}) = {round(self.parameters[key]/1e5,2)} bar")
            elif "Isp" in key:
                print(f"{key} ({self.symbols[i]}) = {round(self.parameters[key]/g,2)} s")
            elif "Thrust" in key:
                print(f"{key} ({self.symbols[i]}) = {round(self.parameters[key]/1e3,2)} kN")
            elif "rate" in key:
                print(f"{key} ({self.symbols[i]}) = {round(self.parameters[key],2)} kg/s")
            elif "power" in key:
                print(f"{key} ({self.symbols[i]}) = {round(self.parameters[key]/1e6,2)} MW")
            elif "temperature" in key:
                print(f"{key} ({self.symbols[i]}) = {round(self.parameters[key],2)} K")
            else:
                print(f"{key} ({self.symbols[i]}) ={round(self.parameters[key],3)}")
            i+=1
            ##print("="*70)

        if print_residual:
            print()
            print(f"--- MAX RESIDUAL = {round(self.max_residual,3)}% ---")
        sys.stdout.close()
        sys.stdout = sys.__stdout__


eng=engine("2400kn")
eng.inner_convergence2()
eng.calculate_all_parameters(print_residuals=True)
eng.get_output()
#print(eng.c_star())
#print(eng.cstar(eng.R))3