from __init__ import *

class engine_sizing:

    def __init__(self,input_path):
        with open(input_path) as file:
            self.data = json.load(file)
        self.c_star()
        self.R_pbf=0.2
        self.R_pbo=2

    def c_star(self):
        self.c_star=((self.data['R_c']*self.data['T_c']/self.data['gamma_c'])**0.5)/(2/(self.data['gamma_c']+1))**((self.data['gamma_c']+1)/2/(self.data['gamma_c']-1))
    
    def equations(self,x,scale_outputs=True):
        if scale_outputs:
            a=1e5
            a2=1e3
        else:
            a=1
            a2=1
        self.F=self.data["thrust"]
        self.A_t=self.data['A_t']
        self.Isp=self.data['Isp']
        self.p_f1=self.data["fuel tank pressure"]       # Inlet pressure of the fuel pump (Pa)
        self.p_o1=self.data["fuel tank pressure"]       # Inlet pressure of the oxidizer pump (Pa)
        self.R=self.data["oxidiser_by_fuel"]          # Mixture ratio of the rocket engine
        C_ff=0.14       # Pressure drop loss coefficient between the fuel pump and the fuel-rich preburner
        C_oo=0.14       # Pressure drop loss coefficient between the oxidizer pump and the oxidizer-rich preburner
        C_f1=0.63       # Pressure drop loss coefficient between the fuel pump and fuel injector of the fuel-rich preburner
        C_o1=0       # Pressure drop loss coefficient between the oxidizer pump and oxidizer injector of the oxidizer-rich preburner
        C_f2=0       # Pressure drop loss coefficient between the fuel-rich preburner and the fuel-rich turbine
        C_o2=0       # Pressure drop loss coefficient between the oxidizer-rich preburner and the oxidizer-rich turbine
        C_f3=0.14       # Pressure drop loss coefficient between the fuel-rich turbine and the combustion chamber
        C_o3=0.43       # Pressure drop loss coefficient between the oxidizer-rich turbine and the combustion chamber
        #C_f4=0        Pressure drop loss coefficient between the fuel-rich preburner and the oxidizer-rich preburner
        #C_o4=0        Pressure drop loss coefficient between the oxidizer-rich preburner and the fuel-rich preburner
        C_fo=0.14       # Pressure drop loss coefficient between the oxidizer pump and the fuel-rich preburner
        C_of=0.45       # Pressure drop loss coefficient between the fuel pump and the oxidizer-rich preburner
        self.eta_tf=0.9     # Fuel-rich turbine efficiency
        self.gamma_pbf=1.13  # Specific heat ratio of fuel-rich preburner gas
        self.R_0=8.314        # Gas constant (8.314 J/(mol·K))
        self.M_pbf=22      # Molecular weight of fuel-rich gas (g/mol)
        self.T_pbf=1500      # Chemical equilibrium temperature of the fuel-rich preburner (K)
        self.rho_f=400      # Fuel density (kg/m³)
        self.eta_pf=0.9     # Fuel pump efficiency
        self.rho_o=1000      # Oxidizer density (kg/m³)
        self.eta_po=0.9     # Oxidizer pump efficiency
        self.gamma_pbo=1.13  # Specific heat ratio of oxidizer-rich preburner gas
        self.eta_to=0.9      #oxidizer-rich turbine efficiency
        self.M_pbo=22        #molecular weight of oxidizer-rich gas (g/mol)
        self.T_pbo=1500      #chemical equilibrium temperature of the oxidizer-rich preburner (K)
        [m_z,       #total propellant flow rate (kg/s)
         P_c,       #chamber pressure (Pa)
         C_f,       #Isp/c*
         m_f,       #fuel pump flow rate (kg/s)
         m_o,       #oxidizer pump flow rate (kg/s)
         m_ff,      #fuel flow rate of the fuel-rich preburner (kg/s)
         m_tf,      #fuel-rich turbine gas flow rate (kg/s)
         m_fo,      #oxidizer flow rate of the fuel-rich preburner (kg/s)
         m_of,      #fuel flow rate of the oxidizer-rich preburner (kg/s)
         m_to,      #oxidizer-rich turbine gas flow rate (kg/s)
         m_oo,      #oxidizer flow rate of the oxidizer-rich preburner (kg/s)
         p_f2,      #outlet pressure of the fuel pump (Pa)
         delta_p_f, #fuel pump head (Pa)
         p_o2,      #outlet pressure of the oxidizer pump (Pa)
         delta_p_o, #oxidizer pump head (Pa)
         p_pbf,     #fuel-rich preburner pressure (Pa), isi me bakchodi h hkhjgjlfhfdytkdhkgftydytsdtyhgfyytcyufytdytfudfyytfytdtykcytdtxgfytdt
         #R_pbo,     mixture ratio of the oxidizer-rich preburner
         #R_pbf,     mixture ratio of the fuel-rich preburner
         p_c,       #chamber pressure (Pa)
         p_pbo,     #oxidizer-rich preburner pressure (Pa)
         p_f3,      #inlet pressure of the fuel-rich turbine (Pa)
         p_o3,      #inlet pressure of the oxidizer-rich turbine (Pa)
         p_f4,      #outlet pressure of the fuel-rich turbine (Pa)
         p_o4,      #outlet pressure of the oxidizer-rich turbine (Pa)
         pi_tf,     #fuel-rich turbine pressure ratio
         pi_to,     #oxidizer-rich turbine pressure ratio
         P_f,       #fuel turbopump power (W)
         g_f,       #fuel-rich turbine power (W)
         p_f,       #fuel pump power (W)
         P_o,       #oxidizer turbopump power (W)
         g_o,       #oxidizer-rich turbine power (W)
         p_o,       #oxidizer pump power (W)
         dummy1,    #dummy variable, ignore
         dummy2,    #dummy variable, ignore
         dummy3,    #dummy variable, ignore
         dummy4,    #dummy variable, ignore
         dummy5,    #dummy variable, ignore
         dummy6,    #dummy variable, ignore
         dummy7,    #dummy variable, ignore
         dummy8,    #dummy variable, ignore
         #dummy9,    dummy variable, ignore
         #dummy10,   dummy variable, ignore
         #dummy11,    dummy variable, ignore
         #dummy12,    dummy variable, ignore         
         ]=x
        R_pbf=self.R_pbf
        R_pbo=self.R_pbo
        eqs=[]
        eqs.append((m_z-self.F/self.Isp/g)*a2)
        eqs.append((m_z-m_f-m_o)*a2)
        eqs.append((P_c-m_z/self.A_t*self.c_star)/a)
        eqs.append((self.F-C_f*self.A_t*P_c)/a)
        eqs.append((m_f-m_z/(1+self.R))*a2)
        eqs.append((m_f-m_ff-m_of)*a2)
        eqs.append((m_o-self.R*m_f)*a2)
        eqs.append((m_o-m_fo-m_oo)*a2)
        eqs.append((m_ff - m_tf / (1 + R_pbf))*a2)
        eqs.append((m_fo - R_pbf * m_tf / (1 + R_pbf))*a2)
        eqs.append((m_of - m_to / (1 + R_pbo))*a2)
        eqs.append((m_oo - R_pbo * m_to / (1 + R_pbo))*a2)
        eqs.append(m_tf / m_f - ((R_pbo - self.R) * (1 + R_pbf)) / (R_pbo - R_pbf))
        eqs.append(m_to / m_o - ((self.R - R_pbf) * (1 + R_pbo)) / ((R_pbo - R_pbf) * self.R))
        eqs.append((p_f2 - (self.p_f1 + delta_p_f))/a)
        eqs.append((p_o2 - (self.p_o1 + delta_p_o))/a)
        eqs.append((p_pbf - (p_f2 - C_f1 * p_c) / (1 + C_ff))/a)
        eqs.append((p_pbo - (p_o2 - C_o1 * p_c) / (1 + C_oo))/a)
        eqs.append((p_f3 - (p_pbf - C_f2 * p_c))/a)
        eqs.append((p_o3 - (p_pbo - C_o2 * p_c))/a)
        eqs.append((p_f4 - p_c * (1 + C_f3))/a)
        eqs.append((p_o4 - p_c * (1 + C_o3))/a)
        eqs.append(pi_tf - p_f3 / p_f4)
        eqs.append(pi_to - p_o3 / p_o4)
        eqs.append((delta_p_f-(1+C_ff)*p_pbf+C_f1*p_c-self.p_f1)/a)
        eqs.append((delta_p_o-(1+C_oo)*p_pbo-self.p_o1)/a)
        #eqs.append(pi_tf-(p_pbf-C_f2*p_c)/(p_c*(1+C_f3)))
        #eqs.append(pi_to-(p_pbo-C_o2*p_c)/(p_c*(1+C_o3)))
        #eqs.append((p_o2-p_pbf*(1+C_fo))/a)
        #eqs.append((p_o2-p_pbo*(1+C_oo)/a))
        eqs.append((p_f2-p_pbo*(1+C_of))/a)
        eqs.append((p_f2-p_pbf*(1+C_ff)-C_f1*p_c)/a)
        eqs.append(P_f-p_f)
        eqs.append(p_f-g_f)
        eqs.append(P_o-p_o)
        eqs.append(p_o-g_o)
        eqs.append(g_f-m_tf*self.eta_tf*self.gamma_pbf/(self.gamma_pbf-1)*(self.R_0/self.M_pbf)*self.T_pbf*(1-(1/pi_tf)**(self.gamma_pbf-1/self.gamma_pbf)))
        eqs.append(g_o-m_to*self.eta_to*self.gamma_pbo/(self.gamma_pbo-1)*(self.R_0/self.M_pbo)*self.T_pbf*(1-(1/pi_tf)**(self.gamma_pbf-1/self.gamma_pbf)))
        eqs.append(p_f-m_f*delta_p_f/self.rho_f/self.eta_pf)
        eqs.append(p_o-m_o*delta_p_o/self.rho_o/self.eta_po)
        eqs.append((delta_p_f-(R_pbo-self.R)*(1+R_pbf)/(R_pbo-R_pbf)*(self.eta_tf*self.rho_f*self.eta_pf)*(self.gamma_pbf/(self.gamma_pbf-1))*(self.R_0/self.M_pbf)*self.T_pbf*(1-(1/pi_tf)**(self.gamma_pbf-1/self.gamma_pbf)))/a)
        eqs.append((delta_p_o-(self.R-R_pbf)*(1+R_pbf)/(R_pbo-R_pbf)*(self.eta_to*self.rho_o*self.eta_po)*(self.gamma_pbo/(self.gamma_pbo-1))*(self.R_0/self.M_pbo)*self.T_pbo*(1-(1/pi_to)**(self.gamma_pbo-1/self.gamma_pbo)))/a)
        return eqs

    def initialise_x(self):
        return [600,       #total propellant flow rate (kg/s)
                200e5,       #chamber pressure (Pa)
                1.3,       #Isp/c*
                200,       #fuel pump flow rate (kg/s)
                400,       #oxidizer pump flow rate (kg/s)
                100,      #fuel flow rate of the fuel-rich preburner (kg/s)
                100,      #fuel-rich turbine gas flow rate (kg/s)
                100,      #oxidizer flow rate of the fuel-rich preburner (kg/s)
                100,      #fuel flow rate of the oxidizer-rich preburner (kg/s)
                100,      #oxidizer-rich turbine gas flow rate (kg/s)
                100,      #oxidizer flow rate of the oxidizer-rich preburner (kg/s)
                20e5,      #outlet pressure of the fuel pump (Pa)
                5e5, #fuel pump head (Pa)
                60e6,      #outlet pressure of the oxidizer pump (Pa)
                50e5, #oxidizer pump head (Pa)
                50e5,     #fuel-rich preburner pressure (Pa)
                #2,     mixture ratio of the oxidizer-rich preburner
                #2,     mixture ratio of the fuel-rich preburner
                200e5,       #chamber pressure (Pa)
                50e5,     #oxidizer-rich preburner pressure (Pa)
                50e5,      #inlet pressure of the fuel-rich turbine (Pa)
                50e5,      #inlet pressure of the oxidizer-rich turbine (Pa)
                50e5,      #outlet pressure of the fuel-rich turbine (Pa)
                50e5,      #outlet pressure of the oxidizer-rich turbine (Pa)
                4,     #fuel-rich turbine pressure ratio
                5,     #oxidizer-rich turbine pressure ratio
                200e3,       #fuel turbopump power (W)
                200e3,       #fuel-rich turbine power (W)
                200e3,       #fuel pump power (W)
                200e3,       #oxidizer turbopump power (W)
                200e3,       #oxidizer-rich turbine power (W)
                200e3,       #oxidizer pump power (W)
                0,    #dummy variable, ignore
                0,    #dummy variable, ignore
                0,    #dummy variable, ignore
                0,    #dummy variable, ignore
                0,    #dummy variable, ignore
                0,    #dummy variable, ignore
                0,    #dummy variable, ignore
                0,    #dummy variable, ignore
                #0,    #dummy variable, ignore
                #0,    #dummy variable, ignore
                #0,    dummy variable, ignore
                #0,    dummy variable, ignore         
                ]
    
    def delta_p_f(self,R_pbo,R_pbf):
        return (R_pbo-self.R)*(1+R_pbf)/(R_pbo-R_pbf)*(self.eta_tf*self.rho_f*self.eta_pf)*(self.gamma_pbf/(self.gamma_pbf-1))*(self.R_0/self.M_pbf)*self.T_pbf*(1-(1/pi_tf)**(self.gamma_pbf-1/self.gamma_pbf))

    def delta_p_o(self,R_pbo,R_pbf):
        return (self.R-R_pbf)*(1+R_pbf)/(R_pbo-R_pbf)*(self.eta_to*self.rho_o*self.eta_po)*(self.gamma_pbo/(self.gamma_pbo-1))*(self.R_0/self.M_pbo)*self.T_pbo*(1-(1/pi_to)**(self.gamma_pbo-1/self.gamma_pbo))
    def Gf(self,R_pbo,R_pbf):
        c1=self.delta_p_f(R_pbo,R_pbf)/self.eta_tf/self.rho_f/self.eta_pf
        c2=self.parameters[14]/self.eta_to/self.rho_o/self.eta_po
        w_f=self.gamma_pbf/(self.gamma_pbf-1)*(self.R_0/self.M_pbf)*self.T_pbf*(1-(1/self.parameters[24])**(self.gamma_pbf-1/self.gamma_pbf))
        w_o=self.gamma_pbo/(self.gamma_pbo-1)*(self.R_0/self.M_pbo)*self.T_pbo*(1-(1/self.parameters[25])**(self.gamma_pbo-1/self.gamma_pbo))
        return (self.R*c2*R_pbo/(1+R_pbo)-self.R*w_o)/(self.R*c2/(1+R_pbo)-w_o)
    
    def Go(self,R_pbo,R_pbf):
        c1=self.parameters[12]/self.eta_tf/self.rho_f/self.eta_pf
        c2=self.parameters[14]/self.eta_to/self.rho_o/self.eta_po
        w_f=self.gamma_pbf/(self.gamma_pbf-1)*(self.R_0/self.M_pbf)*self.T_pbf*(1-(1/self.parameters[24])**(self.gamma_pbf-1/self.gamma_pbf))
        w_o=self.gamma_pbo/(self.gamma_pbo-1)*(self.R_0/self.M_pbo)*self.T_pbo*(1-(1/self.parameters[25])**(self.gamma_pbo-1/self.gamma_pbo))
        return (-c1*R_pbo/(1+R_pbo)+self.R*w_f)/(-c1/(1+R_pbo)+w_f)

    def get_roots(self):
        x=self.initialise_x()
        ans=root(self.equations,np.array(x))
        #print("Solved Parameters")
        #self.print_parameters(ans.x)
        print(ans.message)
        print("Residuals")
        print(self.equations(ans.x))

    def converge_results(self,ho,hf):
        self.parameters=self.initialise_x()
        print(f"R_pbf={self.R_pbf}, R_pbo={self.R_pbo}")
        while True:    
            ans=root(self.equations,self.parameters)
            R_pb=np.array([self.R_pbf,self.R_pbo])
            self.parameters=ans.x
            print(ans.message)
            Fu=np.array([self.R_pbf-self.Gf(self.R_pbo),self.R_pbo-self.Go(self.R_pbf)])
            J=np.array([
                        [1, (-self.Gf(self.R_pbo+ho)+self.Gf(self.R_pbo))/ho],
                        [(-self.Go(self.R_pbf+hf)+self.Go(self.R_pbf))/hf, 1]
                        ])
            delta_R_pb=np.linalg.solve(J,-Fu)
            if np.linalg.norm(delta_R_pb/R_pb)<1e-2:
                break
            self.R_pbf=self.R_pbf+delta_R_pb[0]*hf
            self.R_pbo=self.R_pbo+delta_R_pb[1]*ho
            print(f"R_pbf={self.R_pbf}, R_pbo={self.R_pbo}, Error={np.linalg.norm(delta_R_pb/R_pb)}")

    
    def print_parameters(self):
        x=self.parameters
        labels=['total propellant flow rate (kg/s)',
                'chamber pressure (Pa)',
                'Isp/c*',
                'fuel pump flow rate (kg/s)',
                'oxidizer pump flow rate (kg/s)',
                'fuel flow rate of the fuel-rich preburner (kg/s)',
                'fuel-rich turbine gas flow rate (kg/s)',
                'oxidizer flow rate of the fuel-rich preburner (kg/s)',
                'fuel flow rate of the oxidizer-rich preburner (kg/s)',
                'oxidizer-rich turbine gas flow rate (kg/s)',
                'oxidizer flow rate of the oxidizer-rich preburner (kg/s)',
                'outlet pressure of the fuel pump (Pa)',
                'fuel pump head (Pa)',
                'outlet pressure of the oxidizer pump (Pa)',
                'oxidizer pump head (Pa)',
                'fuel-rich preburner pressure (Pa)',
                'mixture ratio of the oxidizer-rich preburner',
                'mixture ratio of the fuel-rich preburner',
                'chamber pressure (Pa)',
                'oxidizer-rich preburner pressure (Pa)',
                'inlet pressure of the fuel-rich turbine (Pa)',
                'inlet pressure of the oxidizer-rich turbine (Pa)',
                'outlet pressure of the fuel-rich turbine (Pa)',
                'outlet pressure of the oxidizer-rich turbine (Pa)',
                'fuel-rich turbine pressure ratio',
                'oxidizer-rich turbine pressure ratio',
                'fuel turbopump power (W)',
                'fuel-rich turbine power (W)',
                'fuel pump power (W)',
                'oxidizer turbopump power (W)',
                'oxidizer-rich turbine power (W)',
                'oxidizer pump power (W)',
                'dummy variable, ignore',
                'dummy variable, ignore',
                'dummy variable, ignore',
                'dummy variable, ignore',
                'dummy variable, ignore',
                'dummy variable, ignore',
                'dummy variable, ignore',
                'dummy variable, ignore',
                'dummy variable, ignore',
                'dummy variable, ignore',
                ]
        for i in range(len(labels)-6):
            print(f"{labels[i]}={x[i]}")


eng=engine_sizing("input files/bheem.json")
eng.converge_results(0.01,0.01)
sys.stdout = open("output files/engines/bheem.txt", 'w')
eng.print_parameters()
print(eng.equations(eng.parameters,scale_outputs=False)/eng.parameters)
sys.stdout.close()
sys.stdout = sys.__stdout__


class engine_sizing2:

    def __init__(self,input_path):
        with open(input_path) as file:
            self.data = json.load(file)
        self.c_star()
        self.R_pbf=0.2
        self.R_pbo=2

    def c_star(self):
        self.c_star=((self.data['R_c']*self.data['T_c']/self.data['gamma_c'])**0.5)/(2/(self.data['gamma_c']+1))**((self.data['gamma_c']+1)/2/(self.data['gamma_c']-1))
    
    def equations(self,x,scale_outputs=True):
        if scale_outputs:
            a=1e5
            a2=1e3
        else:
            a=1
            a2=1

        [m_z,       #total propellant flow rate (kg/s)
         P_c,       #chamber pressure (Pa)
         #C_f,       #Isp/c*
         m_f,       #fuel pump flow rate (kg/s)
         m_o,       #oxidizer pump flow rate (kg/s)
         m_ff,      #fuel flow rate of the fuel-rich preburner (kg/s)
         m_fo,      #oxidizer flow rate of the fuel-rich preburner (kg/s)
         m_of,      #fuel flow rate of the oxidizer-rich preburner (kg/s)
         m_oo,      #oxidizer flow rate of the oxidizer-rich preburner (kg/s)
         m_tf,      #fuel-rich turbine gas flow rate (kg/s)
         m_to,      #oxidizer-rich turbine gas flow rate (kg/s)
         p_f2,      #outlet pressure of the fuel pump (Pa)
         p_o2,      #outlet pressure of the oxidizer pump (Pa)
         p_pbf,     #fuel-rich preburner pressure (Pa), isi me bakchodi h hkhjgjlfhfdytkdhkgftydytsdtyhgfyytcyufytdytfudfyytfytdtykcytdtxgfytdt
         p_pbo,     #oxidizer-rich preburner pressure (Pa)
         p_f3,      #inlet pressure of the fuel-rich turbine (Pa)
         p_o3,      #inlet pressure of the oxidizer-rich turbine (Pa)
         p_f4,      #outlet pressure of the fuel-rich turbine (Pa)
         p_o4,      #outlet pressure of the oxidizer-rich turbine (Pa)
         pi_tf,     #fuel-rich turbine pressure ratio
         pi_to,     #oxidizer-rich turbine pressure ratio
         delta_p_f, #fuel pump head (Pa)
         delta_p_o, #oxidizer pump head (Pa)
         #R_pbo,     mixture ratio of the oxidizer-rich preburner
         #R_pbf,     mixture ratio of the fuel-rich preburner
         P_f,       #fuel turbopump power (W)
         g_f,       #fuel-rich turbine power (W)
         p_f,       #fuel pump power (W)
         P_o,       #oxidizer turbopump power (W)
         g_o,       #oxidizer-rich turbine power (W)
         p_o,       #oxidizer pump power (W)
         dummy1,    #dummy variable, ignore
         dummy2,    #dummy variable, ignore
         dummy3,    #dummy variable, ignore
         dummy4,    #dummy variable, ignore
         dummy5,    #dummy variable, ignore
         dummy6,    #dummy variable, ignore
         dummy7,    #dummy variable, ignore
         dummy8,    #dummy variable, ignore
         dummy9,    #dummy variable, ignore
         dummy10,    #dummy variable, ignore
         dummy11,    #dummy variable, ignore
         dummy12,    #dummy variable, ignore         
         ]=x
        

        self.F=self.data["thrust"]
        self.A_t=self.data['A_t']
        self.Isp=self.data['Isp']
        self.p_f1=self.data["fuel tank pressure"]       # Inlet pressure of the fuel pump (Pa)
        self.p_o1=self.data["fuel tank pressure"]       # Inlet pressure of the oxidizer pump (Pa)
        self.R=self.data["oxidiser_by_fuel"]          # Mixture ratio of the rocket engine
        C_ff=0.14       # Pressure drop loss coefficient between the fuel pump and the fuel-rich preburner
        C_oo=0.14       # Pressure drop loss coefficient between the oxidizer pump and the oxidizer-rich preburner
        C_f1=0.63       # Pressure drop loss coefficient between the fuel pump and fuel injector of the fuel-rich preburner
        C_o1=0       # Pressure drop loss coefficient between the oxidizer pump and oxidizer injector of the oxidizer-rich preburner
        C_f2=0       # Pressure drop loss coefficient between the fuel-rich preburner and the fuel-rich turbine
        C_o2=0       # Pressure drop loss coefficient between the oxidizer-rich preburner and the oxidizer-rich turbine
        C_f3=0.14       # Pressure drop loss coefficient between the fuel-rich turbine and the combustion chamber
        C_o3=0.43       # Pressure drop loss coefficient between the oxidizer-rich turbine and the combustion chamber
        C_f4=0       # Pressure drop loss coefficient between the fuel-rich preburner and the oxidizer-rich preburner
        C_o4=0       # Pressure drop loss coefficient between the oxidizer-rich preburner and the fuel-rich preburner
        C_fo=0.14       # Pressure drop loss coefficient between the oxidizer pump and the fuel-rich preburner
        C_of=0.45       # Pressure drop loss coefficient between the fuel pump and the oxidizer-rich preburner
        self.eta_tf=0.9     # Fuel-rich turbine efficiency
        self.gamma_pbf=1.13  # Specific heat ratio of fuel-rich preburner gas
        self.R_0=8.314        # Gas constant (8.314 J/(mol·K))
        self.M_pbf=22      # Molecular weight of fuel-rich gas (g/mol)
        self.T_pbf=1500      # Chemical equilibrium temperature of the fuel-rich preburner (K)
        self.rho_f=400      # Fuel density (kg/m³)
        self.eta_pf=0.9     # Fuel pump efficiency
        self.rho_o=1000      # Oxidizer density (kg/m³)
        self.eta_po=0.9     # Oxidizer pump efficiency
        self.gamma_pbo=1.13  # Specific heat ratio of oxidizer-rich preburner gas
        self.eta_to=0.9      #oxidizer-rich turbine efficiency
        self.M_pbo=22        #molecular weight of oxidizer-rich gas (g/mol)
        self.T_pbo=1500      #chemical equilibrium temperature of the oxidizer-rich preburner (K)
        
        R_pbf=self.R_pbf
        R_pbo=self.R_pbo
        eqs=[]
        eqs.append((m_z-self.F/self.Isp/g)*a2)
        eqs.append((m_z-m_f-m_o)*a2)
        eqs.append((P_c-m_z/self.A_t*self.c_star)/a)
        #eqs.append((self.F-C_f*self.A_t*P_c)/a)
        eqs.append((m_f-m_z/(1+self.R))*a2)
        eqs.append((m_f-m_ff-m_of)*a2)
        eqs.append((m_o-self.R*m_f)*a2)
        eqs.append((m_o-m_fo-m_oo)*a2)
        eqs.append((m_ff - m_tf / (1 + R_pbf))*a2)
        eqs.append((m_fo - R_pbf * m_tf / (1 + R_pbf))*a2)
        eqs.append((m_of - m_to / (1 + R_pbo))*a2)
        eqs.append((m_oo - R_pbo * m_to / (1 + R_pbo))*a2)
        eqs.append(m_tf - m_f*((R_pbo - self.R) * (1 + R_pbf)) / (R_pbo - R_pbf))
        eqs.append(m_to - m_o*((self.R - R_pbf) * (1 + R_pbo)) / ((R_pbo - R_pbf) * self.R))
        eqs.append((p_f2 - (self.p_f1 + delta_p_f))/a)
        eqs.append((p_o2 - (self.p_o1 + delta_p_o))/a)
        eqs.append((p_pbf - (p_f2 - C_f1 * p_c) / (1 + C_ff))/a)
        eqs.append((p_pbo - (p_o2 - C_o1 * p_c) / (1 + C_oo))/a)
        eqs.append((p_f3 - (p_pbf - C_f2 * p_c))/a)
        eqs.append((p_o3 - (p_pbo - C_o2 * p_c))/a)
        eqs.append((p_f4 - p_c * (1 + C_f3))/a)
        eqs.append((p_o4 - p_c * (1 + C_o3))/a)
        eqs.append(pi_tf - p_f3 / p_f4)
        eqs.append(pi_to - p_o3 / p_o4)
        eqs.append((delta_p_f-(1+C_ff)*p_pbf+C_f1*p_c-self.p_f1)/a)
        eqs.append((delta_p_o-(1+C_oo)*p_pbo-self.p_o1)/a)
        eqs.append(pi_tf-(p_pbf-C_f3*p_c)/(p_c*(1+C_f4)))
        eqs.append(pi_to-(p_pbo-C_o3*p_c)/(p_c*(1+C_o4)))
        eqs.append((p_o2-p_pbf*(1+C_fo))/a)
        eqs.append((p_o2-p_pbo*(1+C_oo)/a))
        eqs.append((p_f2-p_pbo*(1+C_of))/a)
        eqs.append((p_f2-p_pbf*(1+C_ff)-C_f1*p_c)/a)
        eqs.append(P_f-g_f)
        eqs.append(P_f-p_f)
        eqs.append(P_o-g_o)
        eqs.append(P_o-p_o)
        eqs.append(g_f-m_tf*self.eta_tf*self.gamma_pbf/(self.gamma_pbf-1)*(self.R_0/self.M_pbf)*self.T_pbf*(1-(1/pi_tf)**(self.gamma_pbf-1/self.gamma_pbf)))
        eqs.append(g_o-m_to*self.eta_to*self.gamma_pbo/(self.gamma_pbo-1)*(self.R_0/self.M_pbo)*self.T_pbf*(1-(1/pi_tf)**(self.gamma_pbf-1/self.gamma_pbf)))
        eqs.append(p_f-m_f*delta_p_f/self.rho_f/self.eta_pf)
        eqs.append(p_o-m_o*delta_p_o/self.rho_o/self.eta_po)
        eqs.append((delta_p_f-(R_pbo-self.R)*(1+R_pbf)/(R_pbo-R_pbf)*(self.eta_tf*self.rho_f*self.eta_pf)*(self.gamma_pbf/(self.gamma_pbf-1))*(self.R_0/self.M_pbf)*self.T_pbf*(1-(1/pi_tf)**(self.gamma_pbf-1/self.gamma_pbf)))/a)
        eqs.append((delta_p_o-(self.R-R_pbf)*(1+R_pbf)/(R_pbo-R_pbf)*(self.eta_to*self.rho_o*self.eta_po)*(self.gamma_pbo/(self.gamma_pbo-1))*(self.R_0/self.M_pbo)*self.T_pbo*(1-(1/pi_to)**(self.gamma_pbo-1/self.gamma_pbo)))/a)
        return eqs

    def initialise_x(self):
        return [600,       #total propellant flow rate (kg/s)
                200e5,       #chamber pressure (Pa)
                1.3,       #Isp/c*
                200,       #fuel pump flow rate (kg/s)
                400,       #oxidizer pump flow rate (kg/s)
                100,      #fuel flow rate of the fuel-rich preburner (kg/s)
                100,      #fuel-rich turbine gas flow rate (kg/s)
                100,      #oxidizer flow rate of the fuel-rich preburner (kg/s)
                100,      #fuel flow rate of the oxidizer-rich preburner (kg/s)
                100,      #oxidizer-rich turbine gas flow rate (kg/s)
                100,      #oxidizer flow rate of the oxidizer-rich preburner (kg/s)
                20e5,      #outlet pressure of the fuel pump (Pa)
                5e5, #fuel pump head (Pa)
                60e6,      #outlet pressure of the oxidizer pump (Pa)
                50e5, #oxidizer pump head (Pa)
                50e5,     #fuel-rich preburner pressure (Pa)
                #2,     mixture ratio of the oxidizer-rich preburner
                #2,     mixture ratio of the fuel-rich preburner
                200e5,       #chamber pressure (Pa)
                50e5,     #oxidizer-rich preburner pressure (Pa)
                50e5,      #inlet pressure of the fuel-rich turbine (Pa)
                50e5,      #inlet pressure of the oxidizer-rich turbine (Pa)
                50e5,      #outlet pressure of the fuel-rich turbine (Pa)
                50e5,      #outlet pressure of the oxidizer-rich turbine (Pa)
                4,     #fuel-rich turbine pressure ratio
                5,     #oxidizer-rich turbine pressure ratio
                200e3,       #fuel turbopump power (W)
                200e3,       #fuel-rich turbine power (W)
                200e3,       #fuel pump power (W)
                200e3,       #oxidizer turbopump power (W)
                200e3,       #oxidizer-rich turbine power (W)
                200e3,       #oxidizer pump power (W)
                0,    #dummy variable, ignore
                0,    #dummy variable, ignore
                0,    #dummy variable, ignore
                0,    #dummy variable, ignore
                0,    #dummy variable, ignore
                0,    #dummy variable, ignore
                0,    #dummy variable, ignore
                0,    #dummy variable, ignore
                0,    #dummy variable, ignore
                0,    #dummy variable, ignore
                0,    #dummy variable, ignore
                0,    #dummy variable, ignore         
                ]
    
    def Gf(self,R_pbo):
        c1=self.parameters[12]/self.eta_tf/self.rho_f/self.eta_pf
        c2=self.parameters[14]/self.eta_to/self.rho_o/self.eta_po
        w_f=self.gamma_pbf/(self.gamma_pbf-1)*(self.R_0/self.M_pbf)*self.T_pbf*(1-(1/self.parameters[24])**(self.gamma_pbf-1/self.gamma_pbf))
        w_o=self.gamma_pbo/(self.gamma_pbo-1)*(self.R_0/self.M_pbo)*self.T_pbo*(1-(1/self.parameters[25])**(self.gamma_pbo-1/self.gamma_pbo))
        return (self.R*c2*R_pbo/(1+R_pbo)-self.R*w_o)/(self.R*c2/(1+R_pbo)-w_o)
    
    def Go(self,R_pbo):
        c1=self.parameters[12]/self.eta_tf/self.rho_f/self.eta_pf
        c2=self.parameters[14]/self.eta_to/self.rho_o/self.eta_po
        w_f=self.gamma_pbf/(self.gamma_pbf-1)*(self.R_0/self.M_pbf)*self.T_pbf*(1-(1/self.parameters[24])**(self.gamma_pbf-1/self.gamma_pbf))
        w_o=self.gamma_pbo/(self.gamma_pbo-1)*(self.R_0/self.M_pbo)*self.T_pbo*(1-(1/self.parameters[25])**(self.gamma_pbo-1/self.gamma_pbo))
        return (-c1*R_pbo/(1+R_pbo)+self.R*w_f)/(-c1/(1+R_pbo)+w_f)

    def get_roots(self):
        x=self.initialise_x()
        ans=root(self.equations,np.array(x))
        #print("Solved Parameters")
        #self.print_parameters(ans.x)
        print(ans.message)
        print("Residuals")
        print(self.equations(ans.x))

    def converge_results(self,ho,hf):
        self.parameters=self.initialise_x()
        print(f"R_pbf={self.R_pbf}, R_pbo={self.R_pbo}")
        while True:    
            ans=root(self.equations,self.parameters)
            R_pb=np.array([self.R_pbf,self.R_pbo])
            self.parameters=ans.x
            print(ans.message)
            Fu=np.array([self.R_pbf-self.Gf(self.R_pbo),self.R_pbo-self.Go(self.R_pbf)])
            J=np.array([
                        [1, (-self.Gf(self.R_pbo+ho)+self.Gf(self.R_pbo))/ho],
                        [(-self.Go(self.R_pbf+hf)+self.Go(self.R_pbf))/hf, 1]
                        ])
            delta_R_pb=np.linalg.solve(J,-Fu)
            if np.linalg.norm(delta_R_pb/R_pb)<1e-2:
                break
            self.R_pbf=self.R_pbf+delta_R_pb[0]*hf
            self.R_pbo=self.R_pbo+delta_R_pb[1]*ho
            print(f"R_pbf={self.R_pbf}, R_pbo={self.R_pbo}, Error={np.linalg.norm(delta_R_pb/R_pb)}")

    
    def print_parameters(self):
        x=self.parameters
        labels=['total propellant flow rate (kg/s)',
                'chamber pressure (Pa)',
                'Isp/c*',
                'fuel pump flow rate (kg/s)',
                'oxidizer pump flow rate (kg/s)',
                'fuel flow rate of the fuel-rich preburner (kg/s)',
                'fuel-rich turbine gas flow rate (kg/s)',
                'oxidizer flow rate of the fuel-rich preburner (kg/s)',
                'fuel flow rate of the oxidizer-rich preburner (kg/s)',
                'oxidizer-rich turbine gas flow rate (kg/s)',
                'oxidizer flow rate of the oxidizer-rich preburner (kg/s)',
                'outlet pressure of the fuel pump (Pa)',
                'fuel pump head (Pa)',
                'outlet pressure of the oxidizer pump (Pa)',
                'oxidizer pump head (Pa)',
                'fuel-rich preburner pressure (Pa)',
                'mixture ratio of the oxidizer-rich preburner',
                'mixture ratio of the fuel-rich preburner',
                'chamber pressure (Pa)',
                'oxidizer-rich preburner pressure (Pa)',
                'inlet pressure of the fuel-rich turbine (Pa)',
                'inlet pressure of the oxidizer-rich turbine (Pa)',
                'outlet pressure of the fuel-rich turbine (Pa)',
                'outlet pressure of the oxidizer-rich turbine (Pa)',
                'fuel-rich turbine pressure ratio',
                'oxidizer-rich turbine pressure ratio',
                'fuel turbopump power (W)',
                'fuel-rich turbine power (W)',
                'fuel pump power (W)',
                'oxidizer turbopump power (W)',
                'oxidizer-rich turbine power (W)',
                'oxidizer pump power (W)',
                'dummy variable, ignore',
                'dummy variable, ignore',
                'dummy variable, ignore',
                'dummy variable, ignore',
                'dummy variable, ignore',
                'dummy variable, ignore',
                'dummy variable, ignore',
                'dummy variable, ignore',
                'dummy variable, ignore',
                'dummy variable, ignore',
                ]
        for i in range(len(labels)-10):
            print(f"{labels[i]}={x[i]}")