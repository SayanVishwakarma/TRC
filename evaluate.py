# Import everything from the __init__.py file (likely includes delta_v and stage_optimizer modules)
from __init__ import *

# Define a class called 'evaluate' to perform the rocket evaluation
class evaluate:

    def main(rocket_name="TRC Heavy"):
        out_str=""
        start_time=time.time()
        # Set the rocket name
        #rocket_name = "TRC Heavy"
        output_folder="output files/" + rocket_name + " files"
        os.makedirs(output_folder, exist_ok=True)
        # Define input and output file paths
        input_path = "input files/" + rocket_name + ".json"
        output_path = "output files/" + rocket_name + " files/" + rocket_name + ".txt"
        json_output_path = "output files/" + rocket_name + " files/" + rocket_name + ".json"

        # Open the input JSON file and load rocket parameters into a dictionary
        with open(input_path) as file:
            data = json.load(file)

        old_structural_factors=copy.deepcopy(data["structural factors"])

        out_str=""
        out_str+="**********************************************************************************"+ rocket_name+ "**********************************************************************************"
        out_str+="\n"
        out_str+="---Input Parameters---\n"
        for key,val in data.items():
            out_str+=f"{key} : {val}\n"
        initial_thrust_by_weight = data['thrust by weight'][0]
        Isp = data['average isp']
        orbit = data['target orbit']
        r1 = delta_v.delta_v(initial_thrust_by_weight, Isp, orbit, 2e9, data["Boostback"])
        out_str+=r1.display_breakdown()
        payload = data['payload']
        n_stages = data['number of stages']
        r2 = stage_optimizer.stage_optimizer(
            r1.total_delta_v, r1.boostback_delta_v,
            payload,
            n_stages,
            np.array(data['isps']),
            np.array(data['structural factors']),
            boostback=data["Boostback"]
        )
        rocket = r2.rocket
        out_str+=r2.tabulate_mass_breakdown()
        rocket.add_rocket_data(data)
        rocket.thrust_by_weight = np.array(data["thrust by weight"])
        out_str+=rocket.tabulate_rocket_stages()


        while data["modify factors"]:
            for key,val in data.items():
                if key=="updated structural factors":
                    print(f"Updated structural factors already exist. Using them for calculations")
                    data["structural factors"]=data["updated structural factors"]
            out_str=""
            out_str+="**********************************************************************************"+ rocket_name+ "**********************************************************************************"
            out_str+="\n"
            out_str+="---Input Parameters---\n"
            for key,val in data.items():
                out_str+=f"{key} : {val}\n"
            initial_thrust_by_weight = data['thrust by weight'][0]
            Isp = data['average isp']
            orbit = data['target orbit']
            #r1 = delta_v.delta_v(initial_thrust_by_weight, Isp, orbit, 2e9, data["Boostback"])
            out_str+=r1.display_breakdown()
            payload = data['payload']
            n_stages = data['number of stages']
            r2 = stage_optimizer.stage_optimizer(
                r1.total_delta_v, r1.boostback_delta_v,
                payload,
                n_stages,
                np.array(data['isps']),
                np.array(data['structural factors']),
                boostback=data["Boostback"]
            )
            rocket = r2.rocket
            out_str+=r2.tabulate_mass_breakdown()
            rocket.add_rocket_data(data)
            rocket.thrust_by_weight = np.array(data["thrust by weight"])
            out_str+=rocket.tabulate_rocket_stages()
            stage_1_mass=tank.compute_stage_mass(data["tank diameters"]/2,rocket.step_fuel_tank_height[0],rocket.step_oxidiser_tank_height[0],data["step diameters"][0],rocket.step_height[0],rocket.step_engine_number[0])
            stage_2_mass=tank.compute_stage_mass(data["tank diameters"]/2,rocket.step_fuel_tank_height[1],rocket.step_oxidiser_tank_height[1],data["step diameters"][1],rocket.step_height[1],rocket.step_engine_number[1])
            s_factors=data["structural factors"]
            percent_error_1=np.abs(stage_1_mass-(r2.step_masses[0]*r2.epsilons[0]))/(r2.step_masses[0]*r2.epsilons[0])*100
            percent_error_2=np.abs(stage_2_mass-(r2.step_masses[1]*r2.epsilons[1]))/(r2.step_masses[1]*r2.epsilons[1])*100
            print(f"Stage 1 mass off by {percent_error_1}%")
            print(f"Stage 2 mass off by {percent_error_2}%")    
            if percent_error_1<10 and percent_error_2<10:
                print("\nPercent error in mass acceptable. Proceeding")
                print(f"Updated structural factors are {np.round(s_factors,4)}")
                out_str+=f"Updated structural factors are {np.round(s_factors,4)}\n"
                break
            if stage_1_mass>(r2.step_masses[0]*r2.epsilons[0]):
                s_factors[0]+=0.0001
                print(f"Overestimated Stage 1 Structural Factor, updating it to {s_factors[0]}")
            else:
                s_factors[0]-=0.0001
                print(f"Underestimated Stage 1 Structural Factor, updating it to {s_factors[0]}")
            
            if stage_2_mass>(r2.step_masses[1]*r2.epsilons[1]):
                s_factors[1]+=0.0001
                print(f"Overestimated Stage 2 Structural Factor, updating it to {s_factors[1]}")
            else:
                s_factors[1]-=0.0001
                print(f"Underestimated Stage 2 Structural Factor, updating it to {s_factors[1]}")

        if data['modify factors']:            
            data["updated structural factors"]=[np.round(s_factors,4)[i] for i in range(0,len(np.round(s_factors,4)))]
            data["structural factors"]=copy.deepcopy(old_structural_factors)
        
        with open(input_path, 'w') as file:
            json.dump(data,file,indent=4)

        out_str+="\n********************************************************************************************************************************\n"


        #START TIME HAS BEEN FIXED REMEMBR TO CHANGE IT

        traj=rocket.create_trajectory_object()
        #traj.mins=(5, 60.003144304746236, 0.1237625780263158, 420.4619684100563)
        #traj.mins=(5, 80.29662781038536, 0.11436489662203948, 480.0)
        #traj.mins=(5, 184.48275862068965, 0.03482758620689656, 480.0)
        traj.simulate_trajectory(n=80,n_epochs=1,dynamic_n=True)
                                 #,start_time=20,delta_start_time=15
                                 #,end_time=100,delta_end_time=50
                                 #,beta=0.1,delta_beta=0.09,
                                 #thrust_cutoff=5*60,delta_thrust_cutoff=3*60)

        #traj.gp_optimiser(n_calls=100)
                                
        if not(traj.mins==(0,0,0,0) or traj.mins==[0,0,0,0]):
            if "optimal launch parameters" in data.keys():
                data["optimal launch parameters"]+=[traj.mins]
            else:
                data["optimal launch parameters"]=[traj.mins]
        with open(input_path, 'w') as file:
            json.dump(data,file,indent=4)
        #print("Simulating final trajectory\n")
        out_str+="									---------------Launch Trajectories---------------\n"
        out_str+="********************************************************************************************************************************\n"
        for vals in data["optimal launch parameters"]:
            traj.__init__(traj.data_copy)
            traj.mins=vals
            out_str+="\n\nOptimal Launch Parameters = " + str(traj.mins) + "\n"
            out_str+=traj.model(vector_start_time=traj.mins[0],vector_end_time=traj.mins[1],beta_max=traj.mins[2],thrust_cutoff_time=traj.mins[3],
                                post_burnout=False,
                                total_time=1.7*60*60,
                                display_breakdown=True)
            #traj.plot_altitudes()
            #traj.plotter()
            #print(traj.total_rotation*180/pi)
            traj.save_trajectory()
            #traj.plot_magnitudes()#aoa=True)
            #traj.dyn_plotter(animation_speed=5000)
            print()
                
        
        # Redirect stdout to the output text file
        sys.stdout = open(output_path, 'w')
        print(out_str)
        sys.stdout.close()
        sys.stdout = sys.__stdout__

        print(f"\nExecution time taken : {time.time()-start_time} s")
        print(f"Results saved to {output_path}")
        #traj.model()
        #traj.plot_altitudes()
        #traj.plotter()
        #traj.plotter_inertial_frame()
        
        

    def sweep(rocket_name):
        mass_difference_threshold=3
        out=[]
        start_time=time.time()
        # Set the rocket name
        #rocket_name = "TRC Heavy"
        output_folder="output files/" + rocket_name + " files"
        os.makedirs(output_folder, exist_ok=True)
        # Define input and output file paths
        input_path = "input files/" + rocket_name + ".json"
        output_excel = "output files/" + rocket_name + " files/" + "sweep data.xlsx"

        # Open the input JSON file and load rocket parameters into a dictionary
        with open(input_path) as file:
            data = json.load(file)

        payload_range=[50e3]#np.arange(50e3,100e3,10e3)
        tank_diameter_range=[7,7.25,7.5,7.75,8,8.25,8.5]#np.arange(6,9,0.25)
        boostback_1_range=[1500]#[1500,2000]
        boostback_2_range=[2000]#[2000,2500]
        eng_thrust_range=[[[2171.39e3,2600.39e3],[1.23]],[[2395e3,2800e3],[1.32]],[[280e4,3285e3],[1.41]]]
        #eng_exit_dia_range=[1.23,1.32]
        n_iter=len(payload_range)*len(tank_diameter_range)*len(boostback_1_range)*len(boostback_2_range)*len(eng_thrust_range)
        sr_no=1
        initial_thrust_by_weight = data['thrust by weight'][0]
        Isp = data['average isp']
        orbit = data['target orbit']
        r1 = delta_v.delta_v(initial_thrust_by_weight, Isp, orbit, 2e9, data["Boostback"])
        for payload in payload_range:
            for tank_diameter in tank_diameter_range:
                for boostback_1 in boostback_1_range:
                    for boostback_2 in boostback_2_range:
                        for eng_data in eng_thrust_range:
                            eng_thrust=eng_data[0]
                            eng_dia=eng_data[1][0]
                            out_dict={}
                            out_dict["Sr. No."]=sr_no
                            out_dict["Payload"]=payload
                            out_dict["Step Diameters"]=tank_diameter+0.5
                            data["payload"]=payload
                            data["tank diameters"]=tank_diameter
                            data["step diameters"]=[tank_diameter+0.5,tank_diameter+0.5]
                            data["engine thrust"]=eng_thrust
                            out_dict["Thrust"]=eng_thrust[0]
                            out_dict["TWR Initial"]=data["thrust by weight"][0]
                            out_dict["Engine exit dia"]=eng_dia
                            output_path = "output files/" + rocket_name + " files/" + str(sr_no) + ".txt"
                            out_str=""
                            out_str+="**********************************************************************************"+ rocket_name+ "**********************************************************************************"
                            out_str+="\n"
                            out_str+="---Input Parameters---\n"
                            for key,val in data.items():
                                out_str+=f"{key} : {val}\n"
                            #initial_thrust_by_weight = data['thrust by weight'][0]
                            #Isp = data['average isp']
                            #orbit = data['target orbit']
                            #r1 = delta_v.delta_v(initial_thrust_by_weight, Isp, orbit, 2e9, data["Boostback"])
                            r1.boostback_delta_v=[boostback_1,boostback_2]
                            out_str+=r1.display_breakdown()
                            payload = data['payload']
                            n_stages = data['number of stages']
                            r2 = stage_optimizer.stage_optimizer(
                                r1.total_delta_v, r1.boostback_delta_v,
                                payload,
                                n_stages,
                                np.array(data['isps']),
                                np.array(data['structural factors']),
                                boostback=data["Boostback"]
                            )
                            rocket = r2.rocket
                            out_str+=r2.tabulate_mass_breakdown()
                            rocket.add_rocket_data(data)
                            rocket.thrust_by_weight = np.array(data["thrust by weight"])
                            out_str+=rocket.tabulate_rocket_stages()

                            while data["modify factors"]:
                                out_str=""
                                out_str+="**********************************************************************************"+ rocket_name+ "**********************************************************************************"
                                out_str+="\n"
                                out_str+="---Input Parameters---\n"
                                for key,val in data.items():
                                    out_str+=f"{key} : {val}\n"
                                initial_thrust_by_weight = data['thrust by weight'][0]
                                Isp = data['average isp']
                                orbit = data['target orbit']
                                #r1 = delta_v.delta_v(initial_thrust_by_weight, Isp, orbit, 2e9, data["Boostback"])
                                #r1.boostback_delta_v=[boostback_1,boostback_2]
                                out_str+=r1.display_breakdown()
                                payload = data['payload']
                                n_stages = data['number of stages']
                                r2 = stage_optimizer.stage_optimizer(
                                    r1.total_delta_v, r1.boostback_delta_v,
                                    payload,
                                    n_stages,
                                    np.array(data['isps']),
                                    np.array(data['structural factors']),
                                    boostback=data["Boostback"]
                                )
                                rocket = r2.rocket
                                out_str+=r2.tabulate_mass_breakdown()
                                rocket.add_rocket_data(data)
                                rocket.thrust_by_weight = np.array(data["thrust by weight"])
                                out_str+=rocket.tabulate_rocket_stages()
                                stage_1_mass=tank.compute_stage_mass(data["tank diameters"]/2,rocket.step_fuel_tank_height[0],rocket.step_oxidiser_tank_height[0],data["step diameters"][0],rocket.step_height[0],rocket.step_engine_number[0])
                                stage_2_mass=tank.compute_stage_mass(data["tank diameters"]/2,rocket.step_fuel_tank_height[1],rocket.step_oxidiser_tank_height[1],data["step diameters"][1],rocket.step_height[1],rocket.step_engine_number[1])
                                s_factors=data["structural factors"]
                                percent_error_1=np.abs(stage_1_mass-(r2.step_masses[0]*r2.epsilons[0]))/(r2.step_masses[0]*r2.epsilons[0])*100
                                percent_error_2=np.abs(stage_2_mass-(r2.step_masses[1]*r2.epsilons[1]))/(r2.step_masses[1]*r2.epsilons[1])*100
                                print(f"Stage 1 mass off by {percent_error_1}%")
                                print(f"Stage 2 mass off by {percent_error_2}%")    
                                if percent_error_1<mass_difference_threshold and percent_error_2<mass_difference_threshold:
                                    print("Percent error in mass acceptable. Proceeding")
                                    print(f"Updated structural factors are {np.round(s_factors,4)}")
                                    out_str+=f"Updated structural factors are {np.round(s_factors,4)}\n"
                                    break
                                if stage_1_mass>(r2.step_masses[0]*r2.epsilons[0]):
                                    s_factors[0]+=0.0001
                                    print(f"Overestimated Stage 1 Structural Factor, updating it to {s_factors[0]}")
                                else:
                                    s_factors[0]-=0.0001
                                    print(f"Underestimated Stage 1 Structural Factor, updating it to {s_factors[0]}")
                                
                                if stage_2_mass>(r2.step_masses[1]*r2.epsilons[1]):
                                    s_factors[1]+=0.0001
                                    print(f"Overestimated Stage 2 Structural Factor, updating it to {s_factors[1]}")
                                else:
                                    s_factors[1]-=0.0001
                                    print(f"Underestimated Stage 2 Structural Factor, updating it to {s_factors[1]}")
                                    
                                data["structural factors"]=np.round(s_factors,4)
                            out_str+="********************************************************************************************************************************\n"
                            out_dict["Total Rocket Mass"]=np.sum(r2.step_masses)
                            out_dict["Rocket Height"]=rocket.rocket_height
                            out_dict["Stage 1 Mass"]=r2.step_masses[0]
                            out_dict["Stage 2 Mass"]=r2.step_masses[1]
                            out_dict["Stage 1 Structural Factor"]=data["structural factors"][0]
                            out_dict["Stage 2 Structural Factor"]=data["structural factors"][1]
                            out_dict["Boostback delta v stage 1"]=boostback_1
                            out_dict["Boostback delta v stage 2"]=boostback_2
                            out_dict["Engine count (step 1)"]=rocket.step_engine_number[0]
                            out_dict["Engine count (step 2)"]=rocket.step_engine_number[1]

                            """traj=rocket.create_trajectory_object()
                            
                            best_time,best_beta=trajectory.optimize_betas_2(traj.data,n_betas=50,n_ts=50,t1=6,t2=300,beta_max=1.5)
                            out_str+=traj.model(vector_end_time=best_time,beta_max=best_beta,total_time=3*60*60,post_burnout=False,)
                            
                            out_str+=traj.model(beta_max=0,post_burnout=False)
                            """
                            
                            # Redirect stdout to the output text file
                            sys.stdout = open(output_path, 'w')
                            print(out_str)
                            sys.stdout.close()
                            sys.stdout = sys.__stdout__
                            out.append(out_dict)
                            print(f"{round(sr_no/n_iter*100,2)}% completed")
                            sr_no+=1
                            time.sleep(0.5)
        df = pd.DataFrame(out)
        df.to_excel(output_excel, index=False)
        print(f"Sweep complete. Results saved to {output_excel}")
    
    #main("TRC Heavy 2")
    main("TRC Superheavy v3")
    #main("TRC Heavy")
    #sweep("TRC Heavy tank dia sweep")
