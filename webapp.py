from __init__ import *
from taipy.gui import Gui, Icon, navigate
import taipy.gui.builder as tgb

cats=["lol","lol2"]
rocket_name="Name of the rocket"
target_orbit=500
average_isp=300
number_of_stages=2
thrust_by_weight=[1.55,0.56]
isps=[311,348]
structural_factors=[0.05136,0.0434]
payload=22800
Boostback=[True,False]
step_diameters=[3.7,3.7]
fuel_density=820
oxidiser_density=1141
tank_diameters=3.66
o_f=2.6
engine_height=2.3
engine_thrust=[845e3,915e3]
data_dict={}

data={"names": ["falcon 9", "falcon", "saturn v", " starship", "TRC Heavy 2", "TRC Heavy"]}
rocket_list = data["names"]
selected_rocket = "TRC"

def change_selected_rocket(state, action, info):
    print("lol")

def menu_option_selected(state, action, info):
    page = info["args"][0]
    navigate(state, to=page)

def go_to_add_rocket(state):
    navigate(state, to='page3')

def go_to_display_rocket_data(state):
    global selected_rocket
    selected_rocket = state["selected_rocket"]
    print(f"Selected Rocket: {selected_rocket}")
    navigate(state, to='page2')

def create_rocket():
    data_dict["rocket name"] = rocket_name
    data_dict["target orbit"] = target_orbit
    data_dict["number of stages"] = number_of_stages
    data_dict["payload"] = payload
    data_dict["thrust by weight"] = thrust_by_weight
    data_dict["average isp"] = average_isp
    data_dict["isps"] = isps
    data_dict["structural factors"] = structural_factors
    data_dict["Boostback"] = Boostback
    data_dict["step diameters"] = step_diameters
    data_dict["fuel density"] = fuel_density
    data_dict["oxidiser density"] = oxidiser_density
    data_dict["tank diameters"] = tank_diameters
    data_dict["o/f"] = o_f
    data_dict["engine height"] = engine_height
    data_dict["engine thrust"] = engine_thrust
    print("Rocket data received:", data_dict)
    nam=data_dict["rocket name"]
    #output_folder="input files/" + nam + " files"
    #os.makedirs(output_folder, exist_ok=True)
    with open(f"input files/{nam}.json", "w") as file:
        json.dump(data_dict, file, indent=4)
    print(f"Rocket {nam} created successfully!")

with tgb.Page() as root_page:
    tgb.menu(
        label="Menu",
        lov=[
            ("page1", Icon("output files/TRC Heavy 2 files/3d trajectory.png", "Home")),
            ("page2", Icon("output files/TRC Heavy 2 files/3d trajectory.png", "Rockets")),
        ],
        on_action=menu_option_selected
    )

with tgb.Page() as homepage:
    with tgb.part(class_name="container"):
        tgb.text("# Welcome to TRC's Rocket Browser", mode="md")
    tgb.html("br")
    with tgb.part(class_name="card"):
        #with tgb.layout(columns="3 1"):    
            tgb.button("Add Rocket", on_action=go_to_add_rocket)
            tgb.html("br")
            tgb.selector(label="Rocket",lov=rocket_list,value="{selected_rocket}",dropdown=True, on_action=go_to_display_rocket_data)
            #if selected_rocket!="":
            tgb.text(f"Selected Rocket: {selected_rocket}", mode="md")

with tgb.Page() as rocket_data_page:
    with tgb.part(class_name="card"):
        with open(f"output files/{selected_rocket} files/{selected_rocket}.txt") as file:
            contents= file.read()
        tgb.text(contents, mode="txt")
    
with tgb.Page() as add_rocket:

    tgb.text("Creating a Rocket, enter the inputs.")
    tgb.html("br")
    tgb.input(bind="{rocket_name}", label="Enter Rocket's Name",value=rocket_name)
    tgb.input(bind="{target_orbit}", label="Target Orbit (in km)",value=target_orbit)
    tgb.input(bind="{number_of_stages}", label="Number of stages",value=number_of_stages)
    tgb.input(bind="{payload}", label="Payload",value=payload)
    tgb.input(bind="{thrust_by_weight}", label="Thrust by weight, array of required size",value=thrust_by_weight)
    tgb.input(bind="{average_isp}", label="Average ISP",value=average_isp)
    tgb.input(bind="{isps}", label="ISPs, array of required size",value=isps)
    tgb.input(bind="{structural_factors}", label="Structural Factors, array of required size",value=structural_factors)
    tgb.input(bind="{Boostback}", label="Boostback, array of required size",value=Boostback)
    tgb.input(bind="{step_diameters}", label="Step Diameters, array of required size",value=step_diameters)
    tgb.input(bind="{fuel_density}", label="Fuel Density",value=fuel_density)
    tgb.input(bind="{oxidiser_density}", label="Oxidiser Density",value=oxidiser_density)
    tgb.input(bind="{tank_diameters}", label="Tank Diameters",value=tank_diameters)
    tgb.input(bind="{o_f}", label="Oxidiser to Fuel Ratio",value=o_f)
    tgb.input(bind="{engine_height}", label="Engine Height",value=engine_height)
    tgb.input(bind="{engine_thrust}", label="Engine Thrust, array of required size",value=engine_thrust)
    
    tgb.button("Create Rocket", on_action=create_rocket)

pages = {"/": root_page, "page1": homepage, "page2": rocket_data_page, "page3": add_rocket}

Gui(pages=pages).run(title="TRC", port=5000, use_reloader=True)



