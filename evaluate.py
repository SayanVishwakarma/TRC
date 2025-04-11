# Import everything from the __init__.py file (likely includes delta_v and stage_optimizer modules)
from __init__ import *

# Define a class called 'evaluate' to perform the rocket evaluation
class evaluate:
    # Set the rocket name
    rocket_name = "falcon 9"

    # Define input and output file paths
    input_path = "input files/" + rocket_name + ".json"
    output_path = "output files/" + rocket_name + ".txt"

    # Open the input JSON file and load rocket parameters into a dictionary
    with open(input_path) as file:
        data = json.load(file)

    # Redirect stdout to the output text file
    sys.stdout = open(output_path, 'w')

    # Print a header with the rocket name
    print("**********************************************************************************", rocket_name, "**********************************************************************************")

    # Extract initial thrust-to-weight ratio and average specific impulse
    initial_thrust_by_weight = data['thrust by weight'][0]
    Isp = data['average isp']

    # Extract the target orbit altitude
    orbit = data['target orbit']

    # Calculate the required delta-v using the delta_v module
    # Parameters: T/W, Isp, target orbit (km), payload mass (fixed here as 2e9), and boostback flag
    r1 = delta_v.delta_v(initial_thrust_by_weight, Isp, orbit, 2e9, data["Boostback"])

    # Display the breakdown of the delta-v calculation
    r1.display_breakdown()

    # Extract payload mass and number of stages from the input data
    payload = data['payload']
    n_stages = data['number of stages']

    # Check if boostback is enabled in the input data
    if data["Boostback"] == "True":
        # Use stage_optimizer with boostback capability
        r2 = stage_optimizer.stage_optimizer(
            [r1.total_delta_v, r1.boostback_delta_v],
            payload,
            n_stages,
            np.array(data['isps']),
            np.array(data['structural factors']),
            boostback=True
        )
    else:
        # Use stage_optimizer without boostback
        r2 = stage_optimizer.stage_optimizer(
            r1.total_delta_v,
            payload,
            n_stages,
            np.array(data['isps']),
            np.array(data['structural factors'])
        )

    # Get the resulting rocket object from the optimizer
    rocket = r2.rocket

    # Add additional rocket configuration data (e.g., tank, engine, fuel) to the rocket object
    rocket.add_rocket_data(data)

    # Assign thrust-to-weight ratios for each stage
    rocket.thrust_by_weight = np.array(data["thrust by weight"])

    # Optionally, display simulation and sizing (currently commented out)
    # rocket.display_states()
    # rocket.size_rocket()
    # rocket.simulate_burn_2()

    # Tabulate and print the stage breakdown of the rocket
    rocket.tabulate_rocket_stages()

    # Print a footer to mark end of output
    print("********************************************************************************************************************************")

'''
# Old interactive code block for manual input (now replaced by file-based approach)

initial_thrust_by_weight = float(input("ENTER INITIAL THRUST BY WEIGHT:"))
Isp = float(input("ENTER AVERAGE ISP:"))
orbit = float(input("ENTER TARGET ORBIT:"))
r1 = delta_v.delta_v(initial_thrust_by_weight, Isp, orbit, 2e8)
print("REQUIRED DELTA_V=", r1.total_delta_v)
payload = float(input("ENTER PAYLOAD MASS:"))
n_stages = int(input("ENTER NUMBER OF STAGES:"))
r2 = stage_optimizer.stage_optimizer(r1.total_delta_v, payload, n_stages)
rocket = r2.rocket
print(rocket)
rocket.display_states()
'''
