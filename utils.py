import math
def clean_fuel_data2(data,i):
    if i==len(data)-1:
        return data[i-1]
    elif math.isnan(data[i+1]) and math.isnan(data[i-1]):
        if math.isnan(data[i+1]):
            return clean_fuel_data2(data,i+1)
        else:
            return data[i]
    else:
        return (data[i+1]+data[i-1])/2
    return
def clean_fuel_data(path="data files/LOXMETHANE.txt"):
    fuel_data=np.genfromtxt(path)[1:].T
    fuel_data_heads=["O/F","temp","isp","mw","cstar"]
    for j in range(len(fuel_data)):
        data=fuel_data[j]
        for i in range(len(data)):
            if math.isnan(data[i]):
                data[i]=clean_fuel_data2(data,i)
        fuel_data[j]=data
    sys.stdout = open("data files/LOXMETHANE cleaned.txt", 'w')
    print(f"O/F\ttemp\tisp\tmw\tcstar\tgamma")
    for i in range(len(data)):
        print(f"{fuel_data[0][i]}\t{fuel_data[1][i]}\t{fuel_data[2][i]}\t{fuel_data[3][i]}\t{fuel_data[4][i]}\t{fuel_data[5][i]}")
    sys.stdout.close()
    sys.stdout = sys.__stdout__