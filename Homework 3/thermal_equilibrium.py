import os
import pandas as pd
import matplotlib.pyplot as plt


# Scott D. Hull, 2017
# Models thermal equilibrium of a subducting slab into the mantle as a function of time.
# Assumption is a vertical slab of 100km thickness (can adjust) going into an isotherm mantle (can also adjust to add
#   a temperature gradient.
# Spits out a pretty plot and a rather large output file.



if "therm_eq.csv" in os.listdir(os.getcwd()):
    os.remove("therm_eq.csv")

Kappa = 0.000005 * 3.154*10**7 * 1*10**-6 # convert seconds to years, m^2 to km^2, the thermal diffusivity
deltaX = 1 # km, change in position per iteration
deltaTime = (0.2 * deltaX**2) / Kappa # years, the time in between iterations
boundary_T = 1800 # degK, isotherm, can add geotherm below in the for loop
slab_T = 300 # degK, downgoing slab initial temperature
slab_thickness = 100 #km, the thickness of the downgoing slab into the mantle
max_time_interations = 1000 # number of model iterations
curr_time_iteration = 1 # tracks current model iteraiton
print("\nModel Parameters:\nKappa: {} km^2/yr\ndeltaTime: {} years\ndeltaX: {} km\nSlab Thickness: {} km\nMax model iterations: {} ({} billion years)".format(
    Kappa, deltaTime, deltaX, slab_thickness, max_time_interations, (max_time_interations*deltaTime) / (1*10**9)))

df = pd.DataFrame({'Depth': [i for i in list(range(slab_thickness + 1))], "Initial Condition": [slab_T for i in list(range(slab_thickness + 1))]})
for i in list(range(max_time_interations)):
    time = str(curr_time_iteration)
    prev_time_iteration = str(curr_time_iteration - 1)
    df[time] = boundary_T
    vals = []
    if i != max_time_interations + 1:
        if curr_time_iteration != 1:
            vals.append(boundary_T)
            for row in df['Depth'].ix[1:slab_thickness - 2]:
                T = (Kappa*deltaTime/deltaX**2)*(df[prev_time_iteration][row + 1] + df[prev_time_iteration][row - 1]
                                              -2*df[str(prev_time_iteration)][row]) + df[str(prev_time_iteration)][row]
                vals.append(T)
            vals.append(boundary_T)
        else:
            time = str(curr_time_iteration)
            prev_time_iteration = str(curr_time_iteration - 1)
            df[time] = ''
            df2 = pd.DataFrame({time: []})
            for row in df.index:
                T = (Kappa * deltaTime / deltaX ** 2) * (slab_T + slab_T - 2 * slab_T) + slab_T
                vals.append(T)
        df2 = pd.DataFrame({time: vals})
        df[time] = df2[time]
        curr_time_iteration += 1
        boundary_T += 0 # add a geotherm
        # print("New iteration: {}".format(curr_time_iteration))
        if curr_time_iteration % 10 == 0:
            print("New iteration: {}".format(curr_time_iteration))
    else:
        curr_time_iteration += 1

df.to_csv("therm_eq.csv")

df2 = df
df2['Depth'] = df2.index
plt.Figure()
new = df2.iloc[0:slab_thickness + 1, 3:]
temp = []
for row in new.iterrows():
    index, data = row
    temp.append(data.tolist())
for index, i in enumerate(temp):
    if index != 0:
        if slab_thickness % index == 0:
            plt.plot(list(range(len(i))), i, label=str(index) + 'km in slab')
plt.grid()
plt.xlabel("Model Iterations (1 iteration = {} years)".format(round(deltaTime, 2)))
plt.ylabel('Slab temperature (degK)')
plt.title("Thermal Equilibrium")
plt.legend(loc='lower right')
plt.show()
plt.close()