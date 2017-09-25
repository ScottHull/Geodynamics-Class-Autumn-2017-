import os
import pandas as pd
import matplotlib.pyplot as plt
import time


# Scott D. Hull, 2017
# Models thermal equilibrium of a subducting body into the mantle as a function of time.
# Assumption is a vertical body of 100km thickness (can adjust) going into an isotherm mantle (can also adjust to add
#   a temperature gradient.
# Spits out a pretty plot and a rather large output file.



if "therm_eq.csv" in os.listdir(os.getcwd()):
    os.remove("therm_eq.csv")

Kappa = 0.000005 * 3.154*10**7 * 1*10**-6 # convert seconds to years, m^2 to km^2, the thermal diffusivity
deltaX = 1 # km, change in position per iteration
deltaTime = (0.2 * deltaX**2) / Kappa # years, the time in between iterations
boundary_T = 0 # degK, isotherm, can add geotherm below in the for loop
body_T = 1800 # degK, downgoing body initial temperature
body_thickness = 100 #km, the thickness of the downgoing body into the mantle
max_time_interations = 26000 # number of model iterations
curr_time_iteration = 1 # tracks current model iteraiton
print("\nModel Parameters:\nKappa: {} km^2/yr\ndeltaTime: {} years\ndeltaX: {} km\nSlab Thickness: {} km\nMax model iterations: {} ({} billion years)\n".format(
    Kappa, deltaTime, deltaX, body_thickness, max_time_interations, (max_time_interations*deltaTime) / (1*10**9)))


depth = [0]
half_depth = list(range(round((body_thickness + 2)/2)))
for i in half_depth:
    depth.append(i)
adjusted_depth = []
for i in list(reversed(depth))[:-2]:
    adjusted_depth.append(i)
for i in half_depth:
    adjusted_depth.append(i)


df = pd.DataFrame({'Depth': adjusted_depth, "Initial Condition": [body_T for i in list(range(len(adjusted_depth)))]})
for i in list(range(max_time_interations)):
    time = str(curr_time_iteration)
    prev_time_iteration = str(curr_time_iteration - 1)
    df[time] = boundary_T
    vals = []
    if i != max_time_interations + 1:
        if curr_time_iteration != 1:
            vals.append(boundary_T)
            # for row in df['Depth'].ix[1:len(adjusted_depth) - 1]:
            for row in df.index[1:len(adjusted_depth) - 1]:
                curr_radius = int(df['Depth'][row])
                # print("Index: {}\nCurr_Radius: {}\nt-deltat: {}\nprev_row: {}\nrow+1: {}\nrow-1: {}".format(
                #     row, curr_radius, prev_time_iteration, df[str(prev_time_iteration)][row],
                #     df[str(prev_time_iteration)][row + 1], df[str(prev_time_iteration)][row - 1]))
                if int(df['Depth'][row]) == 0:
                    T = df[prev_time_iteration][row] + ((Kappa * deltaTime) * (((df[str(prev_time_iteration)][row + 1] +
                                                                             df[str(prev_time_iteration)][row - 1] - (
                                                                             2 * df[str(prev_time_iteration)][
                                                                                 row])) / deltaX ** 2)))
                    vals.append(T)

                else:
                    T = df[prev_time_iteration][row] + ((Kappa*deltaTime)*(((df[str(prev_time_iteration)][row + 1] +
                    df[str(prev_time_iteration)][row - 1] - (2*df[str(prev_time_iteration)][row]))/deltaX**2)) +
                                                        ((2/curr_radius) * ((df[str(prev_time_iteration)][row + 1] -
                    df[str(prev_time_iteration)][row - 1])/(2*deltaTime))))

                    # term1 = df[prev_time_iteration][row]
                    # term2 = (Kappa*deltaTime)*(((df[str(prev_time_iteration)][row + 1] +
                    # df[str(prev_time_iteration)][row - 1] - (2*df[str(prev_time_iteration)][row]))/deltaX**2))
                    # term3 = (Kappa*deltaTime)*((2/curr_radius) * ((df[str(prev_time_iteration)][row + 1] -
                    # df[str(prev_time_iteration)][row - 1])/(2*deltaTime)))
                    # print('\nTerm 1: {}     Term 2: {}      Term 3: {}\n'.format(term1, term2, term3))
                    vals.append(T)
            vals.append(boundary_T)
        else:
            time = str(curr_time_iteration)
            prev_time_iteration = str(curr_time_iteration - 1)
            df[time] = ''
            df2 = pd.DataFrame({time: []})
            vals.append(boundary_T)
            for row in df.index[:-2]:
                T = (Kappa * deltaTime / deltaX ** 2) * (body_T + body_T - 2 * body_T) + body_T
                vals.append(T)
            vals.append(boundary_T)
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
fig1 = plt.figure()
fig2 = plt.figure()
ax1 = fig1.add_subplot(111)
ax2 = fig2.add_subplot(111)
new = df2.iloc[0:body_thickness + 1, 3:]
temp = []
for row in new.iterrows():
    index, data = row
    temp.append(data.tolist())
for index, i in enumerate(temp):
    if index != 0:
        if body_thickness % index == 0:
            ax1.plot(list(range(len(i))), i, label=str(index) + 'km in body')
for i in list(range(max_time_interations)):
    if i % 5000 == 0 and i != 0:
        temps = df[str(i)]
        depths = df['Depth']
        ax2.plot(temps, depth, label='{} years'.format(round(i*deltaTime, 2)))
plt.grid()
ax1.set_xlabel("Model Iterations (1 iteration = {} years)".format(round(deltaTime, 2)))
ax2.set_xlabel("Temperature (degK)")
ax2.set_ylabel("Radius (km)")
ax1.set_ylabel('Body temperature (degK)')
ax1.set_title("Thermal Equilibrium (Time vs Temperature)")
ax2.set_title("Thermal Equilibrium (Radius vs Temperature)")
ax1.legend(loc='upper right')
ax2.legend(loc='upper right')
plt.show()
plt.close()
