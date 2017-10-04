import os
import pandas as pd
import matplotlib.pyplot as plt
from math import cos, pi, exp, sqrt


# Scott D. Hull, 2017





if "therm_eq.csv" in os.listdir(os.getcwd()):
    os.remove("therm_eq.csv")

Kappa = 0.000005 # convert seconds to years, m^2 to km^2, the thermal diffusivity
deltaX = 1 # m, change in position per iteration
deltaTime = (0.2 * deltaX**2) / Kappa # years, the time in between iterations
period = (2*pi)/(pi*10**7)
surface_maxT = 800
def boundary_T(N_timesteps, deltaT=deltaTime, period=period, T_not_surf=surface_maxT):
    boundary_t = T_not_surf + float(100*cos(period*(deltaT*N_timesteps)))
    return boundary_t
model_depth = 100 # m
T_o = 800
max_time_iterations = 1000
curr_time_iteration = 1

def given_formula(depth, time, T_o=T_o, period=period, K=Kappa, deltaT=100):
    T = T_o + deltaT*exp(-depth*sqrt(period/(2*K)))*cos((period*time)+depth*sqrt(period/(2*K)))
    return T


df = pd.DataFrame({'Depth': [i for i in list(range(model_depth + 1))], "Initial Condition": [T_o for i in list(range(model_depth + 1))]})
for i in list(range(max_time_iterations)):
    time = str(curr_time_iteration)
    prev_time_iteration = str(curr_time_iteration - 1)
    df[time] = boundary_T(N_timesteps=float(i))
    vals = []
    if i != max_time_iterations + 1:
        if curr_time_iteration != 1:
            vals.append(boundary_T(N_timesteps=i))
            for row in df['Depth'].ix[1:model_depth - 2]:
                T = (Kappa*deltaTime/deltaX**2)*(df[prev_time_iteration][row + 1] + df[prev_time_iteration][row - 1] - 2*df[str(prev_time_iteration)][row]) + df[str(prev_time_iteration)][row]
                vals.append(T)
            vals.append(T_o)
        else:
            time = str(curr_time_iteration)
            prev_time_iteration = str(curr_time_iteration - 1)
            df[time] = ''
            df2 = pd.DataFrame({time: []})
            for row in df.index:
                T = (Kappa * deltaTime / deltaX ** 2) * (T_o + T_o - 2 * T_o) + T_o
                vals.append(T)
        df2 = pd.DataFrame({time: vals})
        df[time] = df2[time]
        curr_time_iteration += 1
        print("New iteration: {}".format(curr_time_iteration))
    else:
        curr_time_iteration += 1

df.to_csv("therm_eq.csv")

df2 = df
df2['Depth'] = df2.index
plt.Figure()
new = df2.iloc[0:model_depth + 1, 3:]
temp = []
for row in new.iterrows():
    index, data = row
    temp.append(data.tolist())
for i in df2:
    if i != 'Depth' and i != "" and i != "Initial Condition":
        if float(i) % 500 == 0:
            formula_temps = []
            time = float(i) * deltaTime
            plt.plot(df2[i], df2['Depth'], label='Time: {:.2e} yrs (Numerical)'.format(round(time, 1)))
            for z in df2['Depth']:
                formula_temps.append(given_formula(depth=z, time=float(i)*deltaTime))
            plt.plot(formula_temps, df2['Depth'], label='Time: {:.2e} yrs (Analytical)'.format(round(time, 1)))


# for index, i in enumerate(temp):
#     if index != 0:
#         if model_depth % index == 0:
#             plt.plot(list(range(len(i))), i, label=str(index) + 'depth')
plt.vlines(T_o, ymin=0, ymax=model_depth, linestyles=':', label='Initial Mantle T', color='r', linewidth=2)
plt.gca().invert_yaxis()
plt.grid()
plt.xlabel("Temperature (degK)")
plt.ylabel('Depth (km)')
plt.title("Thermal Equilibrium In Oscillating Surface Temperature Scenario")
plt.legend(loc='lower right')
plt.show()
plt.close()

time2 = list(range(round(max_time_iterations+1)))[1:]
plt.figure()
for row in df2.index:
    if df2['Depth'][row] % 5 == 0 and df2['Depth'][row] <= 15:
        temps = []
        formula_temps = []
        for z in df2:
            temps.append(df2[z][row])
        for r in time2:
            formula_temps.append(given_formula(depth=df2['Depth'][row], time=r*deltaTime))
        plt.plot(time2, temps[2:], label='Depth: {} (Numerical)'.format(df2['Depth'][row]))
        plt.plot(time2, formula_temps, label='Depth: {} (Analytical)'.format(df2['Depth'][row]))
plt.grid()
plt.xlabel("Time Iterations (1 iteration = {} sec)".format(deltaTime))
plt.ylabel('Temperature (degK)')
plt.title("Thermal Equilibrium In Oscillating Surface Temperature Scenario")
plt.legend(loc='lower right')
plt.show()
plt.close()

