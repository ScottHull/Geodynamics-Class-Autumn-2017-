import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from math import exp, log
from operator import add


# question 1a
df = pd.read_csv('hw3_table4-2.csv', index_col='Isotope')
print("\nRead dataframe:")
print(df)

conc_238U = float(df['Concentration']['238-U']) / float(df['Concentration']['U'])
conc_235U = float(df['Concentration']['235-U']) / float(df['Concentration']['U'])
conc_232Th = float(df['Concentration']['232-Th']) / float(df['Concentration']['232-Th'])
conc_40K = float(df['Concentration']['40-K']) / float(df['Concentration']['K'])
conc_U_Total = float(df['Concentration']['U'])
conc_Th_Total = float(df['Concentration']['232-Th'])
conc_K_Total = float(df['Concentration']['K'])
timespan = list(reversed([i for i in np.linspace(0, (5*10**9), 10000, endpoint=True)]))

def heat_prod_isotope(conc, halflife, H_0, time):
    H_dict = []
    for i in time:
        H = conc * H_0 * exp((i * log(2)) / halflife)
        H_dict.append(H)
    return H_dict

heat_238U = heat_prod_isotope(conc=(conc_238U * conc_U_Total), halflife=float(df['Half-Life']['238-U']),
                        H_0 = float(df['H']['238-U']), time=timespan)
heat_235U = heat_prod_isotope(conc=(conc_235U * conc_U_Total), halflife=float(df['Half-Life']['235-U']),
                        H_0 = float(df['H']['235-U']), time=timespan)
heat_U = list(map(add, heat_235U, heat_238U))
heat_Th = heat_prod_isotope(conc=(conc_232Th * conc_Th_Total), halflife=float(df['Half-Life']['232-Th']),
                        H_0 = float(df['H']['232-Th']), time=timespan)
heat_K = heat_prod_isotope(conc=(conc_40K * conc_K_Total), halflife=float(df['Half-Life']['40-K']),
                        H_0 = float(df['H']['40-K']), time=timespan)
heat_Total = list(map(add, list(map(add, heat_U, heat_Th)), heat_K))


plt.plot(timespan, heat_U, label='U')
plt.plot(timespan, heat_Th, label='Th')
plt.plot(timespan, heat_K, label='K')
plt.plot(timespan, heat_Total, label='Total')
plt.title("Homework 3, Question 1a")
plt.xlabel('Time (Yr)')
plt.ylabel('H (W*Kg-1)')
plt.legend(loc='upper right')
plt.ylim(ymin=0, ymax=(40*10**-12))
plt.xlim(xmin=0, xmax=(5*10**9))
plt.gca().invert_xaxis()
plt.grid()
plt.show()
plt.close()



# # question 1d
timespan = [i for i in np.linspace(0, (1.25*10**9), 10000, endpoint=True)]
def heat_as_time(H, halflife, time):
    H_dict = []
    for i in time:
        H = (7000*10**-6)*(H)*exp((i*log(2))/halflife)
        H_dict.append(H)
    return H_dict

print(heat_as_time(H=2.92*10**-5, halflife=1.25*10**9, time=timespan))

plt.plot(timespan, heat_as_time(H=2.92*10**-5, halflife=1.25*10**9, time=timespan), label='K in core')
plt.title("Homework 3, Question 1d")
plt.xlabel('Time (Yr)')
plt.ylabel('H (W*Kg-1)')
plt.legend(loc='upper right')
plt.grid()
plt.show()
plt.close()




# question 3

if "hw3_prob3.csv" in os.listdir(os.getcwd()):
    os.remove("hw3_prob3.csv")

deltaX = 0.01
Kappa = 0.000005
deltaT = (0.2 * deltaX**2) / Kappa
boundary_T = 1800
slab_T = 300
print("Delta T is: {}".format(deltaT))
max_time_interations = 250
curr_time_iteration = 1
curr_depth = 0
start_temp = 25 + 273.15
print(deltaT)

df = pd.DataFrame({'Depth': [i for i in list(range(11))], "Initial Condition": [slab_T for i in list(range(11))]})
for i in list(range(max_time_interations)):
    time = str(curr_time_iteration)
    prev_time_iteration = str(curr_time_iteration - 1)
    df[time] = boundary_T
    vals = []
    if i != max_time_interations + 1:
        if curr_time_iteration != 1:
            vals.append(boundary_T)
            for row in df['Depth'].ix[1:9]:
                T = (Kappa*deltaT/deltaX**2)*(df[prev_time_iteration][row + 1] + df[prev_time_iteration][row - 1]
                                              -2*df[str(prev_time_iteration)][row]) + df[str(prev_time_iteration)][row]
                vals.append(T)
            vals.append(boundary_T)
        else:
            time = str(curr_time_iteration)
            prev_time_iteration = str(curr_time_iteration - 1)
            df[time] = ''
            df2 = pd.DataFrame({time: []})
            for row in df.index:
                T = (Kappa * deltaT / deltaX ** 2) * (slab_T + slab_T - 2 * slab_T) + slab_T
                vals.append(T)
        df2 = pd.DataFrame({time: vals})
        df[time] = df2[time]
        curr_depth += 1
        curr_time_iteration += 1
        boundary_T += 0
        print("New iteration: {}".format(curr_time_iteration))
    else:
        curr_time_iteration += 1

df.to_csv("hw3_prob3.csv")

df2 = df
df2['Depth'] = df2.index
plt.Figure()
new = df2.iloc[0:11, 3:]
temp = []
for row in new.iterrows():
    index, data = row
    temp.append(data.tolist())
for index, i in enumerate(temp):
    plt.plot(list(range(len(i))), i, label=str(index) + 'km in slab')
plt.grid()
plt.xlabel("Model Iterations (1 iteration = 4 years)")
plt.ylabel('Slab temperature (degK)')
plt.title("HW 3 Problem 3")
plt.legend(loc='lower right')
plt.show()
plt.close()





