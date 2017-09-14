import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from math import exp, log
from operator import add


#Question 1a
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


