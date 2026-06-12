from array import array
import matplotlib.pyplot as plt
import csv
import numpy as np

X = []
Uinit = []
Ufinal = []

# opening the CSV file
with open('initialCondition.csv', mode ='r')as file:
   
  # reading the CSV file
  csvFile = csv.reader(file)
  next(csvFile)  # Skip header
 
  # displaying the contents of the CSV file
  for lines in csvFile:
        X.append(float(lines[0]))
        Uinit.append(float(lines[1]))


with open('final.csv', mode ='r')as file:
   
  # reading the CSV file
  csvFile = csv.reader(file)
  next(csvFile)  # Skip header
 
  # displaying the contents of the CSV file
  for lines in csvFile:
        Ufinal.append(float(lines[1]))
        
#print(X, Uinit)
plt.figure(figsize=(10, 6))
plt.plot(X, Uinit, '.-', color = 'k', label="Initial (t=0)", linewidth=2, markersize=4)
plt.plot(X, Ufinal, '.-', color = 'r', label="Final (t=0.1)", linewidth=2, markersize=4)
plt.xlabel('Position (x)', fontsize=12)
plt.ylabel('Density (ρ)', fontsize=12)
plt.title('Euler 1D: Density Evolution', fontsize=14)
plt.legend(fontsize=11)
plt.grid(True, alpha=0.3)
plt.savefig('euler1d_solution.png', dpi=150, bbox_inches='tight')
print("Graph saved as euler1d_solution.png")
plt.show()