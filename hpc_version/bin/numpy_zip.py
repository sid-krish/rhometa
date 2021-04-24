import numpy as np

sites1 = np.char.array(['A', 'C', 'G', 'T'])
sites2 = np.char.array(['A', 'C', 'G', 'T'])

m1, m2 = np.meshgrid(sites1, sites2)

m1 = m1.flatten()
m2 = m2.flatten()

m3 = np.stack((m2, m1), axis=1)  # Similar to zip

m3 = m3[m2 != m1]  # remove [A,A], [C,C] etc
