import sys
import gzip
import numpy as np


# Read in values for proportion_blw_thrhld for Real Data (Vernot 2016)
with gzip.open(sys.argv[2],'r') as true_values:
        tvs = [float(i) for i in true_values.readline().strip().split()]

tvs_array = np.array(tvs)

#print(tvs_array)

# Read in values for proportion_blw_thrhld (pbt) for simulated data
## The file should look like:
## Model winsize count prop admix1 admix2 admix1_admix2
## Tenn     5     100  0.1   0.01   0.0     0.01_0.0
## Tenn     6     10   0.01  0.01   0.0     0.01_0.0
## Tenn     7     10   0.001 0.01   0.0     0.01_0.0
## Tenn     8     0    0.0   0.01   0.0     0.01_0.0
## ...

with open(sys.argv[1],'r') as simulated_values:
        n1_n2 = -1
        for line in simulated_values:
        #        print(line,file=sys.stdout)
                ls = line.strip().split('\t')
                winsize = ls[1]
                pbt = ls[3]
                n1 = ls[4]
                n2 = ls[5]
                new_n1_n2 = ls[6]

                # We need to store the results for a single set of parameters (i.e. admixture levels)
                # in a single line...so this loop will read the file line by line
                # and add to a single string the values of prop_blw_thrhld
                # for each window size

                ## The output should look like:
                ## admix1 admix2 prop_win1 prop_win2 prop_win3
                ## 0.01   0.0    0.1       0.01      0.001

                if new_n1_n2 != n1_n2:
                        n1_n2 = new_n1_n2
                        newls = n1+' '+n2+' '+pbt
                        #continue
                elif new_n1_n2 == n1_n2:
                        newls = newls+' '+pbt

                if len(newls.strip().split(' ')) == 8:
                        # Once we've collected the pbt values for all the window sizes,
                        # transform into an np.array to perform mse calculation

                        pvs = [float(i) for i in newls.strip().split(' ')]
                        pvs_array = np.array(pvs[2:])

                        mse = np.square(np.subtract(tvs_array,pvs_array)).mean()
                        print(' '.join(map(str,[float(pvs[0]),float(pvs[1]),mse])), file=sys.stdout)
