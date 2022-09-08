#### Matthew Lollar 
#### Latest Update: December 9th, 2021
#### Tested in: Python 3.x
#### Power analysis simulator for Uniparental by X/Autosome two-locus recessive incompatibilities

import argparse
import numpy as np
import pandas as pd
import random as rd
import statsmodels.stats.contingency_tables as sm
np.seterr(all='raise') # Set all numpy errors to raise for exception handling in tests

### Input arguments
parser = argparse.ArgumentParser(description='Power Analysis Simulator, last update: 10/2021')
parser.add_argument('--i', help='Input File', required=True)
parser.add_argument('--o', help='Output File', required=True)
parser.add_argument('--g', help='Focal window 1 Genotype', type=int, required=True)
parser.add_argument('--gg', help='Focal window 2 Genotype', type=int, required=True)
parser.add_argument('--s', help='Number of Sterile Focal', type=int, required=True)
parser.add_argument('--f', help='Number of Fertile', type=int, required=True)
parser.add_argument('--bs', help='Background Sterility (As percent, currenty only accetps integer value)', type=int, required=True)
args = parser.parse_args()
number_sterile = args.s
number_fertile = args.f/2

### Chromosome window boundaries (zero-indexed)
X_end = 545
Chr3_end = 2579

### Load Array, Transpose (Inidividuals read in row-wise individual and converted to column-wise)
num_array = np.loadtxt(args.i, dtype=int, delimiter='\t')
num_array = num_array.transpose()
df = pd.DataFrame(num_array)

### Randomly assign uniparental loci to individuals
uni_choice = [0,2]   #Focal locus of Y/Mitochondria takes on value of 0 or 2
uni_windows = []   #List focal loci for individuals
for x in range(0,df.shape[1]):
	uni_windows.append(rd.choice(uni_choice))

### Randomly selected incompatibility locus
x_inc = rd.randint(1,(X_end-1)) 
a_inc = rd.randint(X_end, (Chr3_end-1))

### Initialize list of IDs for groups X and A
sterile_ids_X = []     # Sterile IDs X-Uni
fertile_f_ids_X = []   # Fertile Focal IDs X-Uni, equal number as fertile NF
fertile_nf_ids_X = []  # Fertile Non-focal IDs X-Uni, equal number as fertile F
sterile_ids_A = []     # Sterile A-Uni
fertile_f_ids_A = []   # Fertile Focal A-Uni
fertile_nf_ids_A = []  # Fertile Non-focal A-Uni

Ymito_focal = args.g # Focal locus of uniparental chromosome (Window 1 focal)
win2_focal = args.gg # Focal genotype of X/A chromosome (Window 2 focal)

#### Identify individuals to fill sterile/fertile groups
### Uni-X incompatibility
i=0
while len(sterile_ids_X) < number_sterile or len(fertile_f_ids_X) < number_fertile and len(fertile_nf_ids_X) < number_fertile:
	if i == int(df.shape[1]):
		break # Break if groups are not filled and individuals are exhausted
	x_loci = df.iat[x_inc,i]
	uni_loci = uni_windows[i]
	if uni_loci == Ymito_focal:   # If focal window 1
		if x_loci == win2_focal:   # If focal window 2
			prob = rd.randint(0,3) # 25% penatrence
			if prob == 0:
				if len(sterile_ids_X) < number_sterile: # Add to ID list if not equal to input sterile number
					sterile_ids_X.append(i)
				else:
					pass
			else: # If penatrence draw not success
				if rd.randrange(0,100) < args.bs: # Background sterility draw
					if len(sterile_ids_X) < number_sterile:
						sterile_ids_X.append(i)
					else:
						pass
				else:
					if len(fertile_f_ids_X) < number_fertile:
						fertile_f_ids_X.append(i)
					else:
						pass
		else: # If not focal window 2
			if rd.randrange(0,100) < args.bs:
				if len(sterile_ids_X) < number_sterile:
					sterile_ids_X.append(i)
				else:
					pass
			else:
				if len(fertile_f_ids_X) < number_fertile:
					fertile_f_ids_X.append(i)
	else: # If not focal window 1
		if rd.randrange(0,100) < args.bs:
			if len(sterile_ids_X) < number_sterile:
				sterile_ids_X.append(i)
			else:
				pass
		else:
			if len(fertile_nf_ids_X) < number_fertile:
				fertile_nf_ids_X.append(i)
			else:
				pass
	i = i + 1

### Initialize Breslow-Day cell count lists
###                  W2F                         W2NF
### Table1:        S      F       Table2:        S       F
###        W1F    bd1    bd2              W1F   bd5     bd6
###        W1NF   bd3    bd4              W1NF  bd7     bd8
###
bd1 = []
bd2 = []
bd3 = []
bd4 = []
bd5 = []
bd6 = []
bd7 = []
bd8 = []

### Count Breslow-Day cells
for index in sterile_ids_X:
	if uni_windows[index] == Ymito_focal:
		if df.iat[x_inc, index] == win2_focal: # Win1 F, Win2 F
			bd1.append(1)
		else: # Win1 F, Win2 NF
			bd5.append(1)
	else:
		if df.iat[x_inc, index] == win2_focal: # Win1 NF, Win2 F
			bd3.append(1)
		else: # Win1 NF, Win2 NF
			bd7.append(1)
for index in fertile_f_ids_X:
	if df.iat[x_inc, index] == win2_focal: # Win1 F, Win2 F
		bd2.append(1)
	else: # Win1 F, Win2 NF
		bd6.append(1)
for index in fertile_nf_ids_X:
	if df.iat[x_inc, index] == win2_focal: # Win1 NF, Win2 F
		bd4.append(1)
	else: # Win1 NF, Win2 NF
		bd8.append(1)
### Total counts
bd1 = len(bd1)
bd2 = len(bd2)
bd3 = len(bd3)
bd4 = len(bd4)
bd5 = len(bd5)
bd6 = len(bd6)
bd7 = len(bd7)
bd8 = len(bd8)

### Put counts into two 2x2 arrays
array_one = np.array([[bd1, bd2], [bd3, bd4]])
array_two = np.array([[bd5, bd6], [bd7, bd8]])

### Maximum odds filter
### Odds must be greatest in the comparison between sterile and fertile W1F,W2F group (O.R. = bd1/bd1 + bd2)
try:
	odds_one = bd1/(bd1 + bd2)
	odds_two = bd3/(bd3 + bd4)	
	odds_three = bd5/(bd5 + bd6)	
	odds_four = bd7/(bd7 + bd8)
except: # Catch divide by zero errors. Note: cases of 0/0 will default to an O.R = 0.50
	odds_one = (0.5 + bd1)/(bd1 + bd2 + 1)
	odds_two = (0.5 + bd3)/(bd3 + bd4 + 1)	
	odds_three = (0.5 + bd5)/(bd5 + bd6 + 1)	
	odds_four = (0.5 + bd7)/(bd7 + bd8 + 1)		

odds_list = [odds_one, odds_two, odds_three, odds_four]
max_odds = max(odds_list)
if max_odds == odds_one:
	array_one = array_one + 0.5 # 0.5 added as Fisher correction to avoid divide by zero
	array_two = array_two + 0.5
	tester = [array_one, array_two]
	results = sm.StratifiedTable(tester, shift_zeros=False) # Generate stratified table for input to BD test
	try:
		odds_test = results.test_equal_odds()  # Breslow-Day Test of equal odds (See: statsmodels.stats.contingency_tables.StratifiedTable.test_equal_odds)
		#odds_test = results.test_equal_odds(adjust=bool) # Tarone's adjustment alternate
		emp_pvalue_X = odds_test.pvalue
	except:
		emp_pvalue_X = 'NaN' # Return NaN float (no value) if BD test fails
else:
	emp_pvalue_X = 1 # Return 1 if odds filter fails

#### Uni-A incompatibility
i=0
while len(sterile_ids_A) < number_sterile or len(fertile_f_ids_A) < number_fertile and len(fertile_nf_ids_A) < number_fertile:
	if i == int(df.shape[1]):
		break
	a_loci = df.iat[a_inc,i]
	uni_loci = uni_windows[i]
	if uni_loci == Ymito_focal:   
		if a_loci == win2_focal:  
			prob = rd.randint(0,3)
			if prob == 0:
				if len(sterile_ids_A) < number_sterile:
					sterile_ids_A.append(i)
				else:
					pass
			else:
				if rd.randrange(0,100) < args.bs:
					if len(sterile_ids_A) < number_sterile:
						sterile_ids_A.append(i)
					else:
						pass
				else:
					if len(fertile_f_ids_A) < number_fertile:
						fertile_f_ids_A.append(i)
					else:
						pass
		else:
			if rd.randrange(0,100) < args.bs:
				if len(sterile_ids_A) < number_sterile:
					sterile_ids_A.append(i)
				else:
					pass
			else:
				if len(fertile_f_ids_A) < number_fertile:
					fertile_f_ids_A.append(i)
	else:
		if rd.randrange(0,100) < args.bs:
			if len(sterile_ids_A) < number_sterile:
				sterile_ids_A.append(i)
			else:
				pass
		else:
			if len(fertile_nf_ids_A) < number_fertile:
				fertile_nf_ids_A.append(i)
			else:
				pass
	i = i + 1

bd1 = []
bd2 = []
bd3 = []
bd4 = []
bd5 = []
bd6 = []
bd7 = []
bd8 = []

### Count Breslow-Day cells
for index in sterile_ids_A:
	if uni_windows[index] == Ymito_focal:
		if df.iat[a_inc, index] == win2_focal: # Win1 F, Win2 F
			bd1.append(1)
		else: # Win1 F, Win2 NF
			bd5.append(1)
	else:
		if df.iat[a_inc, index] == win2_focal: # Win1 NF, Win2 F
			bd3.append(1)
		else: # Win1 NF, Win2 NF
			bd7.append(1)
for index in fertile_f_ids_A:
	if df.iat[a_inc, index] == win2_focal: # Win1 F, Win2 F
		bd2.append(1)
	else: # Win1 F, Win2 NF
		bd6.append(1)
for index in fertile_nf_ids_A:
	if df.iat[a_inc, index] == win2_focal: # Win1 NF, Win2 F
		bd4.append(1)
	else: # Win1 NF, Win2 NF
		bd8.append(1)
bd1 = len(bd1)
bd2 = len(bd2)
bd3 = len(bd3)
bd4 = len(bd4)
bd5 = len(bd5)
bd6 = len(bd6)
bd7 = len(bd7)
bd8 = len(bd8)

array_one = np.array([[bd1, bd2], [bd3, bd4]])
array_two = np.array([[bd5, bd6], [bd7, bd8]])
### Maximum odds filter
try:
	odds_one = bd1/(bd1 + bd2)
	odds_two = bd3/(bd3 + bd4)	
	odds_three = bd5/(bd5 + bd6)	
	odds_four = bd7/(bd7 + bd8)
except:
	odds_one = (0.5 + bd1)/(bd1 + bd2 + 1)
	odds_two = (0.5 + bd3)/(bd3 + bd4 + 1)	
	odds_three = (0.5 + bd5)/(bd5 + bd6 + 1)	
	odds_four = (0.5 + bd7)/(bd7 + bd8 + 1)		

odds_list = [odds_one, odds_two, odds_three, odds_four]
max_odds = max(odds_list)
if max_odds == odds_one:
	array_one = array_one + 0.5 # 0.5 added as Fisher correction to avoid divide by zero
	array_two = array_two + 0.5
	tester = [array_one, array_two]
	results = sm.StratifiedTable(tester, shift_zeros=False)
	try:
		odds_test = results.test_equal_odds()  # Breslow-Day Test of equal odds (See: statsmodels.stats.contingency_tables.StratifiedTable.test_equal_odds)
		#odds_test = results.test_equal_odds(adjust=bool) # Tarone's adjustment alternate
		emp_pvalue_A = odds_test.pvalue
	except:
		emp_pvalue_A = 'NaN' # Return NaN float (no value) if BD test fails
else:
	emp_pvalue_A = 1 # Return 1 if odds filter fails

	
#### Null calculations


### Randomize sterile/fertile identity of samples
### X-A focal incompatibility 
null_ids_X = sterile_ids_X + fertile_f_ids_X + fertile_nf_ids_X # Put IDs used in expiremental test into list
rd.shuffle(null_ids_X) # Randomly shuffle the order of the list
null_sterile_ids_X = null_ids_X[0:number_sterile] # Assign first N number of IDs as sterile, where N=arg.s
null_fertile_ids_X = null_ids_X[number_sterile:]  # Assign remaining IDs as fertile, where len(null_fertile_ids_X)== arg.f

### A-A focal incompatibility
null_ids_A = sterile_ids_A + fertile_f_ids_A + fertile_nf_ids_A
rd.shuffle(null_ids_A)
null_sterile_ids_A = null_ids_A[0:number_sterile]
null_fertile_ids_A = null_ids_A[number_sterile:]

### Initialize count list for X-A (X) and A-A (A) incompatibility tests
bd1_X = []
bd2_X = []
bd3_X = []
bd4_X = []
bd5_X = []
bd6_X = []
bd7_X = []
bd8_X = []
bd1_A = []
bd2_A = []
bd3_A = []
bd4_A = []
bd5_A = []
bd6_A = []
bd7_A = []
bd8_A = []

### Uni-X window counts
len_X = list(range(0, X_end))
len_A = list(range(X_end, Chr3_end))
for window in len_X:
	f1 = 0
	f3 = 0
	f5 = 0
	f7 = 0
	for index in null_sterile_ids_X:
		if uni_windows[index] == Ymito_focal: # W1 F
			if df.iat[window, index] == win2_focal: # W2 F
				f1 += 1
			else: # W2 NF
				f5 += 1
		else: #W1 NF
			if df.iat[window, index] == win2_focal: # W2 F
				f3 += 1
			else: # W2 NF
				f7 += 1
	bd1_X.append(f1)
	bd3_X.append(f3)
	bd5_X.append(f5)
	bd7_X.append(f7)
for window in len_X:
	f2 = 0
	f4 = 0
	f6 = 0
	f8 = 0
	for index in null_fertile_ids_X:
		if uni_windows[index] == Ymito_focal: #W1 F
			if df.iat[window, index] == win2_focal: # W2 F
				f2 += 1
			else: # W2 NF
				f6 += 1
		else: # W1 NF
			if df.iat[window, index] == win2_focal: # W2 F:
				f4 += 1
			else: # W2 NF
				f8 += 1
	bd2_X.append(f2)
	bd4_X.append(f4)
	bd6_X.append(f6)
	bd8_X.append(f8)

### Uni-A window counts
for window in len_A:
	f1 = 0
	f3 = 0
	f5 = 0
	f7 = 0
	for index in null_sterile_ids_A:
		if uni_windows[index] == Ymito_focal: # W1 F
			if df.iat[window, index] == win2_focal: # W2 F
				f1 += 1
			else: # W2 NF
				f5 += 1
		else: #W1 NF
			if df.iat[window, index] == win2_focal: # W2 F
				f3 += 1
			else: # W2 NF
				f7 += 1
	bd1_A.append(f1)
	bd3_A.append(f3)
	bd5_A.append(f5)
	bd7_A.append(f7)
for window in len_A:
	f2 = 0
	f4 = 0
	f6 = 0
	f8 = 0
	for index in null_fertile_ids_A:
		if uni_windows[index] == Ymito_focal: # W1 F
			if df.iat[window, index] == win2_focal: # W2 F
				f2 += 1
			else: # W2 NF
				f6 += 1
		else: # W1 NF
			if df.iat[window, index] == win2_focal: # W2 F:
				f4 += 1
			else: # W2 NF
				f8 += 1
	bd2_A.append(f2)
	bd4_A.append(f4)
	bd6_A.append(f6)
	bd8_A.append(f8)
### Ensure equal number of comparisons between count groups
assert len(bd1_X) == len(bd2_X) == len(bd3_X) == len(bd4_X) == len(bd5_X) == len(bd6_X) == len(bd7_X) == len(bd8_X), "Error in BD group counting X"
assert len(bd1_A) == len(bd2_A) == len(bd3_A) == len(bd4_A) == len(bd5_A) == len(bd6_A) == len(bd7_A) == len(bd8_A), "Error in BD group counting A"

### Null Uni-X BD tests
#nan_count_X = [1]
p_values_X = []
for i in range(0, len(bd1_X)):
	array_one = np.array([[bd1_X[i], bd2_X[i]], [bd3_X[i], bd4_X[i]]])
	array_two = np.array([[bd5_X[i], bd6_X[i]], [bd7_X[i], bd8_X[i]]])
	try:
		odds_one = bd1_X[i]/(bd1_X[i] + bd2_X[i])
		odds_two = bd3_X[i]/(bd3_X[i] + bd4_X[i])
		odds_three = bd5_X[i]/(bd5_X[i] + bd6_X[i])
		odds_four = bd7_X[5,i]/(bd7_X[i] + bd8_X[i])
	except:
		odds_one = (0.5 + bd1_X[i])/(bd1_X[i] + bd2_X[i] + 1)
		odds_two = (0.5 + bd3_X[i])/(bd3_X[i] + bd4_X[i] + 1)
		odds_three = (0.5 + bd5_X[i])/(bd5_X[i] + bd6_X[i] + 1)
		odds_four = (0.5 + bd7_X[i])/(bd7_X[i] + bd8_X[i] + 1)
	odds_list = [odds_one, odds_two, odds_three, odds_four]
	max_odds = max(odds_list)
	if max_odds == odds_one:
		array_one = array_one + 0.5
		array_two = array_two + 0.5
		tester = [array_one, array_two]
		results = sm.StratifiedTable(tester, shift_zeros=False)
		try:
			odds_test = results.test_equal_odds()
			#odds_test = results.test_equal_odds(adjust=bool) # Tarone's adjustment alternate
			p_values_X.append(odds_test.pvalue)
		except:
			p_values_X.append(float('NaN'))
			#nan_count_X = []		
	else:
		p_values_X.append(1)

### Null Uni-A BD tests
#nan_count_A = [1]
p_values_A = []
for i in range(0, len(bd1_A)):
	array_one = np.array([[bd1_A[i], bd2_A[i]], [bd3_A[i], bd4_A[i]]])
	array_two = np.array([[bd5_A[i], bd6_A[i]], [bd7_A[i], bd8_A[i]]])
	try:
		odds_one = bd1_A[i]/(bd1_A[i] + bd2_A[i])
		odds_two = bd3_A[i]/(bd3_A[i] + bd4_A[i])
		odds_three = bd5_A[i]/(bd5_A[i] + bd6_A[i])
		odds_four = bd7_A[i]/(bd7_A[i] + bd8_A[i])
	except:
		odds_one = (0.5 + bd1_A[i])/(bd1_A[i] + bd2_A[i] + 1)
		odds_two = (0.5 + bd3_A[i])/(bd3_A[i] + bd4_A[i] + 1)
		odds_three = (0.5 + bd5_A[i])/(bd5_A[i] + bd6_A[i] + 1)
		odds_four = (0.5 + bd7_A[i])/(bd7_A[i] + bd8_A[i] + 1)
	odds_list = [odds_one, odds_two, odds_three, odds_four]
	max_odds = max(odds_list)
	if max_odds == odds_one:
		array_one = array_one + 0.5
		array_two = array_two + 0.5
		tester = [array_one, array_two]
		results = sm.StratifiedTable(tester, shift_zeros=False)
		try:
			odds_test = results.test_equal_odds()
			#odds_test = results.test_equal_odds(adjust=bool) # Tarone's adjustment alternate
			p_values_A.append(odds_test.pvalue)
		except:
			p_values_A.append(float('NaN'))
			#nan_count.append(1)		
	else:
		p_values_A.append(1)

### Converts NaN values to zeros, then replaces with 999 pvalue result to leave lowest pvalue unaffected by NaN results
p_values_X = np.nan_to_num(p_values_X, copy=False)
p_values_A = np.nan_to_num(p_values_A, copy=False) 
p_values_X = [x or 999 for x in p_values_X]
p_values_A = [x or 999 for x in p_values_A]

### Find the lowest null pvalue for Uni-X/A comparisons
lowest_p_X = str(min(p_values_X))
lowest_p_A = str(min(p_values_A))
emp_print_A = str(emp_pvalue_A)
emp_print_X = str(emp_pvalue_X)

### Output to file
### Prints one row, four columns in order:
### Uni-X tested pvalue, Uni-A tested pvalue, null Uni-X pvalue, null Uni-A pvalue
lister = open(args.o, 'w')
lister.write(emp_print_X + ',' + emp_print_A + ',' + lowest_p_X + ',' + lowest_p_A)
lister.close()

### Output to stdout
#nan_out_X = len(nan_count_X) - 1
#nan_out_A = len(nan_count_A) - 1
#print("There were " + str(nan_out_X) + " NaN values in null X comparisons") 
#print("There were " + str(nan_out_A) + " NaN values in null A comparisons")
