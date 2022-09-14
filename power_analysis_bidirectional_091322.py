#### Latest Update: September 13th, 2022 MJL
#### Power analysis simulator for XA and AA, two-locus recessive incompatibilities
#### Currently scripted to work with a 25% Penetrance of BDMI locus
#### Requires input have calls in order X -> Autosome

### Required Python libraries
import argparse
import numpy as np
import pandas as pd
import random as rd
import statsmodels.stats.contingency_tables as sm

### Input arguments
parser = argparse.ArgumentParser(description='Power Analysis Simulator, last update: 10/2021')
parser.add_argument('--i', help='Input File', required=True, type=str)
parser.add_argument('--o', help='Output File Prefix', required=True, type=str)
parser.add_argument('--s', help='Number of Sterile individuals', type=int, required=True)
parser.add_argument('--f', help='Number of Fertile individuals', type=int, required=True)
parser.add_argument('--bs', help='Background Sterility (As integer value of percent)', type=int, required=True)
parser.add_argument('--xa', help='Run X-A bidirectional scan', action='store_true')
parser.add_argument('--aa', help='Run A-A bidirectional scan', action='store_true')
parser.add_argument('--rev', help='Run Reverse Focal Genotype (default: window 1=0 and window2=2))', action='store_true')
args = parser.parse_args()
number_sterile = args.s
number_fertile = args.f
np.seterr(all='raise') #Catch-all for numpy error handling

### Chromosome window boundaries (zero-indexed)
## Adjust as needed
X_end = 545
Chr2_end = 1524
Chr3_end = 2579

### Load Array, Transpose (Inidividuals read in row-wise individual and converted to column-wise)
### Needed for current sibsam output format from sibsam_flies_2579windows_300_ind_Dec2021.txt
num_array = np.loadtxt(args.i, dtype=int)
num_array = num_array.transpose()
df = pd.DataFrame(num_array)

### Empirical test
def emp_test(s_id, f_id, test):
	### Define focal incompatibility
	if test == 'xa':
		inc_1 = x_inc
		inc_2 = a_inc
		global emp_pvalue_XA
	else:
		inc_1 = a_inc2
		inc_2 = a_inc3
		global emp_pvalue_AA
	
	### Initialize Breslow-Day cell count lists
	###                   W2F                            W2NF
	###              W1F   W1NF                 W1F    W1NF
	###        S    bd1     bd2              S   bd5     bd6
	###        F    bd3     bd4              F    bd7     bd8

	bd1 = 0 # S, W1F, W2F
	bd2 = 0 # S, W1NF, W2F
	bd3 = 0 # F, W1F, W2F
	bd4 = 0 # F, W1NF, W2F
	bd5 = 0 # S, W1F, W2NF
	bd6 = 0 # S, W1NF, W2NF
	bd7 = 0 # F, W1F, W2NF
	bd8 = 0 # F, W1NF, W2NF

	### Count Breslow-Day cells
	for index in s_id:
		if df.at[inc_1, index] == focal_1:
			if df.at[inc_2, index] == focal_2:
				bd1 += 1 #Sterile, W1F/W2F
			else:
				bd5 += 1 #Sterile, W1F/W2NF
		else:
			if df.at[inc_2, index] == focal_2:
				bd2 += 1 #Sterile, W1NF/WF2
			else:
				bd6 += 1 #Sterile, W1NF/W2NF
	for index in f_id:
		if df.at[inc_1, index] == focal_1:
			if df.at[inc_2, index] == focal_2:
				bd3 += 1 #Fertile, W1F/W2NF
			else:
				bd7 += 1 #Fertile, W1F/W2NF
		else:
			if df.at[inc_2, index] == focal_2:
				bd4 += 1 #Fertile W1NF/W2F
			else:
				bd8 += 1 #Fertile W1NF/W2NF

	### Arrange counts into array for BD test
	array_one = np.array([[bd1, bd2], [bd3,bd4]])
	array_two = np.array([[bd5, bd6], [bd7,bd8]])
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
			if test == 'xa':
				emp_pvalue_XA = odds_test.pvalue
			else:
				emp_pvalue_AA = odds_test.pvalue
		except:
			if test == 'xa':
				emp_pvalue_XA = -999 # Return 999 if BD test fails
			else:
				emp_pvalue_AA = -999
	else:
		if test == 'xa':
				emp_pvalue_XA = 999 # Return -999 if no test runs
		else:
				emp_pvalue_AA = 999	

### Null test
def null_test(s_id, f_id, test):
	### Define chromosome boundaries
	if test == 'xa':
		len_1 = list(range(0, X_end))
		len_2 = list(range(X_end, Chr3_end))
	else:
		len_1 = list(range(X_end, Chr2_end))
		len_2 = list(range(Chr2_end, Chr3_end))
	## Shuffle groups randomly
	null_ids = s_id + f_id
	rd.shuffle(null_ids)
	null_s_id = null_ids[0:number_sterile]
	null_f_id = null_ids[number_sterile:]

	bd1 = [] # S, W1F, W2F
	bd2 = [] # S, W1NF, W2F
	bd3 = [] # F, W1F, W2F
	bd4 = [] # F, W1NF, F2F
	bd5 = [] # S, W1F, W2NF
	bd6 = [] # S, W1NF, W2NF
	bd7 = [] # F, W1F, W2NF
	bd8 = [] # F, W1F, W2NF

	#Forward count
	for window_1 in len_1:
		for window_2 in len_2:
			b1 = 0
			b2 = 0
			b5 = 0
			b6 = 0
			for index in null_s_id: #Sterile
				if df.at[window_1, index] == focal_1: #W1F
					if df.at[window_2, index] == focal_2: #W2F
						b1 += 1
					else: #W2NF
						b5 += 1
				else: #W1NF
					if df.at[window_2, index] == focal_2: #W2F
						b2 += 1
					else: #W2NF
						b6 += 1
			bd1.append(b1)
			bd2.append(b2)
			bd5.append(b5)
			bd6.append(b6)
	for window_1 in len_1:
		for window_2 in len_2:
			b3 = 0
			b4 = 0
			b7 = 0
			b8 = 0
			for index in null_f_id: #Fertile
				if df.at[window_1, index] == focal_1:
					if df.at[window_2, index] == focal_2:
						b3 += 1
					else:
						b7 += 1
				else:
					if df.at[window_2, index] == focal_2:
						b4 += 1
					else:
						b8 += 1
			bd3.append(b3)
			bd4.append(b4)
			bd7.append(b7)
			bd8.append(b8)
	##Reverse Count
	for window_1 in len_2:
		for window_2 in len_1:
			b1 = 0
			b2 = 0
			b5 = 0
			b6 = 0
			for index in null_s_id: #Sterile
				if df.at[window_1, index] == focal_1: #W1F
					if df.at[window_2, index] == focal_2: #W2F
						b1 += 1
					else: #W2NF
						b5 += 1
				else: #W1NF
					if df.at[window_2, index] == focal_2: #W2F
						b2 += 1
					else: #W2NF
						b6 += 1
			bd1.append(b1)
			bd2.append(b2)
			bd5.append(b5)
			bd6.append(b6)
	for window_1 in len_2:
		for window_2 in len_1:
			b3 = 0
			b4 = 0
			b7 = 0
			b8 = 0
			for index in null_f_id: #Fertile
				if df.at[window_1, index] == focal_1:
					if df.at[window_2, index] == focal_2:
						b3 += 1
					else:
						b7 += 1
				else:
					if df.at[window_2, index] == focal_2:
						b4 += 1
					else:
						b8 += 1
			bd3.append(b3)
			bd4.append(b4)
			bd7.append(b7)
			bd8.append(b8)
	### Sanity Assert
	assert len(bd1) == len(bd2) == len(bd3) == len(bd4) == len(bd5) == len(bd6) == len(bd7) == len(bd8)
	### Null BD test
	p_values = []
	error_count = 0
	end_range = len(bd1)
	for i in range(0, end_range):
		array_one = np.array([[bd1[i], bd2[i]], [bd3[i], bd4[i]]])
		array_two = np.array([[bd5[i], bd6[i]], [bd7[i], bd8[i]]])
		try:
			odds_one = bd1[i]/(bd1[i] + bd2[i])
			odds_two = bd3[i]/(bd3[i] + bd4[i])
			odds_three = bd5[i]/(bd5[i] + bd6[i])
			odds_four = bd7[i]/(bd7[i] + bd8[i])
		except:
			odds_one = (0.5 + bd1[i])/(bd1[i] + bd2[i] + 1)
			odds_two = (0.5 + bd3[i])/(bd3[i] + bd4[i] + 1)
			odds_three = (0.5 + bd5[i])/(bd5[i] + bd6[i] + 1)
			odds_four = (0.5 + bd7[i])/(bd7[i] + bd8[i] + 1)
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
				p_values.append(odds_test.pvalue)
			except:
				p_values.append(999)
				print("Encountered BD test error at test: " + str(i))
				error_count += 1
		else:
			p_values.append(999)
	print(str(error_count))
	#Get lowest p_value
	if test == 'xa':
		global null_pvalue_XA
		null_pvalue_XA = str(min(p_values))
	else:
		global null_pvalue_AA
		null_pvalue_AA = str(min(p_values))


### Define focal genotype
if args.rev == True:
	focal_1 = 2
	focal_2 = 0
else:
	focal_1 = 0
	focal_2 = 2

### X-Autosome incompatibility testing
if args.xa == True:
	sterile_ids_XA = []
	fertile_ids_XA = []

	#Randomly define X-A incompatibility
	x_inc = rd.randint(1,(X_end-1))
	a_inc = rd.randint(X_end, (Chr3_end-1))
	
	### Score sterile and fertile individuals for X-A incompatibility
	i=0
	while len(sterile_ids_XA) < number_sterile or len(fertile_ids_XA) < number_fertile:
		if i == int(df.shape[1]):
			print("Number of individuals in input exhausted before filling sterile/fertile groups.")
			sys.exit("Exiting program...") #Exit if groups are not filled after individuals are exhausted
		site_1 = df.at[x_inc, i]
		site_2 = df.at[a_inc, i]
		if site_1 == focal_1 and site_2 == focal_2: #If focal window1/window2
			 if rd.randint(0,3) == 0: #25% penetrance
			 	if len(sterile_ids_XA) < number_sterile:
			 		sterile_ids_XA.append(i)
			 	else:
			 		pass
			 else:
			 	if rd.randrange(0,100) < args.bs: #Background sterility
			 		if len(sterile_ids_XA) < number_sterile:
			 			sterile_ids_XA.append(i)
			 		else:
			 			pass
			 	else:
			 		if len(fertile_ids_XA) < number_fertile:
			 			fertile_ids_XA.append(i)
			 		else:
			 			pass
		else:
			if len(fertile_ids_XA) < number_fertile:
				fertile_ids_XA.append(i)
			else:
				pass
		i = i + 1
	
	### Call tests
	emp_test(sterile_ids_XA, fertile_ids_XA, 'xa')
	null_test(sterile_ids_XA, fertile_ids_XA, 'xa')
	
	### Print to output
	out_name = args.o + '_xa_power_pvalues_bidirectional.csv'
	lister = open(out_name, 'w')
	lister.write(str(emp_pvalue_XA) + ',' + str(null_pvalue_XA))
	lister.close()

## Autosome-Autosome incompatibility testing
if args.aa == True:
	sterile_ids_AA = []
	fertile_ids_AA = []
	
	### Randomly assign A-A incompatibility
	a_inc2 = rd.randint(X_end, (Chr2_end-1))
	a_inc3 = rd.randint(Chr2_end, (Chr3_end-1))

	### Score sterile and fertile individuals for X-A incompatibility
	i=0
	while len(sterile_ids_AA) < number_sterile or len(fertile_ids_AA) < number_fertile:
		if i == int(df.shape[1]):
			print("Number of individuals in input exhausted before filling sterile/fertile groups.")
			sys.exit("Exiting program...") #Exit if groups are not filled after individuals are exhausted
		site_1 = df.at[a_inc2, i]
		site_2 = df.at[a_inc3, i]
		if site_1 == focal_1 and site_2 == focal_2: #If focal window1/window2
			 if rd.randint(0,3) == 0: #25% penetrance
			 	if len(sterile_ids_AA) < number_sterile:
			 		sterile_ids_AA.append(i)
			 	else:
			 		pass
			 else:
			 	if rd.randrange(0,100) < args.bs: #Background sterility
			 		if len(sterile_ids_AA) < number_sterile:
			 			sterile_ids_AA.append(i)
			 		else:
			 			pass
			 	else:
			 		if len(fertile_ids_AA) < number_fertile:
			 			fertile_ids_AA.append(i)
			 		else:
			 			pass
		else:
			if len(fertile_ids_AA) < number_fertile:
				fertile_ids_AA.append(i)
			else:
				pass
		i = i + 1
	
	### Call tests
	emp_test(sterile_ids_AA, fertile_ids_AA, 'aa')
	null_test(sterile_ids_AA, fertile_ids_AA, 'aa')
	
	### Print to output
	out_name = args.o + '_aa_power_pvalues_bidirectional.csv'
	lister = open(out_name, 'w')
	lister.write(str(emp_pvalue_AA) + ',' + str(null_pvalue_AA))
	lister.close()
