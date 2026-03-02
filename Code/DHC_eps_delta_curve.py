import numpy as np
import csv, os
import itertools
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
import math
import mpmath as mp
from mpmath import mpf, exp, findroot, almosteq
from collections import Counter
warnings.filterwarnings("ignore", category=RuntimeWarning)
import time
import sys
from fractions import Fraction
import argparse

# ------------------ parse args early (single parser; includes --k, --age, --n) ------------------
parser = argparse.ArgumentParser(description="index of directed path along the US Census Bureau allocation paths")
parser.add_argument("--k", type=int, default=13, help="Path number k")
arg = parser.parse_args()
path_num = int(arg.k)

def format_func(value, tick_number):
    if value < 10:
        return f"{value:.2f}"
    elif value < 100:
        return f"{value:.1f}"
    else:
        return f"{int(value)}"

start_time = time.time()

############# Precision #############
accuracy = 50
mp.dps = accuracy

mp.mp.dps = 50  # Set 50 decimal places of precision

## Test Machine Error

s = 1.0

for k in range(1, 101):
    s /= 2.0
    if (1.0 + s) <= 1.0:
        s *= 2.0
        print('Machine Epsilon by default (not privacy budget):', f"k={k-1}, eps={s}")
        break

s = mp.mpf('1.0')

for k in range(1, 1001):
    s /= mp.mpf('2.0')
    if (mp.mpf('1.0') + s) <= mp.mpf('1.0'):
        s *= mp.mpf('2.0')
        print('Machine Epsilon by mpmath (not privacy budget):', f"k={k-1}, eps={s}")
        break

## one record may change from path_num_1 to path_num_2 due to the choice of sensitivity
path_num_1 = path_num
path_num_2 = path_num

## Preparation

### Parameters
zero = mp.mpf('0')
one = mp.mpf('1')
two = mp.mpf('2')
three = mp.mpf('3')
four = mp.mpf('4')
five = mp.mpf('5')
six = mp.mpf('6')
seven = mp.mpf('7')
eight = mp.mpf('8')
nine = mp.mpf('9')

ten = mp.mpf('10')
eleven = mp.mpf('11')
twelve = mp.mpf('12')
thirteen = mp.mpf('13')
fourteen = mp.mpf('14')
fifteen = mp.mpf('15')
sixteen = mp.mpf('16')
seventeen = mp.mpf('17')
eighteen = mp.mpf('18')
nineteen = mp.mpf('19')

twenty = mp.mpf('20')
twenty_one = mp.mpf('21')
twenty_two = mp.mpf('22')
twenty_three = mp.mpf('23')
twenty_four = mp.mpf('24')
twenty_five = mp.mpf('25')
twenty_six = mp.mpf('26')
twenty_seven = mp.mpf('27')
twenty_eight = mp.mpf('28')
twenty_nine = mp.mpf('29')

thirty_two = mp.mpf('32')
forty_five = mp.mpf('45')
sixty_four = mp.mpf('64')
two_over_45 = two / forty_five


ten = mp.mpf('10')
hundred = mp.mpf('100')
one_fifty = mp.mpf('150')
thousand = mp.mpf('1000')

def float_to_mpf(x):
    return mp.mpf(x)

def fraction_to_mpf(frac):
    """Convert a Fraction to an mpmath.mpf with high precision."""
    num = mp.mpf(frac.numerator)
    denom = mp.mpf(frac.denominator)
    return num / denom

### read allocation path
file_path_1 = f"allocation_path/dhc_allocation_path_{path_num_1}.csv"
df_1 = pd.read_csv(file_path_1)

allocation_run_1 = []
for col in df_1.columns:
    allocation_run_1.extend(df_1[col].tolist())

file_path_2 = f"allocation_path/dhc_allocation_path_{path_num_2}.csv"
df_2 = pd.read_csv(file_path_2)

allocation_run_2 = []
for col in df_2.columns:
    allocation_run_2.extend(df_2[col].tolist())

allocation_run = [entry for entry in allocation_run_1 if entry != '0/1'] + \
                 [entry for entry in allocation_run_2 if entry != '0/1']

allocation_run[:22], len(allocation_run)

from fractions import Fraction

allocation_run_sorted = sorted(
    allocation_run,
    key=lambda x: Fraction(x.strip("'"))
)

allocation_run = allocation_run_sorted
allocation_run[:22], len(allocation_run)


rho_1 = [Fraction(s) for s in allocation_run_1]
rho_2 = [Fraction(s) for s in allocation_run_2]
if rho_1 != rho_2:
    print("Incorrect Allocations")


## Determine Parameters

### Delta

Delta = mp.mpf('1e-36')

### rho

rho_path = [Fraction(s) for s in allocation_run]
rho = sum(rho_path)/2

### sigma2_list

sigma2_list = [1/val for val in rho_path]

### psi2

psi2 = 1/rho

### a_list

a_list = [val / rho for val in rho_path]

### L

# Extract the denominators of all percentage
denominators = [frac.denominator for frac in a_list]
numerators = [frac.numerator for frac in a_list]

# Compute LCM of all denominators using math.lcm
lcm_value = math.lcm(*denominators)
gcd_value = math.gcd(*numerators)

L = Fraction(lcm_value, gcd_value)

a_list_L = [val * L for val in a_list]
# print(a_list_L)

### m

m = len(a_list)

### convert to numerical

rho = fraction_to_mpf(rho)
a_list = [fraction_to_mpf(frac) for frac in a_list]
a_list_L = [fraction_to_mpf(frac) for frac in a_list_L]
m = float_to_mpf(m)
L = fraction_to_mpf(L)
psi2 = fraction_to_mpf(psi2)
sigma2_list = [fraction_to_mpf(frac) for frac in sigma2_list]

print("tolerance:", Delta, "rho", rho, "L:", L)


# Numerical Computation

### U

summ = mp.mpf('0')
for i in range(len(a_list)):
    summ += a_list[i]**2 * sigma2_list[i]

U = mp.ceil( L * mp.sqrt(- two * summ * mp.log(Delta / eight)) )

# optional only for nice expression of U, Find the smallest multiple of 100000 greater than U
step = mp.mpf('250000')
U = step * mp.ceil(U / step)

print(f"U must be larger than", U)

### N

summ = mp.mpf('0')
for i in range(len(a_list)):
    summ += a_list[i]**2 * sigma2_list[i]
log_term = mp.sqrt( - two  * summ * mp.log(Delta / eight) )

# ==== Final RHS expression ====
N = mp.ceil( U + L * log_term )

step = mp.mpf('250000')
N = step * mp.ceil(N / step)

print(f"N must be larger than", N)

## Characteristic Functions

# Step 1: Get unique sigma2 values and their first occurrence indices
sigma2_index_map = {}
sigma2_unique_list = []
index_unique_list = []

for idx, sigma2 in enumerate(sigma2_list):
    if sigma2 not in sigma2_index_map:
        sigma2_index_map[sigma2] = len(sigma2_unique_list)
        sigma2_unique_list.append(sigma2)
        index_unique_list.append(idx)

# Step 2: Count appearances of each unique sigma2 value
appearance_counter = Counter(sigma2_list)
appearance_list = [appearance_counter[sigma2] for sigma2 in sigma2_unique_list]


print(appearance_counter)
print()
print(appearance_list)

############# Characteristic Functions #############
def char_func_num(sigma_i2, a_i_L, zeta):
    # Precompute constants
    trunc = int(mp.ceil(twenty_four * two * sigma_i2) + one)

    denom_coeff = one / (two * sigma_i2)

    # Numerator sum
    numerator = one
    for x in range(1, trunc + 1):
        x_mpf = mp.mpf(x)
        numerator += two * mp.cos(a_i_L * zeta * x_mpf) * mp.exp(-x_mpf**2 * denom_coeff)

    return numerator


def char_func_num_prod(zeta):
    product = mp.mpf("1.0")
    for i in range(len(sigma2_unique_list)):
        sigma_i2 = sigma2_unique_list[i]
        a_i_L = a_list_L[index_unique_list[i]]
        char_val = char_func_num(sigma_i2, a_i_L, zeta)
        product *= char_val ** appearance_list[i]

    return product


char_func_denom_prod = char_func_num_prod(zero)


def char_func_values(zeta):
    values = []
    for i in range(len(sigma2_unique_list)):
        sigma_i2 = sigma2_unique_list[i]
        a_i_L = a_list_L[index_unique_list[i]]
        char_num_val = char_func_num(sigma_i2, a_i_L, zeta)
        char_denom_val = char_func_num(sigma_i2, a_i_L, zero)
        values.append( (char_num_val / char_denom_val) ** appearance_list[i] )

    return values


############# Weight Function #############
def weight_first(zeta, t_0, L):
    lower_1 = mp.ceil(t_0 * L)
    stable_cut = Delta
    if abs(zeta) >= stable_cut:
        first = mp.cos(lower_1 * zeta) + mp.cos(U * zeta)
        factor = mp.cos(zeta / two) / mp.sin(zeta / two)
        second = mp.sin(U * zeta) - mp.sin(lower_1 * zeta)
        return (one / two) * (one / mp.pi) * (first + factor * second)

    elif 0 < abs(zeta) < stable_cut:
        first = mp.cos(lower_1 * zeta) + mp.cos(U * zeta)
        factor = mp.cos(zeta / two) / (zeta / two)
        second = mp.sin(U * zeta) - mp.sin(lower_1 * zeta)
        return (one / two) * (one / mp.pi) * (first + factor * second)

    else:
        return (one / mp.pi) * (U - lower_1 + one)

def weight_second(zeta, T_0, L):
    lower_2 = mp.ceil(T_0 * L)
    stable_cut = Delta
    if abs(zeta) >= stable_cut:
        first = mp.cos(lower_2 * zeta) + mp.cos(U * zeta)
        factor = mp.cos(zeta / two) / mp.sin(zeta / two)
        second = mp.sin(U * zeta) - mp.sin(lower_2 * zeta)
        return (one / two) * (one / mp.pi) * (first + factor * second)

    elif 0 < abs(zeta) < stable_cut:
        first = mp.cos(lower_2 * zeta) + mp.cos(U * zeta)
        factor = mp.cos(zeta / two) / (zeta / two)
        second = mp.sin(U * zeta) - mp.sin(lower_2 * zeta)
        return (one / two) * (one / mp.pi) * (first + factor * second)

    else:
        return (one / mp.pi) * (U - lower_2 + one)

############# Integrant #############
def char_first(zeta):
    return char_func_num_prod(zeta) / char_func_denom_prod  # one/two here because even extension

def char_second(zeta):
    return mp.exp(eps_zcdp) * char_func_num_prod(zeta) / char_func_denom_prod


## Determine Threshold

results = []  # Initialize a list to store results

for i in range(len(sigma2_unique_list)):
    sigma_i2 = sigma2_unique_list[i]
    a_i_L = a_list_L[index_unique_list[i]]

    half_period = mp.pi / a_i_L
    zeta_list = np.linspace(0, 1, 10000) * half_period

    # Calculate char_denom_val once outside the loop for efficiency
    char_denom_val = char_func_num(sigma_i2, a_i_L, zero)

    for zeta in zeta_list:
        zeta_mpf = mp.mpf(zeta)

        char_num_val = char_func_num(sigma_i2, a_i_L, zeta_mpf)

        char_single = (char_num_val / char_denom_val) ** appearance_list[i]

        # Check stopping criterion
        if char_single < Delta/U:
            # Use formatted string for clearer output
            result_dict = {
                "sigma2": sigma_i2,
                "a_i_L": a_i_L,
                "zeta": zeta,
                "char_single": char_single
            }

            results.append(result_dict)

            # Print nicely formatted message
            print(f"Stopped at sigma2: {result_dict['sigma2']}, "
                  f"a_i_L: {result_dict['a_i_L']}, "
                  f"zeta: {result_dict['zeta']}, "
                  f"char_single: {result_dict['char_single']}")

            break  # break inner loop upon condition met

# Save the results to CSV for later analysis
df_results = pd.DataFrame(results)
df_results.to_csv(f'results_DHC_v0/eps_delta_curve/tmp/char_func_half_period_threshold_path_{path_num_1}_to_{path_num_2}.csv', index=False)

# Load previously computed data
df_results = pd.read_csv(f'results_DHC_v0/eps_delta_curve/tmp/char_func_half_period_threshold_path_{path_num_1}_to_{path_num_2}.csv')

# Ensure sigma2 values are precisely handled as mpf
df_results['sigma2'] = df_results['sigma2'].apply(mp.mpf)
df_results['zeta'] = df_results['zeta'].apply(mp.mpf)

midpoints = {}
zeta_interval = {}

for i in range(len(sigma2_unique_list)):
    sigma_i2 = sigma2_unique_list[i]
    a_i_L = a_list_L[index_unique_list[i]]
    half_period = mp.pi / a_i_L

    # Check presence with full precision
    matched_rows = df_results[df_results['a_i_L'] == int(a_i_L)]

    if matched_rows.empty:
        print(f"{sigma_i2} is not found")
        continue  # Skip if sigma_i2 not found

    print(f"Processing sigma2: {sigma_i2}, a_i_L: {a_i_L}")

    zeta_val = matched_rows['zeta'].iloc[0]

    # Step 1: Compute midpoints array
    two_pi_over_a = two * mp.pi / a_i_L
    midpoint = [two_pi_over_a * mp.mpf(j) for j in range(int(mp.floor(a_i_L)) + 1)]

    # Step 2: Keep only midpoints in [0, π]
    filtered_midpoints = [x for x in midpoint if zero <= x <= mp.pi]
    midpoints[i] = filtered_midpoints

    # Save zeta
    zeta_interval[i] = zeta_val

# Step 1: Sort the midpoints by the length of their value list
sorted_keys = sorted(midpoints.keys(), key=lambda k: len(midpoints[k]))

# Step 2: Reconstruct midpoints and zeta_interval with sorted keys
midpoints_sorted = {new_idx: midpoints[k] for new_idx, k in enumerate(sorted_keys)}
zeta_interval_sorted = {new_idx: zeta_interval[k] for new_idx, k in enumerate(sorted_keys)}

# overwrite the originals
midpoints = midpoints_sorted
zeta_interval = zeta_interval_sorted

# Number of intervals
n_intervals = int(N)
print("Number of intervals:", n_intervals)

# Output CSV path
csv_path = f"results_DHC_v0/eps_delta_curve/threshold/char_fnc_threshold_point_path_{path_num_1}_to_{path_num_2}.csv"

# Open CSV and write header
with open(csv_path, mode='w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(['zeta'])  # header

    for idx in range(n_intervals):
        i = mp.mpf(idx)
        zeta = two * mp.pi * i / n_intervals

        # Check for each j
        passed_all_j = True
        for j in midpoints:
            midpoint_list = midpoints[j]
            zeta_thresh = zeta_interval[j]

            # Find the closest midpoint to zeta
            closest = min(midpoint_list, key=lambda x: abs(x - zeta))

            # Check if abs(zeta - closest) < zeta_thresh
            if abs(zeta - closest) >= zeta_thresh:
                passed_all_j = False
                break  # early exit if any j fails

        if passed_all_j:
            writer.writerow([str(zeta)])  # only write valid zeta

n_intervals = int(N)
# print("Number of intervals:", n_intervals)

# Input and output paths
input_csv = f"results_DHC_v0/eps_delta_curve/threshold/char_fnc_threshold_point_path_{path_num_1}_to_{path_num_2}.csv"
output_csv = f"results_DHC_v0/eps_delta_curve/char_fnc_eval/first_char_fnc_eval_path_{path_num_1}_to_{path_num_2}.csv"

# Step 1: Read all zeta values
zeta_list = []
with open(input_csv, mode='r') as infile:
    reader = csv.reader(infile)
    next(reader)  # skip header
    for row in reader:
        zeta_list.append(mp.mpf(row[0]))

# Step 3: Write result row by row with zeta, increment, cumulative sum
with open(output_csv, mode='w', newline='') as outfile:
    writer = csv.writer(outfile)
    writer.writerow(['zeta', 'char_fnc'])  # header

    for zeta in zeta_list:
        print("zeta:", zeta)
        char = char_first(zeta)
        print("char_fnc:", char)
        writer.writerow([str(zeta), str(char)])
        print("*"*45)


delta_str_list = [
    '1e-1', '5e-1',
    '1e-2', '5e-2',
    '1e-3', '5e-3',
    '1e-4', '5e-4',
    '1e-5', '5e-5',
    '1e-6', '5e-6',
    '1e-7', '5e-7',
    '1e-8', '5e-8',
    '1e-9', '5e-9',
    '1e-10', '5e-10',
    '1e-11', '5e-11',
    '1e-12', '5e-12',
    '1e-13', '5e-13',
    '1e-14', '5e-14',
    '1e-15', '5e-15',
    '1e-16', '5e-16',
    '1e-17', '5e-17',
    '1e-18', '5e-18',
    '1e-19', '5e-19',
    '1e-20', '5e-20',
    '1e-21', '5e-21',
    '1e-22', '5e-22',
    '1e-23', '5e-23',
    '1e-24', '5e-24',
    '1e-25', '5e-25'
]


for delta_str in delta_str_list:
    print("delta:", delta_str)

    """## delta"""

    delta = mp.mpf(delta_str)

    """## eps_zcdp"""

    eps_zcdp = rho + two * mp.sqrt(- rho * mp.log(delta))
    print("eps_zcdp:", eps_zcdp)

    """## T(eps), t(eps)"""

    t_0 = psi2 * eps_zcdp - sum(a_list) / two
    T_0 = psi2 * eps_zcdp + sum(a_list) / two
    print("t_0, T_0:", t_0, T_0)

    # Input and output paths
    input_csv = f"results_DHC_v0/eps_delta_curve/threshold/char_fnc_threshold_point_path_{path_num_1}_to_{path_num_2}.csv"
    output_csv = f"results_DHC_v0/eps_delta_curve/weight_eval/eps_specific_delta_{delta_str}_first_weight_path_{path_num_1}_to_{path_num_2}.csv"

    # Step 1: Read all zeta values
    zeta_list = []
    with open(input_csv, mode='r') as infile:
        reader = csv.reader(infile)
        next(reader)  # skip header
        for row in reader:
            zeta_list.append(mp.mpf(row[0]))

    # Step 3: Write result row by row with zeta, increment, cumulative sum
    with open(output_csv, mode='w', newline='') as outfile:
        writer = csv.writer(outfile)
        writer.writerow(['zeta', 'weight'])  # header

        for zeta in zeta_list:
            # print("zeta:", zeta)
            weight = weight_first(zeta, t_0, L)
            # print("weight:", weight)
            writer.writerow([str(zeta), str(weight)])
            # print("*"*45)


    # Input and output paths
    input_csv = f"results_DHC_v0/eps_delta_curve/threshold/char_fnc_threshold_point_path_{path_num_1}_to_{path_num_2}.csv"
    output_csv = f"results_DHC_v0/eps_delta_curve/weight_eval/eps_specific_delta_{delta_str}_second_weight_path_{path_num_1}_to_{path_num_2}.csv"

    # Step 1: Read all zeta values
    zeta_list = []
    with open(input_csv, mode='r') as infile:
        reader = csv.reader(infile)
        next(reader)  # skip header
        for row in reader:
            zeta_list.append(mp.mpf(row[0]))


    # Step 3: Write result row by row with zeta, increment, cumulative sum
    with open(output_csv, mode='w', newline='') as outfile:
        writer = csv.writer(outfile)
        writer.writerow(['zeta', 'weight'])  # header

        for zeta in zeta_list:
            # print("zeta:", zeta)
            weight = weight_second(zeta, T_0, L)
            # print("weight:", weight)
            writer.writerow([str(zeta), str(weight)])
            # print("*"*45)


# Load character function and weights
first_char_csv = f"results_DHC_v0/eps_delta_curve/char_fnc_eval/first_char_fnc_eval_path_{path_num_1}_to_{path_num_2}.csv"
with open(first_char_csv, newline='') as f:
    reader = csv.DictReader(f)
    char_vals = [mp.mpf(row['char_fnc']) for row in reader]

# Prepare output storage
output_rows = []

# Loop over delta values
for delta_str in delta_str_list:
    # delta
    delta = mp.mpf(delta_str)

    # eps_zcdp
    eps_zcdp = rho + two * mp.sqrt(-rho * mp.log(delta))
    print("eps_zcdp:", eps_zcdp)
    print("delta_str:", delta_str)
    
    # Read weight CSV for this delta
    first_weight_csv = f"results_DHC_v0/eps_delta_curve/weight_eval/eps_specific_delta_{delta_str}_first_weight_path_{path_num_1}_to_{path_num_2}.csv"
    with open(first_weight_csv, newline='') as f:
        reader = csv.DictReader(f)
        weight_vals = [mp.mpf(row['weight']) for row in reader]

    # Compute integrand = char * weight * 1/2
    tobeint_first = [(c * w * one / two) for c, w in zip(char_vals, weight_vals)]

    # Compute increment = Δzeta * tobeint_first
    delta_zeta = two * mp.pi / N
    increment = [delta_zeta * val for val in tobeint_first]

    # Compute correction
    correction = delta_zeta * tobeint_first[0]

    # Compute final integral result
    result_first = two * mp.fsum(increment) - correction
    delta_fdp = result_first

    print("delta_fdp_first_term:", delta_fdp)

    # Store result
    output_rows.append([delta_str, mp.nstr(eps_zcdp, n=50), mp.nstr(delta_fdp, n=50)])

# Output file path
output_csv = f"results_DHC_v0/eps_delta_curve/eps_delta_results/eps_delta_curve_path_{path_num_1}_to_{path_num_2}.csv"
os.makedirs(os.path.dirname(output_csv), exist_ok=True)

# Save to CSV
with open(output_csv, mode='w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['delta_zcdp', 'epsilon', 'delta_fdp_first_term'])
    writer.writerows(output_rows)

print("✅ Summary saved")

input_char_path = f"results_DHC_v0/eps_delta_curve/char_fnc_eval/first_char_fnc_eval_path_{path_num_1}_to_{path_num_2}.csv"

# Step 1: Read zeta and char_fnc values
zeta_vals = []
char_vals = []
with open(input_char_path, mode='r') as infile:
    reader = csv.reader(infile)
    next(reader)  # skip header
    for row in reader:
        zeta_vals.append(mp.mpf(row[0]))
        char_vals.append(mp.mpf(row[1]))

# Loop over delta values
for delta_str in delta_str_list:
    # delta
    delta = mp.mpf(delta_str)

    # eps_zcdp
    eps_zcdp = rho + two * mp.sqrt(-rho * mp.log(delta))

    # Input and output file paths for the second stage
    output_char_path = f"results_DHC_v0/eps_delta_curve/char_fnc_eval/second_char_fnc_eval_delta_{delta_str}_path_{path_num_1}_to_{path_num_2}.csv"

    # Step 2: Apply exp(eps_zcdp) scaling
    exp_eps = mp.exp(eps_zcdp)
    scaled_char_vals = [val * exp_eps for val in char_vals]

    # Step 3: Write output to new CSV
    with open(output_char_path, mode='w', newline='') as outfile:
        writer = csv.writer(outfile)
        writer.writerow(['zeta', 'char_fnc'])  # header
        for zeta, scaled_char in zip(zeta_vals, scaled_char_vals):
            writer.writerow([str(zeta), str(scaled_char)])

print(f"Saved updated char_fnc values with exp(eps_zcdp) scaling to:\n{output_char_path}")

# Inputs
n_intervals = int(N)

# Prepare list to store second term results
second_term_results = []

# Loop over delta values
for delta_str in delta_str_list:
    # delta
    delta = mp.mpf(delta_str)

    # eps_zcdp
    eps_zcdp = rho + two * mp.sqrt(-rho * mp.log(delta))
    print("eps_zcdp:", eps_zcdp)
    print("delta_str:", delta_str)
    
    # Load character function (second term)
    second_char_csv = f"results_DHC_v0/eps_delta_curve/char_fnc_eval/second_char_fnc_eval_delta_{delta_str}_path_{path_num_1}_to_{path_num_2}.csv"
    with open(second_char_csv, newline='') as f:
        reader = csv.DictReader(f)
        char_vals = [mp.mpf(row['char_fnc']) for row in reader]

    # Load corresponding second weight
    second_weight_csv = f"results_DHC_v0/eps_delta_curve/weight_eval/eps_specific_delta_{delta_str}_second_weight_path_{path_num_1}_to_{path_num_2}.csv"
    with open(second_weight_csv, newline='') as f:
        reader = csv.DictReader(f)
        weight_vals = [mp.mpf(row['weight']) for row in reader]

    # Compute integrand
    tobeint_second = [(c * w * one / two) for c, w in zip(char_vals, weight_vals)]

    # Integration
    delta_zeta = two * mp.pi / N
    increment = [delta_zeta * val for val in tobeint_second]
    correction = delta_zeta * tobeint_second[0]
    result_second = two * mp.fsum(increment) - correction

    print("delta_fdp_second_term:", result_second)

    # Store the result
    second_term_results.append(mp.nstr(result_second, n=50))

# === Read existing CSV with full precision ===
df = pd.read_csv(output_csv, dtype=str)  # read all as string to preserve digits

# Add new column using 50-digit string formatting
df['delta_fdp_second_term'] = [mp.nstr(mpf(x), n=50) for x in second_term_results]

# Save CSV back with full-precision values as strings
df.to_csv(output_csv, index=False)
print(f"✅ Updated and saved with full precision: {output_csv}")

# Ensure both columns are strings (if read from CSV)
df['delta_fdp_first_term'] = df['delta_fdp_first_term'].astype(str)
df['delta_fdp_second_term'] = df['delta_fdp_second_term'].astype(str)

# Safely convert each row to mpf and compute the difference
df['delta_fdp'] = df.apply(
    lambda row: mp.nstr(mp.mpf(row['delta_fdp_first_term']) - mp.mpf(row['delta_fdp_second_term']), n=50),
    axis=1
)

# If your DataFrame is named df and delta_zcdp is in string format:
df['delta_zcdp_numeric'] = df['delta_zcdp'].astype(float)

# Sort by delta_zcdp value (ascending)
df = df.sort_values(by='delta_zcdp_numeric', ascending=False).reset_index(drop=True)

# Drop helper column
df.drop(columns=['delta_zcdp_numeric'], inplace=True)

# Save with full precision
df.to_csv(output_csv, index=False)
print("✅ delta_fdp saved with full 50-digit precision.")

end_time = time.time()
elapsed_time = (end_time - start_time)/3600
print(f"\n✅ Total runtime: {elapsed_time:.2f} Hours")
