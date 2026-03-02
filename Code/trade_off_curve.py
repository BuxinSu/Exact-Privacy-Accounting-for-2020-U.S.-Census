import numpy as np
import csv, os
import itertools
import pandas as pd
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)
import math
import mpmath as mp
from mpmath import mpf, exp, findroot, almosteq
from collections import Counter
import time
import sys
from fractions import Fraction
import argparse


# ------------------ parse args early (single parser) ------------------
parser = argparse.ArgumentParser(
    description="index of directed path along the US Census Bureau allocation paths"
)
parser.add_argument(
    "--p1",
    type=int,
    required=True,
    help="First path number (path_num_1)",
)
parser.add_argument(
    "--p2",
    type=int,
    required=True,
    help="Second path number (path_num_2)",
)

args = parser.parse_args()

# one record may change from path_num_1 to path_num_2 due to the choice of sensitivity
path_num_1 = int(args.p1)
path_num_2 = int(args.p2)

def format_func(value, tick_number):
    if value < 10:
        return f"{value:.2f}"
    elif value < 100:
        return f"{value:.1f}"
    else:
        return f"{int(value)}"

start_time = time.time()

############# Precision #############
accuracy = 35
mp.dps = accuracy

mp.mp.dps = 35  # Set 50 decimal places of precision

## Test Machine Accuracy

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
        print('Machine Epsilon by our mpmath (not privacy budget):', f"k={k-1}, eps={s}")
        break

print('*'*45)
print('Finish testing accuracy')
print('*'*45)

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

from fractions import Fraction

allocation_run_sorted = sorted(
    allocation_run,
    key=lambda x: Fraction(x.strip("'"))
)

allocation_run = allocation_run_sorted


rho_1 = [Fraction(s) for s in allocation_run_1]
rho_2 = [Fraction(s) for s in allocation_run_2]
if sum(rho_1) != sum(rho_2):
    print("Incorrect Allocations")

print('*'*45)
print('Finish reading allocation path')
print('*'*45)

## Determine Parameters

### Delta

Delta = mp.mpf('1e-25')

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

# optional only for nice expression of U, Find the smallest multiple of 10000 greater than U
step = mp.mpf('10000')
U = step * mp.ceil(U / step)

print(f"U is set to be", U)

### N

summ = mp.mpf('0')
for i in range(len(a_list)):
    summ += a_list[i]**2 * sigma2_list[i]
d = mp.sqrt(- summ / (two * mp.log(Delta / eight)) )
extra_term = d * L * mp.log(two * d * L)

# ==== Final RHS expression ====
N = mp.ceil( two * U + extra_term )

step = mp.mpf('10000')
N = step * mp.ceil(N / step)

print(f"N is set to be", N)

print('*'*45)
print('Finish preparing parameters')
print('*'*45)

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
    return char_func_num_prod(zeta) / char_func_denom_prod

def char_second(zeta):
    return char_func_num_prod(zeta) / char_func_denom_prod

print('*'*45)
print('Finish defining characteristic and weight functions')
print('*'*45)

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
        if abs(char_single) < Delta/(eight * U):
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
df_results.to_csv(f'results_DHC_v2/trade_off_curve/tmp/char_func_half_period_threshold_path_{path_num_1}_to_{path_num_2}.csv', index=False)

# Load previously computed data
df_results = pd.read_csv(f'results_DHC_v2/trade_off_curve/tmp/char_func_half_period_threshold_path_{path_num_1}_to_{path_num_2}.csv', dtype=str)

# Ensure sigma2 values are precisely handled as mpf
df_results['sigma2'] = df_results['sigma2'].apply(mp.mpf)
df_results['zeta'] = df_results['zeta'].apply(mp.mpf)
df_results['a_i_L'] = df_results['a_i_L'].apply(mp.mpf)

midpoints = {}
zeta_interval = {}

for i in range(len(sigma2_unique_list)):
    sigma_i2 = sigma2_unique_list[i]
    a_i_L = a_list_L[index_unique_list[i]]
    half_period = mp.pi / a_i_L

    # Check presence with full precision
    matched_rows = df_results[df_results['a_i_L'] == a_i_L]

    if matched_rows.empty:
        print(f"{sigma_i2} is not used in threshold computation, skipping.")
        continue

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

from bisect import bisect_left

def _nearest_abs_diff(sorted_vals, x):
    """
    Return the absolute difference between x and the nearest value in sorted_vals.
    Assumes sorted_vals is non-empty and sorted ascending (as in your data).
    Works with mpmath.mpf values.
    """
    pos = bisect_left(sorted_vals, x)
    if pos == 0:
        return abs(sorted_vals[0] - x)
    if pos == len(sorted_vals):
        return abs(x - sorted_vals[-1])
    # x lies between neighbors: sorted_vals[pos-1] <= x <= sorted_vals[pos]
    left = sorted_vals[pos - 1]
    right = sorted_vals[pos]
    dl = x - left   # nonnegative distance (mpf)
    dr = right - x  # nonnegative distance (mpf)
    return dl if dl <= dr else dr

# Number of intervals
n_intervals = int(N)

# Output CSV path
csv_path = f"results_DHC_v2/trade_off_curve/threshold/char_fnc_threshold_point_path_{path_num_1}_to_{path_num_2}.csv"

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
            midpoint_list = midpoints[j]      # sorted (as noted)
            zeta_thresh = zeta_interval[j]

            # Distance to the nearest midpoint via binary search (O(log m))
            d = _nearest_abs_diff(midpoint_list, zeta)

            # Check if abs(zeta - closest) < zeta_thresh
            if d >= zeta_thresh:
                passed_all_j = False
                break  # early exit if any j fails

        if passed_all_j:
            writer.writerow([str(zeta)])  # only write valid zeta

print('*'*45)
print('Finish determining nodes')
print('*'*45)

n_intervals = int(N)

# Input and output paths
input_csv = f"results_DHC_v2/trade_off_curve/threshold/char_fnc_threshold_point_path_{path_num_1}_to_{path_num_2}.csv"
output_csv = f"results_DHC_v2/trade_off_curve/char_fnc_eval/char_fnc_eval_path_{path_num_1}_to_{path_num_2}.csv"

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
        # print("zeta:", zeta)
        char = char_first(zeta)
        # print("char_fnc:", char)
        writer.writerow([str(zeta), str(char)])

print('*'*45)
print('Finish evaluating characteristic function')
print('*'*45)

zeta_0_str_list = [str(x) for x in mp.linspace(-mp.floor(U/L), mp.floor(U/L), 901)]

for zeta_0_str in zeta_0_str_list:
    # zeta_0
    zeta_0 = mp.mpf(zeta_0_str)
    print("zeta_0:", zeta_0_str)

    """## T(eps), t(eps)"""
    t_0 = zeta_0 + sum(a_list) / two
    T_0 = - zeta_0 + sum(a_list) / two
    print("t_0, T_0:", t_0, T_0)

    # Input and output paths
    input_csv = f"results_DHC_v2/trade_off_curve/threshold/char_fnc_threshold_point_path_{path_num_1}_to_{path_num_2}.csv"
    output_csv = f"results_DHC_v2/trade_off_curve/weight_eval/trade_off_zeta_0_{zeta_0_str}_alpha_weight_path_{path_num_1}_to_{path_num_2}.csv"
    print("compute alpha weight")

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

    # Input and output paths
    input_csv = f"results_DHC_v2/trade_off_curve/threshold/char_fnc_threshold_point_path_{path_num_1}_to_{path_num_2}.csv"
    output_csv = f"results_DHC_v2/trade_off_curve/weight_eval/trade_off_zeta_0_{zeta_0_str}_beta_weight_path_{path_num_1}_to_{path_num_2}.csv"
    print("compute beta weight")

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

print('*'*45)
print('Finish evaluating weight function')
print('*'*45)

# Load character function and weights
first_char_csv = f"results_DHC_v2/trade_off_curve/char_fnc_eval/char_fnc_eval_path_{path_num_1}_to_{path_num_2}.csv"
with open(first_char_csv, newline='') as f:
    reader = csv.DictReader(f)
    char_vals = [mp.mpf(row['char_fnc']) for row in reader]

# Prepare output storage
output_rows = []

# Loop over zeta_0 values
for zeta_0_str in zeta_0_str_list:
    # zeta_0
    zeta_0 = mp.mpf(zeta_0_str)
    print("zeta_0:", zeta_0_str)
    
    # Read weight CSV for this zeta_0
    first_weight_csv = f"results_DHC_v2/trade_off_curve/weight_eval/trade_off_zeta_0_{zeta_0_str}_alpha_weight_path_{path_num_1}_to_{path_num_2}.csv"
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
    alpha = result_first

    print("alpha_fdp:", alpha)

    # Store result
    output_rows.append([zeta_0_str, mp.nstr(alpha, n=35)])

# Output file path
output_csv = f"results_DHC_v2/trade_off_curve/trade_off_results/trade_off_curve_path_{path_num_1}_to_{path_num_2}.csv"
os.makedirs(os.path.dirname(output_csv), exist_ok=True)

# Save to CSV
with open(output_csv, mode='w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['zeta', 'trade_off_alpha'])
    writer.writerows(output_rows)

print('*'*45)
print('Finish alpha integration and saved')
print('*'*45)

# Inputs
n_intervals = int(N)

# Prepare list to store second term results
beta_results = []

# Loop over zeta_0 values
for zeta_0_str in zeta_0_str_list:
    # zeta_0
    zeta_0 = mp.mpf(zeta_0_str)
    print("zeta_0_str:", zeta_0_str)

    # Load character function (second term)
    second_char_csv = f"results_DHC_v2/trade_off_curve/char_fnc_eval/char_fnc_eval_path_{path_num_1}_to_{path_num_2}.csv"
    with open(second_char_csv, newline='') as f:
        reader = csv.DictReader(f)
        char_vals = [mp.mpf(row['char_fnc']) for row in reader]

    # Load corresponding second weight
    second_weight_csv = f"results_DHC_v2/trade_off_curve/weight_eval/trade_off_zeta_0_{zeta_0_str}_beta_weight_path_{path_num_1}_to_{path_num_2}.csv"
    with open(second_weight_csv, newline='') as f:
        reader = csv.DictReader(f)
        weight_vals = [mp.mpf(row['weight']) for row in reader]

    # Compute integrand
    tobeint_second = [(c * w * one / two) for c, w in zip(char_vals, weight_vals)]

    # Integration
    delta_zeta = two * mp.pi / N
    increment = [delta_zeta * val for val in tobeint_second]
    correction = delta_zeta * tobeint_second[0]
    beta = two * mp.fsum(increment) - correction

    print("beta_fdp:", beta)

    # Store the result
    beta_results.append(mp.nstr(beta, n=35))

print('*'*45)
print('Finish beta integration and saved')
print('*'*45)

# === Read existing CSV with full precision ===
df = pd.read_csv(output_csv, dtype=str)  # read all as string to preserve digits
df['trade_off_beta'] = beta_results

# Ensure both columns are strings (if read from CSV)
df['trade_off_alpha'] = df['trade_off_alpha'].astype(str)
df['trade_off_beta'] = df['trade_off_beta'].astype(str)

# Save CSV back with full-precision values as strings
df.to_csv(output_csv, index=False)
print(f"✅ fdp_curve saved with full {accuracy}-digit precision.")

end_time = time.time()
elapsed_time = (end_time - start_time)/3600
print(f"\n✅ Total runtime: {elapsed_time:.2f} Hours")
