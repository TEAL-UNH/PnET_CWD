'''Jack's recreation of Andy's 'CWD script try2.py'''
import pathlib
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import time
#import math
st = time.perf_counter()
plt.style.use('seaborn-whitegrid')

'''No hardcoded pathways. '''
path = pathlib.Path.cwd()
input_path = path.parent / 'inputs'

'''Read in inputs'''
dynam_in = pd.read_csv(input_path / 'CWD_dynamic_input_params.csv')
disturb_in = pd.read_csv(input_path / 'CWD_disturbance_inputs.csv')

'''Dynamic inputs to dictionary of the same name '''
dynam_in = dynam_in.set_index('Parameter')['Value'].to_dict()


'''Build dataframe with initial conditions...'''
# initialize empty pandas dataframe
age = pd.DataFrame()

# populate with Year, timestep from dynam_in
age['year'] = [
    *range(int(dynam_in['Start_Year']), int(dynam_in['End_Year'])+1)]

age['timestep'] = [
    *range(0, int(dynam_in['End_Year']) - int(dynam_in['Start_Year'])+1)]

# Add stand age col.
age['stand_age'] = age['year'] - dynam_in['Start_Year']

# Add disturbances from disturb_in
age['disturb_mort'] = age['year'].map(
    disturb_in.set_index('Year')['Disturbance_Mortality']).fillna(0)
age['disturb_bio_remove'] = age['year'].map(
    disturb_in.set_index('Year')['Disturbance_Biomass_Removal']).fillna(0)


'''define dynamic functions'''


def lai_scalar_eqn(stand_age, yrs_to_full_can):
    ''' Dynamic LAI Scalar based on stand age'''

    if stand_age >= yrs_to_full_can:
        lai_scalar = 1
    else:
        lai_scalar = (stand_age + 1) / yrs_to_full_can

    return lai_scalar


def fol_n_eqn(stand_age, fol_n_csa, fol_n_cea, fol_n_init, fol_n_end):
    '''Dynamic Foliar N based on stand age, etc... '''

    if stand_age < fol_n_csa:
        fol_n = fol_n_init
    if stand_age >= fol_n_csa and stand_age < fol_n_cea:
        fol_n = fol_n_init - ((stand_age - fol_n_csa) /
                              (fol_n_cea - fol_n_csa)) * (fol_n_init - fol_n_end)
    if stand_age >= fol_n_cea:
        fol_n = fol_n_end

    return fol_n


def npp_decline_scalar_eqn(stand_age, age_start_npp_decline, mat):
    '''dynamic npp decline scalar based on stand age...'''

    if stand_age < age_start_npp_decline:
        npp_decline_scalar = 1
    if stand_age >= age_start_npp_decline:
        npp_decline_scalar = ((156 * pow(stand_age, -0.226) * pow(mat, 1.003)) /
                              (156 * pow(age_start_npp_decline, -0.226) * pow(mat, 1.003)))

    return npp_decline_scalar


'''Use dynamic functions to populate columns'''
# LAI scalar
age['lai_scalar'] = age.apply(lambda row: lai_scalar_eqn(
    row['stand_age'],
    dynam_in['years_until_full_canopy']
), axis=1)


# foliar nitrogen scalar
age['fol_n'] = age.apply(lambda row: fol_n_eqn(
    row['stand_age'],
    dynam_in['folN_change_start_age'],
    dynam_in['folN_change_end_age'],
    dynam_in['folN_init'],
    dynam_in['folN_end']
), axis=1)

# NPP decline scalar
age['npp_decline_scalar'] = age.apply(lambda row: npp_decline_scalar_eqn(
    row['stand_age'],
    dynam_in['age_start_NPP_decline'],
    dynam_in['MAT']
), axis=1)

# Adds wood mortality rate, snag fall rate, snag decay rate, and downed woody debris decay rate
age["wood_mort_rate"] = (dynam_in['wood_mort_rate_a'] + age["stand_age"])/(
    dynam_in['wood_mort_rate_a'] * (dynam_in['wood_mort_rate_b'] + age["stand_age"]))

age["snag_fall_rate"] = (dynam_in['snag_fall_rate_a'] + age["stand_age"])/(
    dynam_in['snag_fall_rate_a'] * (dynam_in['snag_fall_rate_b'] + age["stand_age"]))

age["snag_decay_rate"] = (dynam_in['snag_decay_rate_a'] + age["stand_age"])/(
    dynam_in['snag_decay_rate_a'] * (dynam_in['snag_decay_rate_b'] + age["stand_age"]))

age["dwd_decay_rate"] = (dynam_in['dwd_decay_rate_a'] + age["stand_age"])/(
    dynam_in['dwd_decay_rate_a'] * (dynam_in['dwd_decay_rate_b'] + age["stand_age"]))

# Adds woody production column
age['wood_production'] = ((age['fol_n'] * 200 + 50) * 
    age['lai_scalar'] * age['npp_decline_scalar'])




#                         ///// Matrix stuff \\\\\

# Define number of years between start and stop -- used to build empty matrices
d_yr = int(dynam_in['End_Year']) - int(dynam_in['Start_Year']) + 1
yr = int(dynam_in['Start_Year'])


# Initialize empty lists to be filled in loop
wood_biomass_ls = []
standing_necromass_ls = []
downed_necromass_ls = []
standing_necromass_N_ls = []
standing_necromass_CN_ls = []
downed_necromass_N_ls = []
downed_necromass_CN_ls = []

# Initialize empty Matrices to be filled in loop
# // Initalize with NaN instead? 
snag_cohorts = np.zeros((d_yr,d_yr))
dwd_cohorts = np.zeros((d_yr,d_yr))
snag_cohorts_N = np.zeros((d_yr,d_yr))
dwd_cohorts_N = np.zeros((d_yr,d_yr))
snag_cohorts_CN = np.zeros((d_yr,d_yr))
dwd_cohorts_CN = np.zeros((d_yr,d_yr))

## The loop...where i is row (year) and j is col. 
for i in range(d_yr):
    #print('istep: ',i)
    #print('Year: ', yr)
    j = 0 #need to reset column each loop

# wood biomass has to come before 
    if i == 0:
        wood_biomass = age.iloc[i]['wood_production'] 
    else:
        wood_biomass = wood_biomass * (1 - age.iloc[i-1]['wood_mort_rate']) + age.iloc[i]['wood_production']

    wood_biomass_ls.append(wood_biomass)
    #print('wood_biomass: ', wood_biomass)

####     snag necromass   ####
    while j <= i:
        if i == j: # snag_cohorts[i-1, j] == 0: # I think this way it'll never fail
            snag_cohorts[i,j] = ((age.iloc[i]['wood_mort_rate'] * wood_biomass_ls[i]) * (1 - age.iloc[j]['snag_fall_rate'])) 
        else: 
            snag_cohorts[i,j] = snag_cohorts[i-1,j] - ((age.iloc[j]['snag_fall_rate'] * snag_cohorts[i-1,j]) + (snag_cohorts[i-1,j] * age.iloc[j]['snag_decay_rate']))
    
    
####   downed necromass     ####
        if i == j:
            dwd_cohorts[i,j] = ((age.iloc[i]['wood_mort_rate'] * wood_biomass_ls[i]) - snag_cohorts[i,j])
        else:
            dwd_cohorts[i,j] = (dwd_cohorts[i-1,j] + (snag_cohorts[i,j] * 
                age.iloc[j]['snag_fall_rate']) - (dwd_cohorts[i-1,j] * age.iloc[j]['dwd_decay_rate']))


####     snag Nitrogen ####
        if i == j: # snag_cohorts[i-1, j] == 0: # I think this way it'll never fail
            snag_cohorts_N[i,j] = snag_cohorts[i,j] * dynam_in['wood_perc_N'] 
        else:
            if snag_cohorts_CN[i-1,j] > dynam_in['critical_CN']:
                snag_cohorts_N[i,j] = (snag_cohorts_N[i-1,j] * (1 - age.iloc[j]['snag_fall_rate'])) 
            else: #snag_cohorts_CN[i,j] <= dynam_in['critical_CN']:
                snag_cohorts_N[i,j] = snag_cohorts_N[i-1,j] - ((snag_cohorts[i-1,j]
                - snag_cohorts[i,j]) * (snag_cohorts_N[i-1,j]/snag_cohorts[i-1,j]))
                
              
####    snag C/N ####
        snag_cohorts_CN[i,j] = (snag_cohorts[i,j] * 0.5) / snag_cohorts_N[i,j]


####     dwd Nitrogen ####
        if i == j: # snag_cohorts[i-1, j] == 0: # I think this way it'll never fail
            dwd_cohorts_N[i,j] = dwd_cohorts[i,j] * dynam_in['wood_perc_N'] 
        else:
            if dwd_cohorts_CN[i-1,j] > dynam_in['critical_CN']:
                dwd_cohorts_N[i,j] = dwd_cohorts_N[i-1,j] + snag_cohorts_N[i-1,j] * age.iloc[j]['snag_fall_rate'] 
            else: #snag_cohorts_CN[i,j] <= dynam_in['critical_CN']:
                dwd_cohorts_N[i,j] = (dwd_cohorts_N[i-1,j] * (1 - age.iloc[j]['dwd_decay_rate'])) 
                + (snag_cohorts_N[i-1,j] * age.iloc[j]['snag_fall_rate'])
                
              
####    dwd C/N ####
        dwd_cohorts_CN[i,j] = (dwd_cohorts[i,j] * 0.5) / dwd_cohorts_N[i,j]

# ####     dwd nitrogen  ####
#         if i == j: # snag_cohorts[i-1, j] == 0: # I think this way it'll never fail
#             dwd_cohorts_N[i,j] = dynam_in['wood_N_decay_a']/1000000 * dwd_cohorts[i,j]
#         else: 
#             dwd_cohorts_N[i,j] = (dynam_in['wood_N_decay_a'] * pow(i+1,dynam_in['wood_N_decay_b']) * math.exp(i*dynam_in['wood_N_decay_k']))/1000000*dwd_cohorts[i,j]


# ####     dwd CN ####
#         if i == j:
#             dwd_cohorts_CN[i,j] = 0.5/(dynam_in['wood_N_decay_a']/1000000)
#         else: 
#             dwd_cohorts_CN[i,j] = (dwd_cohorts[i,j] * 0.5) / dwd_cohorts_N[i,j] 


        #while loop(s) column iterator
        j += 1
        ##
    #necromass list append (condensed)
    standing_necromass_ls.append(snag_cohorts[i].sum()) 
    downed_necromass_ls.append(dwd_cohorts[i].sum())
    standing_necromass_N_ls.append(snag_cohorts_N[i].sum())
    standing_necromass_CN_ls.append(snag_cohorts_CN[i].mean())
    downed_necromass_N_ls.append(dwd_cohorts_N[i].sum())
    downed_necromass_CN_ls.append(dwd_cohorts_CN[i].mean())


#lists to dataframe
age['wood_biomass'] = wood_biomass_ls
age['standing_necromass'] = standing_necromass_ls
age['downed_necromass'] = downed_necromass_ls
age['standing_necromass_N'] = standing_necromass_N_ls
age['standing_necromass_CN'] = standing_necromass_CN_ls
age['downed_necromass_N'] = downed_necromass_N_ls
age['downed_necromass_CN'] = downed_necromass_CN_ls


#total necromass column
age['total_necromass'] = age['standing_necromass'] + age['downed_necromass']
age['total_necromass_N'] = age['standing_necromass_N'] + age['downed_necromass_N']

print('Time to run:', time.perf_counter() - st)



''' Plot it up '''
# plt.scatter(x=age['stand_age'], y=age['wood_biomass'])
# plt.scatter(x=age['stand_age'], y=age['downed_necromass'])
# plt.scatter(x=age['stand_age'], y=age['total_necromass'])
# plt.scatter(x=age['stand_age'], y=age['standing_necromass'])
# plt.scatter(x=age['stand_age'], y=age['standing_necromass_CN'])
# plt.scatter(x=age['stand_age'], y=age['wood_mort_rate'])
plt.scatter(x=age['stand_age'], y=age['standing_necromass_N'])
plt.scatter(x=age['stand_age'], y=age['downed_necromass_N'])
plt.scatter(x=age['stand_age'], y=age['total_necromass_N'])

