'''                 Coarse Wood Debris Module for PnET               '''
import pathlib
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import time

# start time 
st = time.perf_counter()



'''No hardcoded pathways. '''
path = pathlib.Path.cwd()
in_path = path.parent / 'inputs' 
out_path = path.parent / 'outputs'

'''Read in inputs'''
dynam_in = pd.read_csv(in_path / 'CWD_dynamic_input_params.csv')
disturb_in = pd.read_csv(in_path / 'CWD_disturbance_inputs.csv')


'''For now specify the model version in script
Set model_version = 
'default'
'no snags'
'no bark'
'no snags or bark'
'''
model_version = 'no snags'

if model_version == 'no snags':
    dynam_in = dynam_in.set_index('Parameter')['No Snags'].to_dict()
    output = out_path / 'no_snags_age_output'
elif model_version == 'no bark':
    dynam_in = dynam_in.set_index('Parameter')['No Bark'].to_dict()
    output = out_path / 'no_bark_age_output'
elif model_version == 'no snags or bark':
    dynam_in = dynam_in.set_index('Parameter')['No Snags or Bark'].to_dict()
    output = out_path / 'no_snags_or_bark_age_output'
else:
    dynam_in = dynam_in.set_index('Parameter')['Default'].to_dict()
    output = out_path / 'default_age_output'



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




#                         ///// Matrix \\\\\

# Define number of years between start and stop -- used to build empty matrices
d_yr = int(dynam_in['End_Year']) - int(dynam_in['Start_Year']) + 1



''' Initialize empty lists to be filled in loop 
For wood unless indicated as bark
bio = living, standing (biomass)
snag = dead, standing
dwd = dead, downed
'''

# biomass/necromass 
bio_ls = []
bark_bio_ls = []
snag_ls = []
dwd_ls = []
bark_snag_ls = []
bark_dwd_ls = []
# nitrogen
snag_n_ls = []
dwd_n_ls = []
bark_snag_n_ls = []
bark_dwd_n_ls = []
# C/N ratio
snag_cn_ls = []
dwd_cn_ls = []
bark_snag_cn_ls = []
bark_dwd_cn_ls = []

'''Initialize empty cohort matrices to be filled in loop'''
# biomass/necromass
snag_coho = np.zeros((d_yr,d_yr))
dwd_coho = np.zeros((d_yr,d_yr))
bark_snag_coho = np.zeros((d_yr,d_yr))
bark_dwd_coho = np.zeros((d_yr,d_yr))
# nitrogen
snag_n_coho = np.zeros((d_yr,d_yr))
dwd_n_coho = np.zeros((d_yr,d_yr))
bark_snag_n_coho = np.zeros((d_yr,d_yr))
bark_dwd_n_coho = np.zeros((d_yr,d_yr))
# C/N ratio
snag_cn_coho = np.zeros((d_yr,d_yr))
dwd_cn_coho = np.zeros((d_yr,d_yr))
bark_snag_cn_coho = np.zeros((d_yr,d_yr))
bark_dwd_cn_coho = np.zeros((d_yr,d_yr))


''' The LOOP '''
# where i is row (year), and j is column
print('Model Options: ')
if dynam_in['snag_fall_rate_a'] == 1 and dynam_in['snag_fall_rate_b'] == 1 :
    print('- Not Including Snags')
else:
    print('- Including Snags')


if dynam_in['bark_frac_wood'] == 0:
    print('- Bark and Wood Together')
else: 
    print('- Bark and Wood Separate \n')

for i in range(d_yr):

    j = 0 # reset column to zero each loop


    # ~~~ biomass ~~~ 
    if i == 0:
        wood_bio = age.iloc[i]['wood_production'] * (1 - dynam_in['bark_frac_wood'])
        bark_bio = age.iloc[i]['wood_production'] * (dynam_in['bark_frac_wood'])
    else:
        wood_bio = (wood_bio * (1 - age.iloc[i-1]['wood_mort_rate']) + 
            age.iloc[i]['wood_production'] * (1 - dynam_in['bark_frac_wood']))
        bark_bio = (bark_bio * (1 - age.iloc[i-1]['wood_mort_rate']) + 
            age.iloc[i]['wood_production'] * (dynam_in['bark_frac_wood']))

    bio_ls.append(wood_bio)
    bark_bio_ls.append(bark_bio)


    # ~~~ snag necromass ~~~
    while j <= i:
        if i == j: 
            snag_coho[i,j] = ((age.iloc[i]['wood_mort_rate'] * bio_ls[i]) * 
                (1 - age.iloc[j]['snag_fall_rate'])) 
            bark_snag_coho[i,j] = ((age.iloc[i]['wood_mort_rate'] * bark_bio_ls[i]) * 
                (1 - age.iloc[j]['snag_fall_rate'])) 
                
        else: 
            snag_coho[i,j] = (snag_coho[i-1,j] - ((age.iloc[j]['snag_fall_rate'] * 
                snag_coho[i-1,j]) + (snag_coho[i-1,j] * age.iloc[j]['snag_decay_rate'])))
            bark_snag_coho[i,j] = (bark_snag_coho[i-1,j] - ((age.iloc[j]['snag_fall_rate'] * 
                bark_snag_coho[i-1,j]) + (bark_snag_coho[i-1,j] * age.iloc[j]['snag_decay_rate'])))
    
    
    # ~~~ downed necromass ~~~  
        if i == j:
            dwd_coho[i,j] = ((age.iloc[i]['wood_mort_rate'] * bio_ls[i]) - snag_coho[i,j])
            bark_dwd_coho[i,j] = ((age.iloc[i]['wood_mort_rate'] * bark_bio_ls[i]) - bark_snag_coho[i,j])
        else:
            dwd_coho[i,j] = (dwd_coho[i-1,j] + (snag_coho[i,j] * 
                age.iloc[j]['snag_fall_rate']) - (dwd_coho[i-1,j] * age.iloc[j]['dwd_decay_rate']))
            bark_dwd_coho[i,j] = (bark_dwd_coho[i-1,j] + (bark_snag_coho[i,j] * 
                age.iloc[j]['snag_fall_rate']) - (bark_dwd_coho[i-1,j] * age.iloc[j]['dwd_decay_rate']))


    # ~~~ snag nitrogen ~~~
        if i == j: 
            snag_n_coho[i,j] = snag_coho[i,j] * dynam_in['wood_perc_N'] 
            bark_snag_n_coho[i,j] = bark_snag_coho[i,j] * dynam_in['bark_perc_N'] 
        else:
            if snag_cn_coho[i-1,j] > dynam_in['critical_CN']:
                snag_n_coho[i,j] = snag_n_coho[i-1,j] * (1 - age.iloc[j]['snag_fall_rate']) 
            else: 
                snag_n_coho[i,j] = (snag_n_coho[i-1,j] - ((snag_coho[i-1,j] - 
                    snag_coho[i,j]) * (snag_n_coho[i-1,j] / snag_coho[i-1,j])))

            if bark_snag_cn_coho[i-1,j] > dynam_in['critical_CN']:
                bark_snag_n_coho[i,j] = bark_snag_n_coho[i-1,j] * (1 - age.iloc[j]['snag_fall_rate']) 
            else: 
                bark_snag_n_coho[i,j] = (bark_snag_n_coho[i-1,j] - ((bark_snag_coho[i-1,j] - 
                    bark_snag_coho[i,j]) * (bark_snag_n_coho[i-1,j] / bark_snag_coho[i-1,j])))
                
              
    # ~~~ snag c/n ~~~ 
        if snag_n_coho[i,j] == 0:
            snag_cn_coho[i,j] = 9999 # to prevent div.by.zero errors. 
        else:
            snag_cn_coho[i,j] = (snag_coho[i,j] * 0.5) / snag_n_coho[i,j]

        if bark_snag_n_coho[i,j] == 0:
            bark_snag_cn_coho[i,j] = 9999 # to prevent div.by.zero errors. 
        else:
            bark_snag_cn_coho[i,j] = (bark_snag_coho[i,j] * 0.5) / bark_snag_n_coho[i,j]

    # ~~~ dwd nitrogen ~~~
        if i == j: 
            dwd_n_coho[i,j] = dwd_coho[i,j] * dynam_in['wood_perc_N'] 
            bark_dwd_n_coho[i,j] = bark_dwd_coho[i,j] * dynam_in['bark_perc_N']
        else:
            if dwd_cn_coho[i-1,j] > dynam_in['critical_CN']:
                dwd_n_coho[i,j] = dwd_n_coho[i-1,j] + snag_n_coho[i-1,j] * age.iloc[j]['snag_fall_rate'] 
            else: 
                dwd_n_coho[i,j] = (dwd_n_coho[i-1,j] - ((dwd_coho[i-1,j] * 0.5 * age.iloc[j]['dwd_decay_rate']) /
                    dwd_cn_coho[i-1,j]) + (snag_n_coho[i-1,j] * age.iloc[j]['snag_fall_rate']))

            if bark_dwd_cn_coho[i-1,j] > dynam_in['critical_CN']:
                bark_dwd_n_coho[i,j] = bark_dwd_n_coho[i-1,j] + bark_snag_n_coho[i-1,j] * age.iloc[j]['snag_fall_rate'] 
            else: 
                bark_dwd_n_coho[i,j] = (bark_dwd_n_coho[i-1,j] - ((bark_dwd_coho[i-1,j] * 0.5 * age.iloc[j]['dwd_decay_rate']) /
                    bark_dwd_cn_coho[i-1,j]) + (bark_snag_n_coho[i-1,j] * age.iloc[j]['snag_fall_rate']))
              
    # ~~~ dwd c/n ~~~
        dwd_cn_coho[i,j] = (dwd_coho[i,j] * 0.5) / dwd_n_coho[i,j]
        bark_dwd_cn_coho[i,j] = (bark_dwd_coho[i,j] * 0.5) / bark_dwd_n_coho[i,j]

        
        j += 1 # while loop column iterator
       
    # necromass list append
    snag_ls.append(snag_coho[i].sum()) 
    dwd_ls.append(dwd_coho[i].sum())
    snag_n_ls.append(snag_n_coho[i].sum())
    snag_cn_ls.append(snag_cn_coho[i].mean())
    dwd_n_ls.append(dwd_n_coho[i].sum())
    dwd_cn_ls.append(dwd_cn_coho[i].mean())

    bark_snag_ls.append(bark_snag_coho[i].sum()) 
    bark_dwd_ls.append(bark_dwd_coho[i].sum())
    bark_snag_n_ls.append(bark_snag_n_coho[i].sum())
    bark_snag_cn_ls.append(bark_snag_cn_coho[i].mean())
    bark_dwd_n_ls.append(bark_dwd_n_coho[i].sum())
    bark_dwd_cn_ls.append(bark_dwd_cn_coho[i].mean())




#lists to dataframe
age['wood_biomass'] = bio_ls
age['standing_necromass'] = snag_ls
age['downed_necromass'] = dwd_ls
age['standing_necromass_N'] = snag_n_ls
age['downed_necromass_N'] = dwd_n_ls
age['standing_necromass_CN'] = snag_cn_ls
age['downed_necromass_CN'] = dwd_cn_ls

age['bark_biomass'] = bark_bio_ls
age['standing_necromass_bark'] = bark_snag_ls
age['downed_necromass_bark'] = bark_dwd_ls
age['standing_necromass_bark_N'] = bark_snag_n_ls
age['downed_necromass_bark_N'] = bark_dwd_n_ls
age['standing_necromass_bark_CN'] = bark_snag_cn_ls
age['downed_necromass_bark_CN'] = bark_dwd_cn_ls


#total necromass column
age['wood_necromass'] = age['standing_necromass'] + age['downed_necromass']
age['wood_necromass_N'] = age['standing_necromass_N'] + age['downed_necromass_N']

age['bark_necromass'] = age['standing_necromass_bark'] + age['downed_necromass_bark']
age['bark_necromass_N'] = age['standing_necromass_bark_N'] + age['downed_necromass_bark_N']

age['total_necromass'] = age['wood_necromass'] + age['bark_necromass']
age['total_necromass_N'] = age['wood_necromass_N'] + age['bark_necromass_N']

age.to_csv(output)


''' Plot it up '''
plt.style.use('seaborn-whitegrid') #plot style

# plt.scatter(x=age['stand_age'], y=age['wood_biomass'])
# plt.scatter(x=age['stand_age'], y=age['downed_necromass'])
# plt.scatter(x=age['stand_age'], y=age['total_necromass'])
# plt.scatter(x=age['stand_age'], y=age['standing_necromass'])
# plt.scatter(x=age['stand_age'], y=age['standing_necromass_CN'])
# plt.scatter(x=age['stand_age'], y=age['wood_mort_rate'])

plt.scatter(x=age['stand_age'], y=age['wood_necromass_N'])
plt.scatter(x=age['stand_age'], y=age['bark_necromass_N'])
plt.scatter(x=age['stand_age'], y=age['snag_fall_rate'])




print('Time to run:', time.perf_counter() - st)