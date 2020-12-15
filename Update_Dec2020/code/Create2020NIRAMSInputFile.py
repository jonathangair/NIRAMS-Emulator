import sys
sys.path.insert(1,'/Library/Python/2.7/site-packages')
import numpy as np
import matplotlib.pyplot as plt
import corner as triangle
import h5py
import os

os.system("cp nirams_ii_input_2017.h5 nirams_ii_input_2020.h5") 
newfile=h5py.File("nirams_ii_input_2020.h5",'r+')
ng2=newfile.get('five_km_grids')
ng3=ng2.get('meteorological_data')
ng4=ng3.get('daily')


updatefile=h5py.File("nirams-inputs-26aug2020_FromAdam.h5",'r')
g2=updatefile.get('five_km_grids')
g3=g2.get('meteorological_data')
g4=g3.get('daily')

daymax= np.array([31,28,31,30,31,30,31,31,30,31,30,31])
variable='rainfall'
ng5=ng4.get(r'%s' % (variable))
g5=g4.get(r'%s' % (variable))
ng6=ng5.get('rainfall_2001')
dset=ng6.get('rainfall_2001_01_01')
#for key in dset.attrs.keys():
#	print(key)
#	print(dset.attrs[key])

for year in range(2016,2019):
	g6=g5.get(r'%s_%i'%(variable,year))
	ng6=ng5.create_group(r'%s_%i' % (variable,year))
	if year == 2016:
		daymax[1]=29
	else:
		daymax[1]=28
	for month in range(1,13):
		for day in range(1,daymax[month-1]+1):
			thisdata=g6.get(r'%s_%i_%i_%i' % (variable,year,month,day))
			dsetname=(r'%s_%i_%02d_%02d' % (variable,year,month,day))
			ng7 = ng6.create_dataset(dsetname,data=np.array(thisdata))		
			for key in dset.attrs.keys():
				ng7.attrs[key]=dset.attrs[key]
			ng7.attrs['TITLE']=dsetname
var_list = ['min_temp', 'max_temp', 'mean_temp']
for variable in var_list:
	ng5=ng4.get(r'%s' % (variable))
	ng6=ng5.get(r'%s_2011' % (variable))
	dset=ng6.get(r'%s_2011_01_01' % (variable))

	for year in range(2016,2019):
		ng6=ng5.create_group(r'%s_%i' % (variable,year))
		if year == 2016:
			daymax[1]=29
		else:
			daymax[1]=28
		for month in range(1,13):
			for day in range(1,daymax[month-1]+1):
				thisdata=np.array(g3.get(r'daily_%s_%s_%i_%s_%i_%i_%i' % (variable,variable,year,variable,year,month,day)))
				dsetname=(r'%s_%i_%02d_%02d' % (variable,year,month,day))
				ng7=ng6.create_dataset(dsetname,data=thisdata)		
				for key in dset.attrs.keys():
					ng7.attrs[key]=dset.attrs[key]
				ng7.attrs['TITLE']=dsetname

ng4=ng3.get('monthly')
ng5=ng4.get('est_pm_pet')
dset=ng5.get('est_pm_pet_2001_01')
for year in range(2016,2019):
	for month in range(1,13):
		thisdata=np.array(g3.get(r'monthly_Thorn_PET_Thorn_PET_%i_%i' % (year,month)))
		ng7=ng5.create_dataset((r'est_pm_pet_%i_%02d' % (year,month)),data=thisdata)
		for key in dset.attrs.keys():
			ng7.attrs[key]=dset.attrs[key]
		ng7.attrs['TITLE']=r'est_pm_pet_%i_%02d' % (year,month)
ng5=ng4.get('pm_pet')
dset=ng5.get('pm_pet_2001_01')
for year in range(2015,2019):
	for month in range(1,13):
		thisdata=np.array(g3.get(r'monthly_PM_PET_PM_PET_%i_%i' % (year,month)))
		ng7=ng5.create_dataset((r'pm_pet_%i_%02d' % (year,month)),data=thisdata)
		for key in dset.attrs.keys():
			ng7.attrs[key]=dset.attrs[key]
		ng7.attrs['TITLE']=r'pm_pet_%i_%02d' % (year,month)
updatefile.close()

# For the nitrate data we will copy old data for the new years 2016 to 2018. We take year 2015 as the reference and 
# copy this to each year in the range 2016 - 2018.
ng2=newfile.get('one_km_grids')
ng3=ng2.get('iacs_pet_facts')
ng4=ng2.get('inorganic_n')
ng5=ng2.get('n_deposition')
ng6=ng2.get('n_uptake')
ng7=ng2.get('organic_n')
pet_data=ng3.get('pet_fact_2015')
in_gr_data=ng4.get('in_gr_2015')
in_ot_data=ng4.get('in_ot_2015')
in_sp_data=ng4.get('in_sp_2015')
in_wi_data=ng4.get('in_wi_2015')
n_dep_data=ng5.get('n_dep_2015')
up_gr_data=ng6.get('up_gr_2015')
up_ot_data=ng6.get('up_ot_2015')
up_sp_data=ng6.get('up_sp_2015')
up_wi_data=ng6.get('up_wi_2015')
or_gr_data=ng7.get('or_gr_2015')
or_ot_data=ng7.get('or_ot_2015')
or_sp_data=ng7.get('or_sp_2015')
or_wi_data=ng7.get('or_wi_2015')
for year in range(2016,2019):
	ndset=ng3.create_dataset(r'pet_fact_%i' % (year),data=np.array(pet_data))
	for key in pet_data.attrs.keys():
		ndset.attrs[key]=pet_data.attrs[key]
	ndset.attrs['TITLE']=r'pet_fact_%i' % (year)
	ndset=ng4.create_dataset(r'in_gr_%i' % (year),data=np.array(in_gr_data))
	for key in in_gr_data.attrs.keys():
		ndset.attrs[key]=in_gr_data.attrs[key]
	ndset.attrs['TITLE']=r'in_gr_%i' % (year)
	ndset=ng4.create_dataset(r'in_ot_%i' % (year),data=np.array(in_ot_data))
	for key in in_ot_data.attrs.keys():
		ndset.attrs[key]=in_ot_data.attrs[key]
	ndset.attrs['TITLE']=r'in_ot_%i' % (year)
	ndset=ng4.create_dataset(r'in_sp_%i' % (year),data=np.array(in_sp_data))
	for key in in_sp_data.attrs.keys():
		ndset.attrs[key]=in_sp_data.attrs[key]
	ndset.attrs['TITLE']=r'in_sp_%i' % (year)
	ndset=ng4.create_dataset(r'in_wi_%i' % (year),data=np.array(in_wi_data))
	for key in in_wi_data.attrs.keys():
		ndset.attrs[key]=in_wi_data.attrs[key]
	ndset.attrs['TITLE']=r'in_wi_%i' % (year)
	ndset=ng5.create_dataset(r'n_dep_%i' % (year),data=np.array(n_dep_data))
	for key in n_dep_data.attrs.keys():
		ndset.attrs[key]=n_dep_data.attrs[key]
	ndset.attrs['TITLE']=r'n_dep_%i' % (year)
	ndset=ng6.create_dataset(r'up_gr_%i' % (year),data=np.array(up_gr_data))
	for key in up_gr_data.attrs.keys():
		ndset.attrs[key]=up_gr_data.attrs[key]
	ndset.attrs['TITLE']=r'up_gr_%i' % (year)
	ndset=ng6.create_dataset(r'up_ot_%i' % (year),data=np.array(up_ot_data))
	for key in up_ot_data.attrs.keys():
		ndset.attrs[key]=up_ot_data.attrs[key]
	ndset.attrs['TITLE']=r'up_ot_%i' % (year)
	ndset=ng6.create_dataset(r'up_sp_%i' % (year),data=np.array(up_sp_data))
	for key in up_sp_data.attrs.keys():
		ndset.attrs[key]=up_sp_data.attrs[key]
	ndset.attrs['TITLE']=r'up_sp_%i' % (year)
	ndset=ng6.create_dataset(r'up_wi_%i' % (year),data=np.array(up_wi_data))
	for key in up_wi_data.attrs.keys():
		ndset.attrs[key]=up_wi_data.attrs[key]
	ndset.attrs['TITLE']=r'up_wi_%i' % (year)
	ndset=ng7.create_dataset(r'or_gr_%i' % (year),data=np.array(or_gr_data))
	for key in or_gr_data.attrs.keys():
		ndset.attrs[key]=or_gr_data.attrs[key]
	ndset.attrs['TITLE']=r'or_gr_%i' % (year)
	ndset=ng7.create_dataset(r'or_ot_%i' % (year),data=np.array(or_ot_data))
	for key in or_ot_data.attrs.keys():
		ndset.attrs[key]=or_ot_data.attrs[key]
	ndset.attrs['TITLE']=r'or_ot_%i' % (year)
	ndset=ng7.create_dataset(r'or_sp_%i' % (year),data=np.array(or_sp_data))
	for key in or_sp_data.attrs.keys():
		ndset.attrs[key]=or_sp_data.attrs[key]
	ndset.attrs['TITLE']=r'or_sp_%i' % (year)
	ndset=ng7.create_dataset(r'or_wi_%i' % (year),data=np.array(or_wi_data))
	for key in or_wi_data.attrs.keys():
		ndset.attrs[key]=or_wi_data.attrs[key]
	ndset.attrs['TITLE']=r'or_wi_%i' % (year)

newfile.close()

