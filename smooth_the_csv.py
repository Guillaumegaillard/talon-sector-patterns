import numpy as np
import pickle
import json
from scipy.interpolate import interp1d

### Parsing the original CSV file, if not already done
try:
    AFvals=pickle.load(open('array_factor/AF_data.dat', 'rb'))
except:
    AF_file="array_factor/array_factor_spherical.csv"

    AFfile = open(AF_file,"r") 
    nline=0
    AFvals={}
    for line in AFfile:
        if nline>0:
            vals=line.split(',')
            if len(vals)!=66:
                print(len(vals))

            if not (vals[0] in AFvals):
                AFvals[vals[0]]={}
            AFvals[vals[0]][vals[1]]=vals[2:]

        nline+=1

    AFfile.close() 

    pickle.dump(AFvals, open('array_factor/AF_data.dat', 'wb'))

### Interpolate it
interp_vals={}
complete_ys=np.arange(-31.5,33.5,2.25) #full range of y angles
cys_str=['-31.5', '-29.25', '-27', '-24.75', '-22.5', '-20.25', '-18', '-15.75', '-13.5', '-11.25', '-9', '-6.75', '-4.5', '-2.25', '0', '2.25', '4.5', '6.75', '9', '11.25', '13.5', '15.75', '18', '20.25', '22.5', '24.75', '27', '29.25', '31.5']

### Complete and interpolate
for azim in AFvals:
    interp_vals[azim]={}
    for col in range(64):
        xs=[-40]
        ys=[0]
        
        ### Complete some missing data
        if azim=='85.5' and col==62:
            xs+=[-27,-24.75]
            ys+=[198,54]
        if azim=='85.5' and col==63:
            xs+=[-27,-24.75]
            ys+=[-12,-189]

        ### Add existing data
        for elev in AFvals[azim]:
            if len(AFvals[azim][elev][col])>1:
                xs.append(float(elev))
                # print(AFvals[azim][elev][col])
                ys.append(float(AFvals[azim][elev][col]))

        ### Complete some missing data
        if azim=='85.5' and col==62:
            xs+=[-9,-6.75,0,11.25,22.5,24.75]
            ys+=[-111,-49,297,51,-72,-356]
        if azim=='85.5' and col==63:
            xs+=[-9,-6.75,0,11.25,22.5,24.75]
            ys+=[-155.5,64.5,-105,17,-59,-115]

        xs.append(40)
        ys.append(0)

        #some prints
        if len(xs)<10:
            print(azim,col,len(xs))

        ### Interpolate the column
        # interp=interp1d(xs,ys, kind='cubic', fill_value="extrapolate")
        # interp=interp1d(xs,ys, kind='quadratic', fill_value="extrapolate")
        interp=interp1d(xs,ys, kind='slinear', fill_value="extrapolate")
        ## Interpolation values for the fixed Y points
        interp_values=interp(complete_ys)
        for cpyid in range(len(complete_ys)):
            if col==0:
                interp_vals[azim][cys_str[cpyid]]=[interp_values[cpyid]]
            else:
                interp_vals[azim][cys_str[cpyid]].append(interp_values[cpyid])

### Save Interpolated CSV
with open('array_factor/AF_spherical_interp.csv','w') as outfile:
    outfile.write("pan,tilt,re00,im00,re01,im01,re02,im02,re03,im03,re04,im04,re05,im05,re06,im06,re07,im07,re08,im08,re09,im09,re10,im10,re11,im11,re12,im12,re13,im13,re14,im14,re15,im15,re16,im16,re17,im17,re18,im18,re19,im19,re20,im20,re21,im21,re22,im22,re23,im23,re24,im24,re25,im25,re26,im26,re27,im27,re28,im28,re29,im29,re30,im30,re31,im31\n")
    for azim in interp_vals:
        rowhead=azim+','
        for elev in interp_vals[azim]:
            row=rowhead+elev+','
            for val in interp_vals[azim][elev]:
                row+=str(val)
                row+=','
            row=row[:-1]+'\n'
            outfile.write(row)



