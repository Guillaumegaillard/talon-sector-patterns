import numpy as np
import pickle
from multiprocessing import Pool,cpu_count
import json
import pathlib
current_path=str(pathlib.Path(__file__).parent.resolve())+"/"

### Express Talons Codebooks:
# in a similar format as in https://github.com/wigig-tools/wigig-module
# as json file
# as returned dict


AZIMUTH_CARDINALITY    = 361
ELEVATION_CARDINALITY  = 181
NB_ELEMENTS_ARRAY      = 32
antennas={}
global_sed=.16#.15#0.01#
reload_AF_File=False
allowed_multi_pool=False

AF_file=current_path+"array_factor/AF_spherical_interp.csv"

### A rapid print
def show_AF():
    AFfile = open(AF_file,"r") 
    nline=0
    for line in AFfile:
        # if nline>0:
            if line.split(',')[0][:4]=='-157':
                print(line[:12], nline)
                # break
                # print((line.split(',')))
                nline+=1
            # else:
            #     print(line.split(',')[0],"wwwwwwwwwwwww")
            #     break
    print(nline)


### load Array Factor. 
# Measurements are packman like (backside of the Talons not measured) 
# => Set complete_packmen to extrapolate a 360-degree-AF
def update_SV(complete_packmen=False):
    global steeringVector_as_array, steeringVector

    if reload_AF_File:

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

        steeringVector={}
        steeringVector_as_array=np.zeros((AZIMUTH_CARDINALITY,ELEVATION_CARDINALITY,NB_ELEMENTS_ARRAY,2))
        for y in range(AZIMUTH_CARDINALITY):
            steeringVector[y]={}
            shifty=y-360 if y>179 else y
            found_pan=False
            pan_key=''
            for k in AFvals.keys():
                if shifty >= float(k)-1.125 and shifty < float(k)+1.125 :
                    found_pan=True
                    # print(y,k)
                    pan_key=k
                    if y>=358 or y <=2:
                        print(y,k)
                    break

            for z in range(ELEVATION_CARDINALITY):
                steeringVector[y][z]={}

                found_tilt=False
                if found_pan:
                    AFtiltkeys=[(k,float(k)+ 90) for k in AFvals[pan_key].keys()]

                    tilt_key=''
                    for k in AFtiltkeys:
                        if z >= k[1]-1.125 and z < k[1]+1.125 :
                            found_tilt=True
                            tilt_key=k[0]
                            break
                
                for x in range(NB_ELEMENTS_ARRAY):
                    if found_pan and found_tilt:
                        steeringVector[y][z][x]={'re':AFvals[pan_key][tilt_key][2*x].strip(),'im':AFvals[pan_key][tilt_key][2*x+1].strip() , 'found':True}
                    else:
                        steeringVector[y][z][x]={'re':'0','im':'0', 'found':False}

        if complete_packmen:
            elev_smoother=10
            for x in range(NB_ELEMENTS_ARRAY):
                for y in range(AZIMUTH_CARDINALITY):
                    boundaries=[0,0] # first y last y
                    for z in range(ELEVATION_CARDINALITY):
                        if (not(steeringVector[y][z][x]['found'])): 
                            if (boundaries[0]==z):
                                boundaries[0]=z+1
                        else:
                            if (boundaries[1]!=z):
                                boundaries[1]=z
                    if boundaries[0]<boundaries[1]:
                        # print(x,y,boundaries, steeringVector[y][boundaries[0]][x])
                        for j in range(boundaries[0]-elev_smoother,boundaries[0]):
                            steeringVector[y][j][x]['re']=str(
                                (j-boundaries[0]+elev_smoother)/elev_smoother*(
                                    float(steeringVector[y][boundaries[0]][x]['re']) \
                                        if (len(steeringVector[y][boundaries[0]][x]['re'])>0) else 0
                                )
                            )
                            steeringVector[y][j][x]['im']=str(
                                (j-boundaries[0]+elev_smoother)/elev_smoother*(
                                    float(steeringVector[y][boundaries[0]][x]['im']) \
                                        if (len(steeringVector[y][boundaries[0]][x]['im'])>0) else 0
                                )
                            )
                            steeringVector[y][j][x]['found']=True
                        for j in range(boundaries[1]+1,boundaries[1]+elev_smoother+1):
                            steeringVector[y][j][x]['re']=str(
                                (elev_smoother-j+boundaries[1])/elev_smoother*(
                                    float(steeringVector[y][boundaries[1]][x]['re']) \
                                        if (len(steeringVector[y][boundaries[1]][x]['re'])>0) else 0
                                )
                            )
                            steeringVector[y][j][x]['im']=str(
                                (elev_smoother-j+boundaries[1])/elev_smoother*(
                                    float(steeringVector[y][boundaries[1]][x]['im']) \
                                        if (len(steeringVector[y][boundaries[1]][x]['im'])>0) else 0
                                )
                            )
                            steeringVector[y][j][x]['found']=True



            azim_smoother=10
            for x in range(NB_ELEMENTS_ARRAY):
                for z in range(ELEVATION_CARDINALITY):

                    boundaries=[150,215] # first y last y
                    there_is_sbd=False
                    for y in range(AZIMUTH_CARDINALITY):
                        if (steeringVector[y][z][x]['found']): 
                            there_is_sbd=True
                            break
                    if there_is_sbd:
                        last_y=180
                        for y in range(180,boundaries[0],-1):
                            if (steeringVector[y][z][x]['found']):
                                last_y=y
                                break
                        first_y=180
                        for y in range(180,boundaries[1]):
                            if (steeringVector[y][z][x]['found']):
                                first_y=y
                                break

                        ############### LINEAR BW the two azim bounds and a 180==0 
                        linvals_re1=np.linspace(
                            float(steeringVector[last_y][z][x]['re']) if (len(steeringVector[last_y][z][x]['re'])>0) else 0,
                            float(steeringVector[180][z][x]['re']) if (len(steeringVector[180][z][x]['re'])>0) else 0,
                            180-last_y+1
                        )
                        linvals_re2=np.linspace(
                            float(steeringVector[180][z][x]['re']) if (len(steeringVector[180][z][x]['re'])>0) else 0,
                            float(steeringVector[first_y][z][x]['re']) if (len(steeringVector[first_y][z][x]['re'])>0) else 0,
                            first_y-180+1
                        )
                        linvals_im1=np.linspace(
                            float(steeringVector[last_y][z][x]['im']) if (len(steeringVector[last_y][z][x]['im'])>0) else 0,
                            float(steeringVector[180][z][x]['im']) if (len(steeringVector[180][z][x]['im'])>0) else 0,
                            180-last_y+1
                        )
                        linvals_im2=np.linspace(
                            float(steeringVector[180][z][x]['im']) if (len(steeringVector[180][z][x]['im'])>0) else 0,
                            float(steeringVector[first_y][z][x]['im']) if (len(steeringVector[first_y][z][x]['im'])>0) else 0,
                            first_y-180+1
                        )
                        
                        for j in range(180-last_y+1):
                            steeringVector[last_y+j][z][x]['re']=str(linvals_re1[j])
                            steeringVector[last_y+j][z][x]['im']=str(linvals_im1[j])
                            steeringVector[last_y+j][z][x]['found']=True
                        
                        for j in range(first_y-180+1):
                            steeringVector[180+j][z][x]['re']=str(linvals_re2[j])
                            steeringVector[180+j][z][x]['im']=str(linvals_im2[j])
                            steeringVector[180+j][z][x]['found']=True                            


        for y in range(AZIMUTH_CARDINALITY):
            for z in range(ELEVATION_CARDINALITY):
                for x in range(NB_ELEMENTS_ARRAY):

                    steeringVector[y][z][x]['amp']=((float(steeringVector[y][z][x]['re']) if (len(steeringVector[y][z][x]['re'])>0) else 0)**2+(float(steeringVector[y][z][x]['im']) if (len(steeringVector[y][z][x]['im'])>0) else 0)**2)**(1/2)
                    steeringVector[y][z][x]['psh']=np.arctan2(float(steeringVector[y][z][x]['im']) if (len(steeringVector[y][z][x]['im'])>0) else 0,float(steeringVector[y][z][x]['re']) if (len(steeringVector[y][z][x]['re'])>0) else 0)   
                    steeringVector[y][z][x]['fre']=float(steeringVector[y][z][x]['re']) if (len(steeringVector[y][z][x]['re'])>0) else 0
                    steeringVector[y][z][x]['fim']=float(steeringVector[y][z][x]['im']) if (len(steeringVector[y][z][x]['im'])>0) else 0

                    steeringVector_as_array[y][z][x][0]=float(steeringVector[y][z][x]['re']) if (len(steeringVector[y][z][x]['re'])>0) else 0
                    steeringVector_as_array[y][z][x][1]=float(steeringVector[y][z][x]['im']) if (len(steeringVector[y][z][x]['im'])>0) else 0

        pickle.dump(steeringVector_as_array, open('steeringVector_as_array.dat', 'wb'))
        pickle.dump(steeringVector, open('steeringVector.dat', 'wb'))
    else:
        steeringVector_as_array=pickle.load(open('steeringVector_as_array.dat','rb'))
        steeringVector=pickle.load(open('steeringVector.dat','rb'))



if not (__name__ == '__main__'):
    try:
        steeringVector_as_array=pickle.load(open('steeringVector_as_array.dat','rb'))
        steeringVector=pickle.load(open('steeringVector.dat','rb'))
    except:
        reload_AF_File=True
        update_SV()


### Given a Talon sector expressed with etype dtypes, return a single amplitude value
# (Arbitrarily homecrafted)
def merge_amplitude(etype,dtype):

    etype_mask=[0,5,1,6,4,2,3,7]
    dtype_mask=[0,1,2,4,5,6,3,7]

    minmaxmedian=0.9#0.5#-.795#0.9#6#0.8639988
    amplif=4.5#2.3#29.5 #0.05

    ret=etype_mask.index(etype)
    rdt=dtype_mask.index(dtype)

    if ret*rdt==0:
        return(0.)
    else:
       return( amplif*(minmaxmedian + (( (ret-1)*7+ (rdt-1)  +1)/49)*(1-minmaxmedian) ))

    # return(etype*dtype/49) #linear
    return((2**etype)*(2**dtype)/1384) #binary
    return((2**etype)*(2**dtype)/16384) #binary
    # return(np.log(etype*dtype+1)/np.log(50)) #log    
    return(np.log(2**etype*2**dtype)/np.log(16384)) # logbin

### codebook expression
def prepare_codebook(
    cb_params,                          # set of sectors to express 
    single_elements=False,              # or express sectors with 1 antenna element on only (e.g. to see single elt patterns)
    use_ns3=True,                       # export a txt file to use in NS3, cf. https://github.com/wigig-tools/wigig-module
    tofile=True,                        # export a json file
    json_file_name="temp_local.json"
    ):
    
    global antennas
    
    random_weights=False # a sticking temporary trial

    CB_file="codebook_ap_parametric_GG.txt" # name of ns3 txt file exported

    if single_elements:
        nsectors=32
        sectors_dic={}

        for s in range(nsectors):
            sectors_dic[s]={'id':s+1}
            sectors_dic[s]['elts']=[(0.,0.)]*NB_ELEMENTS_ARRAY
            sectors_dic[s]['elts'][s]=(1.,0.)

    if not single_elements:
        nsectors=len(cb_params)
        sectors_dic={}

        for s in range(nsectors):
            sectors_dic[s]={'id':cb_params[s]['sid']}
            sectors_dic[s]['elts']=[
                (
                    merge_amplitude(cb_params[s]['etype'][x],cb_params[s]['dtype'][int(x/4)]),
                    -cb_params[s]['psh'][x]*np.pi/2
                ) 
                for x in range(NB_ELEMENTS_ARRAY)]


    if not random_weights:
        m_totalAntennas=1
        antennas={
            1:{
                "azimuthOrientationDegree":15, 
                "elevationOrientationDegree":15, 
                "directivity":[0. for x in np.random.permutation(AZIMUTH_CARDINALITY)], #### 1.0?

                ############ SINGLE ELEMENTS WITH WEIGHT = global_sed     
                "nsectors":nsectors,
                # "nsectors":1,
                "sectors":{
                    s+1:{
                        "sectorType":2,
                        "sectorUsage":2,
                        "elements":{
                            x+1:{
                                    "SectorWeights": sectors_dic[s]['elts'][x], #
                                } for x in range(NB_ELEMENTS_ARRAY)
                            }
                        } for s in range(nsectors)
                    },

                "nb_elements":NB_ELEMENTS_ARRAY,
                "phaseQuantizationBits":2,
                "amplitudeQuantizationBits":3,
                "singleElementDirectivity":
                        np.ones((AZIMUTH_CARDINALITY,ELEVATION_CARDINALITY))*global_sed,
                        # [[global_sed for z in np.random.permutation(ELEVATION_CARDINALITY)] for y in np.random.permutation(AZIMUTH_CARDINALITY)], #### other than 1.0?
                "SVA":steeringVector_as_array,
                "elements":{
                    x+1:{
                        "steeringVector":
                            [[
                            (
                                steeringVector[y][z][x]['amp'],steeringVector[y][z][x]['psh']
                            ) 
                            for z in range(ELEVATION_CARDINALITY)] for y in range(AZIMUTH_CARDINALITY)],
                            # "quasiOmniWeights": (x/25.,2*3.14/(27.*(0.001*x+1))),
                        "quasiOmniWeights": (0.0 if x==1 else 0., 0.0 if x==1 else 0.) #### 1.0? #### 1.0?

                        } for x in range(NB_ELEMENTS_ARRAY)
                    }
                },
            }

    else:
        m_totalAntennas=2
        antennas={
            1:{
                "azimuthOrientationDegree":15, 
                "elevationOrientationDegree":15, 
                "directivity":[x/10. for x in np.random.permutation(AZIMUTH_CARDINALITY)],
                "nsectors":3,
                "sectors":{
                    3:{
                        "sectorType":2,
                        "sectorUsage":2,
                        "elements":{
                            x+1:{
                                "SectorWeights": (x/25.,2*3.14/(27.*(0.001*x+1))),
                                } for x in np.random.permutation(NB_ELEMENTS_ARRAY)
                            }
                        },
                    1:{"sectorType":2,
                        "sectorUsage":2,
                        "elements":{
                            x+1:{
                                "SectorWeights": (x/25.,2*3.14/(27.*(0.001*x+1))),
                                } for x in np.random.permutation(NB_ELEMENTS_ARRAY)
                            }
                        },
                    2:{"sectorType":2,
                        "sectorUsage":2,
                        "elements":{
                            x+1:{
                                "SectorWeights": (x/25.,2*3.14/(27.*(0.001*x+1))),
                                } for x in np.random.permutation(NB_ELEMENTS_ARRAY)
                            }
                        },
                    },
                "nb_elements":NB_ELEMENTS_ARRAY,
                "phaseQuantizationBits":2,
                "amplitudeQuantizationBits":3,
                "singleElementDirectivity":
                        [[z/20. for z in np.random.permutation(ELEVATION_CARDINALITY)] for y in np.random.permutation(AZIMUTH_CARDINALITY)],
                "elements":{
                    x+1:{
                        "steeringVector":
        [[(z/25.,2*3.14/(27.*(0.001*z+1))) for z in np.random.permutation(ELEVATION_CARDINALITY)] for y in np.random.permutation(AZIMUTH_CARDINALITY)],
                        "quasiOmniWeights": (x/25.,2*3.14/(27.*(0.001*x+1))),

                        } for x in np.random.permutation(NB_ELEMENTS_ARRAY)
                    }
                },
            2:{
                "azimuthOrientationDegree":155,
                "elevationOrientationDegree":155,
                "directivity":[x/5. for x in np.random.permutation(AZIMUTH_CARDINALITY)],
                "nsectors":2,
                "sectors":{
                    2:{
                        "sectorType":2,
                        "sectorUsage":2,
                        "elements":{
                            x+1:{
                                "SectorWeights": (x/25.,2*3.14/(27.*(0.001*x+1))),
                                } for x in np.random.permutation(NB_ELEMENTS_ARRAY)
                            }
                        },
                    1:{"sectorType":2,
                        "sectorUsage":2,
                        "elements":{
                            x+1:{
                                "SectorWeights": (x/25.,2*3.14/(27.*(0.001*x+1))),
                                } for x in np.random.permutation(NB_ELEMENTS_ARRAY)
                            }
                        }
                    },
                "nb_elements":NB_ELEMENTS_ARRAY,
                "phaseQuantizationBits":2,
                "amplitudeQuantizationBits":3,
                "singleElementDirectivity":
                        [[z/20. for z in np.random.permutation(ELEVATION_CARDINALITY)] for y in np.random.permutation(AZIMUTH_CARDINALITY)],
                "elements":{
                    x+1:{
                        "steeringVector":
        [[(z/25.,2*3.14/(27.*(0.001*z+1))) for z in np.random.permutation(ELEVATION_CARDINALITY)] for y in np.random.permutation(AZIMUTH_CARDINALITY)],
                        "quasiOmniWeights": (x/25.,2*3.14/(27.*(0.001*x+1))),

                        } for x in np.random.permutation(NB_ELEMENTS_ARRAY)
                    }
                }
            }


    if(use_ns3):
        file = open(CB_file,"w") 
        file.write(str(m_totalAntennas)+'\n')
        for ant_id in antennas:
            file.write(str(ant_id)+'\n')
            file.write(str(antennas[ant_id]["azimuthOrientationDegree"])+'\n')
            file.write(str(antennas[ant_id]["elevationOrientationDegree"])+'\n')    
            file.write(str(antennas[ant_id]["nb_elements"])+'\n')    
            file.write(str(antennas[ant_id]["phaseQuantizationBits"])+'\n')    
            file.write(str(antennas[ant_id]["amplitudeQuantizationBits"])+'\n')    
            for i in range(AZIMUTH_CARDINALITY):
                # file.write(str(antennas[ant_id]["singleElementDirectivity"][i])+'\n')
                file.write(''.join(str(antennas[ant_id]["singleElementDirectivity"][i][j])+',' for j in range(ELEVATION_CARDINALITY))[:-1]+'\n')
            for x in range(NB_ELEMENTS_ARRAY):
                for i in range(AZIMUTH_CARDINALITY):
                    coltojoin=[
                        str(antennas[ant_id]["elements"][x+1]["steeringVector"][i][j][0])+','
                        +str(antennas[ant_id]["elements"][x+1]["steeringVector"][i][j][1])
                        +','
                          for j in range(ELEVATION_CARDINALITY)]
                    # if x==3 and i>350:
                    #     print(coltojoin[100])

                    str_to_wr=''.join(coltojoin)

                    # if x==3 and i>350:
                    #     print(str_to_wr[-1])

                    str_to_wr=str_to_wr[:-1]
                    str_to_wr+='\n'
                    file.write(str_to_wr)  
                    
            file.write(''.join(str(antennas[ant_id]["elements"][x+1]["quasiOmniWeights"][0])+','
                    +str(antennas[ant_id]["elements"][x+1]["quasiOmniWeights"][1])+','
                        for x in range(NB_ELEMENTS_ARRAY))[:-1]+'\n')         

            file.write(str(antennas[ant_id]["nsectors"])+'\n')
            for sectorID in (antennas[ant_id]["sectors"]):
                file.write(str(sectorID)+'\n')
                file.write(str(antennas[ant_id]["sectors"][sectorID]["sectorType"])+'\n')
                file.write(str(antennas[ant_id]["sectors"][sectorID]["sectorUsage"])+'\n')
                file.write(''.join(
                    str(antennas[ant_id]["sectors"][sectorID]["elements"][x+1]["SectorWeights"][0])+','
                    +str(antennas[ant_id]["sectors"][sectorID]["elements"][x+1]["SectorWeights"][1])+',' 
                    for x in range(NB_ELEMENTS_ARRAY))[:-1]+'\n')
        file.close() 
        return({})

    else: 
    
        json_CB={}
        for ant_id in antennas:
            json_CB[ant_id]={
                "Orientation":antennas[ant_id]["azimuthOrientationDegree"],
                "Num_Elements":antennas[ant_id]["nb_elements"],
                "Num_Sectors":antennas[ant_id]["nsectors"],
                "Antenna_QO_Directivities":computeQOD(ant_id) if tofile else {},
                "Sectors":{}
            }
            json_CB[ant_id]["Sectors"]=compute_par_SDs(ant_id,tofile)

        if tofile:
            jsonFile = open(json_file_name,"w")
            json.dump(json_CB, jsonFile)
            jsonFile.close()
        return(json_CB)

## sector directivities
def compute_par_SDs(ant_id,tofile):
    
    list_sector_ids=[]
    for sectorID in (antennas[ant_id]["sectors"]):
        list_sector_ids.append(sectorID)

    if allowed_multi_pool:

        print(cpu_count())
        # with Pool(processes=cpu_count()-2 or 1) as p:
        with Pool(processes=4) as p:       
            if tofile: 
                resdirs=p.map(computeSDP, list_sector_ids)
            else: 
                resdirs=p.map(computeSDP_array, list_sector_ids)

    else:
        resdirs=[]
        if tofile:
            for indsec in list_sector_ids:
                resdirs.append(computeSDP(indsec))
        else:
            for indsec in list_sector_ids:
                resdirs.append(computeSDP_array(indsec))


    # print(resdirs)
    zesectors={}
    ind=0
    for sectorID in (antennas[ant_id]["sectors"]):
        zesectors[str(sectorID)]={"Sector_Directivities":resdirs[ind][ant_id]}
        ind+=1

    return(zesectors)

## old one showing similar expression as in the C++ code
def CalculateDirectivityold (WeightsVector,antennas,ant_id):
    #singleElementDirectivity,steeringVector):
    directivity={}
    for m in range(AZIMUTH_CARDINALITY):
        directivity[m]={}
        for l in range(ELEVATION_CARDINALITY):
            value=0
            j=1
            for weight in WeightsVector:
                steeringVectorJ=antennas[ant_id]["elements"][j]["steeringVector"][m][l][0]+antennas[ant_id]["elements"][j]["steeringVector"][m][l][1]*1j
                value+=(weight[0]+weight[1]*1j)*steeringVectorJ
                j+=1
            value*=antennas[ant_id]["singleElementDirectivity"][m][l]
            directivity[m][l]=10*np.log10(abs(value))
            if (directivity[m][l]<-20.0):
                directivity[m][l]=-20.0

    return(directivity)

### quasi omni dirs
def computeQOD(ant_id):
    # weights_t=[]
    weights=np.zeros((NB_ELEMENTS_ARRAY,2))
    for x in range(NB_ELEMENTS_ARRAY):
        amp=antennas[ant_id]["elements"][x+1]["quasiOmniWeights"][0]
        phase=antennas[ant_id]["elements"][x+1]["quasiOmniWeights"][1]
        real=amp*np.cos(phase)
        imag=amp*np.sin(phase)
        # weights_t.append([real,imag])
        weights[x]=[real,imag]

    # weights=np.array(weights_t)

    return (CalculateDirectivity(weights,ant_id))

### sector dir
def computeSD(ant_id,sec_id):
    # weights_t=[]
    weights=np.zeros((NB_ELEMENTS_ARRAY,2))
    for x in range(NB_ELEMENTS_ARRAY):
        amp=antennas[ant_id]["sectors"][sec_id]["elements"][x+1]["SectorWeights"][0]
        phase=antennas[ant_id]["sectors"][sec_id]["elements"][x+1]["SectorWeights"][1]
        real=amp*np.cos(phase)
        imag=amp*np.sin(phase)
        # weights_t.append([real,imag])
        weights[x]=[real,imag]

    # weights=np.array(weights_t)


    return (CalculateDirectivity(weights,ant_id))

### same par
def computeSDP(sec_id):

    res={}
    for ant_id in antennas:
        # weights_t=[]
        weights=np.zeros((NB_ELEMENTS_ARRAY,2))
        for x in range(NB_ELEMENTS_ARRAY):
            amp=antennas[ant_id]["sectors"][sec_id]["elements"][x+1]["SectorWeights"][0]
            phase=antennas[ant_id]["sectors"][sec_id]["elements"][x+1]["SectorWeights"][1]
            real=amp*np.cos(phase)
            imag=amp*np.sin(phase)
            # weights_t.append([real,imag])
            weights[x]=[real,imag]

        # weights=np.array(weights_t)

        res[ant_id]=CalculateDirectivity(weights,ant_id)

    return (res)

### same parallel and numpyfied
def computeSDP_array(sec_id):


    res={}

    for ant_id in antennas:

        ######################### NEW
        minaray=(np.full((AZIMUTH_CARDINALITY,ELEVATION_CARDINALITY), -20.0))

        complexarray=np.zeros((AZIMUTH_CARDINALITY,ELEVATION_CARDINALITY),dtype='complex128')
        for x in range(NB_ELEMENTS_ARRAY):
            complexarray+=(
                (steeringVector_as_array[:,:,x,0]+steeringVector_as_array[:,:,x,1]*1j)*
                (1 if antennas[ant_id]["sectors"][sec_id]["elements"][x+1]["SectorWeights"][0]>0 else 0)*
                np.exp((antennas[ant_id]["sectors"][sec_id]["elements"][x+1]["SectorWeights"][1])*(-1j))
                )

        res[ant_id]=np.maximum(10*np.log10(np.maximum(0.01,abs(complexarray))),minaray)
            

        # ######################### OLD
        # # weights_t=[]
        # weights=np.zeros((NB_ELEMENTS_ARRAY,2))
        # for x in range(NB_ELEMENTS_ARRAY):
        #     amp=antennas[ant_id]["sectors"][sec_id]["elements"][x+1]["SectorWeights"][0]
        #     phase=antennas[ant_id]["sectors"][sec_id]["elements"][x+1]["SectorWeights"][1]
        #     real=amp*np.cos(phase)
        #     imag=amp*np.sin(phase)
        #     # weights_t.append([real,imag])
        #     weights[x]=[real,imag]

        # # weights=np.array(weights_t)

        # res[ant_id]=CalculateDirectivity(weights,ant_id,as_array=True)

    return (res)

### numpyfied with bound -20dB
def CalculateDirectivity (WeightsVector,ant_id,as_array=False):

    steeringVector_as_array=antennas[ant_id]["SVA"]
    # print(np.nonzero(steeringVector_as_array))
    singleElementDirectivity_as_array=np.array(antennas[ant_id]["singleElementDirectivity"])
    dirsarray=np.zeros((AZIMUTH_CARDINALITY,ELEVATION_CARDINALITY))
    complexarray=np.zeros((AZIMUTH_CARDINALITY,ELEVATION_CARDINALITY),dtype='complex128')
    minaray=(np.full((AZIMUTH_CARDINALITY,ELEVATION_CARDINALITY), -20.0))

    # complexarray+=(WeightsVector[:,0]+WeightsVector[:,1]*1j)*(steeringVector_as_array[:,:,:,0]+steeringVector_as_array[:,:,:,1]*1j)

    for x in range(NB_ELEMENTS_ARRAY):
        complexarray+=(WeightsVector[x][0]+WeightsVector[x][1]*1j)*(steeringVector_as_array[:,:,x,0]+steeringVector_as_array[:,:,x,1]*1j)

    complexarray*=singleElementDirectivity_as_array
    disarray=np.maximum(10*np.log10(np.maximum(0.01,abs(complexarray))),minaray)

    if as_array:
        return(disarray)

    directivity={}
    for m in range(AZIMUTH_CARDINALITY):
        directivity[str(m)]={}
        for l in range(ELEVATION_CARDINALITY):
            directivity[str(m)][str(l)]=disarray[m][l]

    return(directivity)          

if __name__ == '__main__':
    

    reload_AF_File=True
    show_AF()
    update_SV(complete_packmen=True)

    # prepare_codebook({},single_elements=True,use_ns3=True)
    # prepare_codebook({},single_elements=True,use_ns3=False,tofile=True)

