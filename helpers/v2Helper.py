import json
import os
from pickle import FALSE
import psycopg2
import multiprocessing as mp
import tqdm
import random
import string
import math
def id_generator(size=30, chars=string.ascii_uppercase + string.digits):
    return ''.join(random.choice(chars) for _ in range(size))

add_cluster = ("INSERT INTO simulation_clusters "
               "(cluster_id, description) "
               "VALUES (DEFAULT, %(descrip)s) RETURNING cluster_id;")

def clusterExists(cnx, cursor, description):
    cursor.execute("SELECT cluster_id FROM simulation_clusters WHERE description=%s", (description,))
    row = cursor.fetchone()
    if row is not None:
        return row[0]
    else:
        return False

def createCluster(cnx, cursor, description):
    #check if cluster exists
    checkExist = clusterExists(cnx, cursor, description)
    if checkExist is False:
        cluster_data = { 'descrip': description}
        cursor.execute(add_cluster, cluster_data)
        cnx.commit()
        return cursor.fetchone()[0]
    else:
        return checkExist


### CREATE SET
add_set = ("INSERT INTO simulation_sets"
           "(set_id, description, cluster_id)"
           "VALUES (DEFAULT, %(description)s, %(cluster_id)s) RETURNING set_id;")

def setExists(cnx, cursor, description):
    cursor.execute("SELECT set_id FROM simulation_sets WHERE description=%s", (description,))
    row = cursor.fetchone()
    if row is not None:
        return row[0]
    else:
        return False
def createSet(cnx, cursor, description, clusterid):
    checkExist = setExists(cnx, cursor, description)
    if checkExist is False:
        set_data = {'description': description, 'cluster_id': clusterid}
        cursor.execute(add_set, set_data)
        cnx.commit()
        return cursor.fetchone()[0]
    else:
        return checkExist

## CREATE EXPERIMENT
add_experiment = ("INSERT INTO experiments VALUES"
                  "(DEFAULT,"
                  "%(fresh_yield)s, %(juice_yield)s, %(fruit_price)s, %(juice_price)s, %(fixed_costs)s, %(bio_param_id)s, %(set_id)s,"
                  "%(invasionDays)s, %(invasionModalities)s)"
                  "RETURNING experiment_id;"
)

def createExperiment(cnx, cursor, paramList):
    bio_id = 1 #default bio (haven't changed it yet)
    experiment_data = {
        'fresh_yield': paramList['fresh_yield'],
        'juice_yield': paramList['juice_yield'],
        'fruit_price': paramList['fruit_price'],
        'juice_price': paramList['juice_price'],
        'fixed_costs': paramList['fixed_costs'],
        'bio_param_id': paramList['bio_param_id'],
        'set_id': paramList['set_id'],
        'invasionDays': paramList['invasionDays'],
        'invasionModalities': paramList['invasionModalities']
    }
    cursor.execute(add_experiment, experiment_data)
    cnx.commit()
    return cursor.fetchone()[0]

def loadDataFile(cnx, path, tbl):
    fQuery = f"copy {tbl} FROM \'{path}\' DELIMITER \',\' csv header;"
    cursor = cnx.cursor()
    cursor.execute(fQuery)
    cnx.commit()

cnx_a = psycopg2.connect(host="localhost", user="postgres", password="CitrusABM21", database="netbenefits")
cursor = cnx_a.cursor()

cluster_name = "TESTER"
cluster_id = createCluster(cnx_a,cursor,cluster_name)
global_bioid = 1


#Open default config file
config_file = None 

with open("configs/bioConfig.json", "r") as read_file:
    config_file = json.load(read_file)

#Create a list of bio configs
import datetime
import copy 
configList = []
numYears = math.floor(config_file["modelDuration"] / 365)
for numInvasions in ["initial","half","all","double"]:
    for _ in range(0,50):
        filePrefix = f'{datetime.datetime.now()}'.replace(' ', '-').replace(':','.')
        config = json.loads(json.dumps(config_file))
        invasionDays = ""
        invasionDays_lst = []
        if numInvasions == "initial":
            invasionDays_lst = ["81"]
        elif numInvasions == "half":
            invasionDays_lst = [str((x)*365 + 81) for x in range(0,numYears) if ((x) % 2 == 0 or (x)==0)]
        elif numInvasions == "all":
            invasionDays_lst = [str((x)*365+81) for x in range(0,numYears)]
        elif numInvasions == "double":
            invasionDays_lst = [str((x)*365+81) for x in range(0,numYears)] + [str((x)*365+251) for x in range(0,numYears)]
        config["invasionDays"] = ','.join([x for x in invasionDays_lst])
        invasionModality_lst = random.choices(range(1,7),k=len(invasionDays_lst))
        invasionModalities = ','.join([str(x) for x in invasionModality_lst])
        config["invasionModalities"] = invasionModalities
        config["fileName"] = f'{filePrefix}_bio.csv'
        configList.append(config)


econConfig = None 
with open("configs/econConfig.json", "r") as read_file:
    econConfig = json.load(read_file)




#List to hold all of our triples
masterConfigList = []

def getExperimentID(cnx, setid, econConfig, bioConfig, bioid):
    paramList = {
        'fresh_yield': econConfig['freshYield'],
        'juice_yield': econConfig['juiceYield'],
        'fruit_price': econConfig['freshPrice'],
        'juice_price': econConfig['juicePrice'],
        'fixed_costs': econConfig['costs'],
        'invasionModalities': bioConfig['invasionModalities'],
        'invasionDays': bioConfig['invasionDays'],
        'bio_param_id': bioid,
        'set_id': setid,
    }
    cursor = cnx.cursor()
    eID = createExperiment(cnx,cursor, paramList)
    cnx.commit()
    return eID

setid_dict = {}
def appendMasterTriples(cnx, econ, biocons, folder, master, numSets=0, nameSuffix=""):
    if numSets==0:
        rndBio = random.randint(0, len(biocons)-1)
        bio_copy = copy.deepcopy(biocons[rndBio])
        econ_copy = copy.deepcopy(econ) 
        rndmid = id_generator()
        econ_copy["outputFilename"] = f"{rndmid}_econ.csv"
        bio_copy["fileName"] = econ_copy["outputFilename"].replace("econ", "bio")
        master.append([econ_copy, bio_copy, folder])
    else: 
        cursor = cnx.cursor()
        setid = createSet(cnx, cursor, folder, cluster_id)
        setid_dict[folder] = setid
        for ns in range(numSets):
            for bio in biocons:
                econ_copy = copy.deepcopy(econ)
                bio_copy = copy.deepcopy(bio)
                rndmid = id_generator()
                bio_copy["fileName"] = f"{rndmid}_bio.csv"
                econ_copy["outputFilename"] = bio_copy["fileName"].replace("bio", "econ")
                master.append([econ_copy, bio_copy, folder])


#One set
# setid = createSet(cnx_a, cursor, "all", cluster_id)
# setid_dict["all"] = setid
# #Randomize
# for _ in range(0,2000):
#     econ_copy = copy.deepcopy(econConfig)
#     for i in range(0,3):
#         for j in range(0,3):
#             econ_copy[f"g{i}{j}_lambda"] = random.random()
#             econ_copy[f"g{i}{j}_alpha"] = random.random() 
#     rndm_bio_idx = random.randint(0, len(configList) - 1)
#     bio_copy = copy.deepcopy(configList[rndm_bio_idx])
#     masterConfigList.append([econ_copy, bio_copy, "all"])

#REDOING ALL
# for lmbda in [0,0.25,0.5,0.75,1]:
#     for alpha in [0,0.25,0.5,0.75,1]:
#         econ_copy = copy.deepcopy(econConfig)
#         for i in range(0,3):
#             for j in range(0,3):
#                 econ_copy[f"g{i}{j}_alpha"] = alpha
#                 econ_copy[f"g{i}{j}_lambda"] = lmbda
#         appendMasterTriples(cnx_a,econ_copy, configList, f"l{int(lmbda*100)}_a{int(alpha*100)}", masterConfigList)


#PARAMTER FORMAT: removalCost, surveyCost, frequency, radius, efficacy, spraycost, denseCosts, 
#                 yield multiplier, removalcost, surveycost, frequency, width, height
#FIRST NEW TEST
# for e in [600,700,800,900]:
#     econ_copy = copy.deepcopy(econConfig)
#     econ_copy["strategyFlags"] = "0,1,0"
#     econ_copy["strategyParameters"] = f"5,5,0,0,{e/1000},5,5,0"
#     appendMasterTriples(cnx_a, econ_copy, configList, f"sprayVariations",  masterConfigList, numSets=1)
# # # #ROGUE TEST
# for radius in [20,40,60]:
#     for frequency in [1,6,12,18,24]:
#         econ_copy = copy.deepcopy(econConfig)
#         econ_copy["strategyFlags"] = "1,0,0"
#         econ_copy["strategyParameters"] = f"5,5,{frequency},{radius},0,5,5,0"
#         appendMasterTriples(cnx_a, econ_copy, configList, f"rogueVariations", masterConfigList, numSets=1)
# # # # #SPRAY AND ROGUE
# for e in [600,700,800,900]:
#     for frequency in [1,6,12,18,24]:
#         for radius in [1,8,10,20,40]:
#             econ_copy = copy.deepcopy(econConfig)
#             econ_copy["strategyFlags"] = "1,1,0"
#             econ_copy["strategyParameters"] = f"5,5,{frequency},{radius},{e/1000},5,5,0"
#             appendMasterTriples(cnx_a, econ_copy, configList, f"rogueSprayVariations", masterConfigList, numSets=1)

# RECTANGULAR ROGUE
for frequency in [1,6,18]:
    for width in [1,8,20,40]:
        for height in [1,8,20,40]:
            if width == height:
                continue
            econ_copy = copy.deepcopy(econConfig)
            econ_copy["strategyFlags"] = "0,0,0,1"
            econ_copy["strategyParameters"] = f"0,0,0,0,0,0,0,0,5,0,{frequency},{width},{height}"
            appendMasterTriples(cnx_a, econ_copy, configList, f"recRogueVariations", masterConfigList, numSets=1)

# RECTANGULAR ROGUE AND SPRAY
for e in [600,700,800,900]:
    for frequency in [1,6,18]:
        for width in [1,8,20,40]:
            for height in [1,8,20,40]:
                if width == height:
                    continue
                econ_copy = copy.deepcopy(econConfig)
                econ_copy["strategyFlags"] = "0,1,0,1"
                econ_copy["strategyParameters"] = f"0,0,0,0,{e/1000},0,0,0,5,0,{frequency},{width},{height}"
                appendMasterTriples(cnx_a, econ_copy, configList, f"recRogueSprayVariations", masterConfigList, numSets=1)

# # # # DENSE PLANTING
# for density in [20,80,140,200,220]:
#     econ_copy = copy.deepcopy(econConfig)
#     ym = 0.5182*(density/100) + 0.0761
#     econ_copy["strategyFlags"] = "0,0,1"
#     econ_copy["strategyParameters"] = f"5,5,0,0,0,5,5,{ym}"
#     appendMasterTriples(cnx_a, econ_copy, configList, "densePlantingVariations", masterConfigList, numSets=1)

# # # # DENSE PLANTING AND ROGUE
# for density in [20,80,140,200,220]:
#     for frequency in [1,6,12,18,24]:
#         econ_copy = copy.deepcopy(econConfig)
#         ym = 0.5182*(density/100) + 0.0761
#         econ_copy["strategyFlags"] = "0,0,1"
#         econ_copy["strategyParameters"] = f"5,5,{frequency},0,0,5,5,{ym}"
#         appendMasterTriples(cnx_a, econ_copy, configList, "densePlantingVariations", masterConfigList, numSets=1)


# # # # DENSE PLANTING AND SPRAY
# for e in [600,700,800,900]:
#     for density in [20,80,140,200,220]:
#         econ_copy = copy.deepcopy(econConfig)
#         ym = 0.5182*(density/100) + 0.0761
#         econ_copy["strategyFlags"] = "0,1,1"
#         econ_copy["strategyParameters"] = f"5,5,0,0,{e/1000},5,5,{ym}"
#         appendMasterTriples(cnx_a, econ_copy, configList, "denseAndSprayVariations", masterConfigList, numSets=1)

# # # DENSE PLANTING AND SPRAY AND ROGUE
# for density in [20,80,140,200,220]:
#     for e in [600,700,800,900]:
#         for frequency in [1,6,12,18,24]:
#                 econ_copy = copy.deepcopy(econConfig)
#                 ym = 0.5182*(density/100) + 0.0761
#                 econ_copy["strategyFlags"] = "1,1,1"
#                 econ_copy["strategyParameters"] = f"5,5,{frequency},0,{e/1000},5,5,{ym}"
#                 appendMasterTriples(cnx_a, econ_copy, configList, f"denseRogueSprayVariations", masterConfigList, numSets=1)
# #NO ACTION BASELINE
# econ_copy = copy.deepcopy(econConfig)
# econ_copy["strategyFlags"] = "0,0,0"
# econ_copy["strategyParameters"] = f"0,0,0,0,0,0,0,0"
# appendMasterTriples(cnx_a, econ_copy, configList, f"noAction_baseCase",  masterConfigList, numSets=10)

def runInstance(config):
    biofileName = config[1]["fileName"][:-4]
    econFileName = config[0]["outputFilename"][:-4]
    os.makedirs(f'configs/{cluster_name}/{config[2]}', exist_ok=True)
    os.makedirs(f'output/{cluster_name}/{config[2]}', exist_ok=True)
    expid = config[3]
    config[0]["outputFilename"] = f"output/{cluster_name}/{config[2]}/{econFileName}.csv"
    config[1]["fileName"] = f"output/{cluster_name}/{config[2]}/{biofileName}.csv"
    config[0]["experimentID"] = expid
    with open(f"configs/{cluster_name}/{config[2]}/{biofileName}.json", "w") as write_file:
        json.dump(config[1], write_file)
    with open(f"configs/{cluster_name}/{config[2]}/{econFileName}.json", "w") as write_file:
        json.dump(config[0], write_file)   
    os.system(f"./serial.out \"configs/{cluster_name}/{config[2]}/{econFileName}.json\" \"configs/{cluster_name}/{config[2]}/{biofileName}.json\"")


revisedMaster = []
for config in masterConfigList:
    expid = getExperimentID(cnx_a, setid_dict[config[2]], config[0], config[1], global_bioid)
    newConfig = copy.deepcopy(config)
    newConfig.append(expid)
    revisedMaster.append(newConfig)



pool = mp.Pool(mp.cpu_count())

for _ in tqdm.tqdm(pool.imap_unordered(runInstance, revisedMaster), total=len(revisedMaster)):
    pass 

pool.close()

basePath = f"/home/instr1/repo/EconABM_accounting/EconABM/output/{cluster_name}"
print("Starting data loading...")
for config in tqdm.tqdm(revisedMaster):
    filename = f"{basePath}/{config[2]}/{config[0]['outputFilename']}"
    loadDataFile(cnx_a, filename, "econ")
    os.remove(filename)
    #biofilename = f"{basePath}/{config[2]}/{config[1]['fileName']}"
    #loadDataFile(cnx_a, biofilename, "bio")
    #os.remove(biofilename)




cnx_a.close()



