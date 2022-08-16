import json
import os
import tqdm
import psycopg2 

def loadDataFile(cnx, path, tbl):
    fQuery = f"COPY {tbl} FROM \'{path}\' DELIMITER \',\' csv header;"
    cursor = cnx.cursor()
    cursor.execute(fQuery)
    cnx.commit()


cnx_a = psycopg2.connect(host="localhost", user="postgres", password="CitrusABM21", database="citrus")
cursor = cnx_a.cursor()

basepath = "/home/instr1/repo/EconABM_accounting/EconABM/output/initial_info"
filesList = []
for subfolder in os.listdir(basepath):
    for file in os.listdir(f"{basepath}/{subfolder}"):
        tbl = ""
        if "econ" in file:
            tbl = "econ"
        else:
            tbl = "bio"
        filesList.append([f"{basepath}/{subfolder}/{file}", tbl])

for file in tqdm.tqdm(filesList):
    loadDataFile(cnx_a, file[0], file[1])

cnx_a.close()

