import os
import psycopg2
import tqdm
def loadDataFile(cnx, path, tbl):
    fQuery = f"copy {tbl} FROM \'{path}\' DELIMITER \',\' csv header;"
    cursor = cnx.cursor()
    cursor.execute(fQuery)
    cnx.commit()

cnx_a = psycopg2.connect(host="localhost", user="postgres", password="CitrusABM21", database="citrus")
cursor = cnx_a.cursor()

folder = "/home/instr1/repo/EconABM_accounting/EconABM/output/randomized_runs2"

for subfolder in os.listdir(folder):
    for file in tqdm.tqdm(os.listdir(f"{folder}/{subfolder}")):
        if "econ" in file:
            loadDataFile(cnx_a, f"{folder}/{subfolder}/{file}", "econ")            
