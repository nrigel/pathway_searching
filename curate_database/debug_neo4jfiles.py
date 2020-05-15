debug_row = 31

import csv, json

with open('pathbank_data/neo4jnodes_Pathway.csv') as csvin:
    csvin.readline()
    row_counter = -1
    for row in csv.reader(csvin):
        row_counter += 1
        if row_counter == debug_row:
            print(row)
            exit()