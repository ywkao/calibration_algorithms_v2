#!/usr/bin/env python2

#import pandas as pd
#df1 = pd.read_csv(r'./meta_conditions/output_EDAnalyzer_pedestals.csv')
#df2 = pd.read_csv(r'./meta_conditions/output_EDAnalyzer_cm_parameters_ped_subtracted.csv')
#print( df1 )

output_csv_file = "./meta_conditions/calibration_parameters.csv"

f1 = open("./meta_conditions/output_EDAnalyzer_pedestals.csv", "r")
f2 = open("./meta_conditions/output_EDAnalyzer_cm_parameters_ped_subtracted.csv", "r")

def retrieve_data(f):
    data = {}
    for i, line in enumerate(f.readlines()):
        if "columns" in line:
            columns = line.strip().split(': ')[1]
            items = columns.split(', ')
            print items
            for item in items: data[item] = []
            continue
        elif "#" in line: continue
        
        row = line.strip().split(',')
        for k in range(len(items)):
            data[items[k]].append(row[k])
    
    return data

def print_header(fout, keys):
    output = "# "
    for i, key in enumerate(keys):
        output += key
        if not i+1==len(keys): output += ", "
        else: output += "\n"

    fout.write("#" + "-"*140 + "\n")
    fout.write(output)
    fout.write("#" + "-"*140 + "\n")

def create_merged_csv_file(data):
    keys = ["channel", "pedestal", "slope", "intercept", "kappa", "charge", "residual_offset", "conversion_ADC_to_fC", "conversion_ToT_to_fC", "conversion_ToA_to_ns"]

    initial_values = [0.] * len(data["channel"])

    with open(output_csv_file, 'w') as fout:
        print_header(fout, keys)
        for i in range(len(data["channel"])):
            row = "%s,%s,%s,%s,0.,0.,0.,0.,0.,0.\n" % (
                    data["channel"][i],
                    data["pedestal"][i],
                    data["slope"][i],
                    data["intercept"][i]
                  )
            fout.write(row)

if __name__ == "__main__":
    d1 = retrieve_data(f1)
    d2 = retrieve_data(f2)

    d3 = d1.copy()
    for key, value in d2.items(): d3[key] = value

    create_merged_csv_file(d3)
