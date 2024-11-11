import scipy.io
import matplotlib.pyplot as plt
import numpy as np
import itertools
import copy
import pprint

import csv

"""
Measure the sum of the resistance vector at every point and plot that as well.
"""


def main():
    mat = scipy.io.loadmat('HeatmapWorkspace.mat')

    TD_Heatmaps = np.zeros((20,5,4))
    TD_Heatmaps[:,:,0] = mat['D0D2'][:,:]
    TD_Heatmaps[:,:,1] = mat['D2D4'][:,:]
    TD_Heatmaps[:,:,2] = mat['D4D6'][:,:]
    TD_Heatmaps[:,:,3] = mat['D6D8'][:,:]

    possible_drugs = [0, 1, 2, 3, 4]

    resistance_final = []
    growth_final = []
    totRes_final = []

        
    #for dose_seq in itertools.combinations_with_replacement(possible_drugs, 4):
    countuh = 0
    for dose_seq in list(itertools.product(range(5), repeat=4)):
        doseCount = [0, 0, 0, 0, 0]
        
        resistance_vector = []
        resistance_vector_new = []

        growth_vector = []
        totRes_vector = [[0],[0],[0],[0]]

        switch_counter = 0

        print(countuh,' ',dose_seq)
        countuh += 1
        for dose in dose_seq:
            resistance_vector = copy.deepcopy(resistance_vector_new)
            resistance_vector_new = []

            if (sum(doseCount) == 0):
                growth_vector = [[0],[0],[0],[0]]
                for i in range(4):
                    resistance_vector_new.append(TD_Heatmaps[dose*4+i, : ,0]) 
                
                for i in range(4):
                    totRes_vector[i].append(np.sum(resistance_vector_new[i]))

                doseCount[dose] += 1
                
            else:
                if doseCount[dose] == 0:
                    #Growth Vector Stuff
                    growth_vector_new = [[] for x in range(len(growth_vector)*4)]

                    for j in range(len(growth_vector)):
                        for i in range(4):
                            growth_vector_new[j*4 + i] = copy.deepcopy(growth_vector[j])        

                    for j in range(len(growth_vector_new)):
                        growth_vector_new[j].append(resistance_vector[j % 4][dose])

                    growth_vector = copy.deepcopy(growth_vector_new) 

                    # Resistance Vector Stuff
                    for j in range(len(resistance_vector)):
                        for i in range(4):
                            resistance_vector_new.append(resistance_vector[j] + TD_Heatmaps[dose*4+i, :, 0])

                    # Total Resistance Stuff
                    totRes_vector_new = [[] for x in range(len(totRes_vector)*4)]

                    for j in range(len(totRes_vector)):
                        for i in range(4):
                            totRes_vector_new[j*4 + i] = copy.deepcopy(totRes_vector[j])
                    
                    for j in range(len(totRes_vector_new)):
                        totRes_vector_new[j].append(np.sum(resistance_vector_new[j]))

                    totRes_vector = copy.deepcopy(totRes_vector_new)

                else:
                    #Growth Vector Stuff
                    for j in range(len(growth_vector)):
                        growth_vector[j].append(resistance_vector[j][dose])

                    #Resistance Vector Stuff
                    for j in range(len(resistance_vector)):
                        resistance_vector_new.append(resistance_vector[j] + TD_Heatmaps[dose*4 + (j % 4), :, doseCount[dose]])

                    #Total resistance stuff
                    for j in range(len(totRes_vector)):
                        totRes_vector[j].append(np.sum(resistance_vector_new[j]))

                doseCount[dose] += 1

        resistance_final.append(resistance_vector_new)
        growth_final.append(growth_vector)
        totRes_final.append(totRes_vector)   

    dummy = 0
    for i in range(len(growth_final)):
        dummy += len(growth_final[i])
    print(dummy)

    growth_final_norm = copy.deepcopy(growth_final)
    totRes_final_norm = copy.deepcopy(totRes_final)

    for i in range(len(growth_final_norm)):
        for j in range(len(growth_final_norm[i])):
            growth_final_norm[i][j][1] = sum(growth_final[i][j][0:2])
            growth_final_norm[i][j][2] = sum(growth_final[i][j][0:3])
            growth_final_norm[i][j][3] = sum(growth_final[i][j][0:4])


    
    
    """
    Calculate the number of unique drugs used in the sequence so we can plot by drugs used.
    """
    unq_drugs = []
    for dose_seq in list(itertools.product(range(5), repeat=4)):
        dose_seq_set = set(dose_seq)
        unq_drugs.append(len(dose_seq_set))
        
    print(unq_drugs)
    
    plt.figure(figsize=(14,7))
    
    two_drugs_growth = []
    three_drugs_growth = []
    four_drugs_growth = []
    
    for i in range(len(growth_final_norm)):
        for j in range(len(growth_final_norm[i])):
            plt.subplot(2,4,unq_drugs[i])
            if unq_drugs[i] == 1:
                plt.plot([0, 2, 4, 6], growth_final_norm[i][j], color='red', alpha=0.5, linewidth=2.0)
            if unq_drugs[i] == 2:
                plt.plot([0, 2, 4, 6], growth_final_norm[i][j], color='black', alpha=0.30, linewidth=0.30)
                two_drugs_growth.append(growth_final_norm[i][j][3])
            if unq_drugs[i] == 3:
                plt.plot([0, 2, 4, 6], growth_final_norm[i][j], color='black', alpha=0.10, linewidth=0.10)
                three_drugs_growth.append(growth_final_norm[i][j][3])
            if unq_drugs[i] == 4:
                plt.plot([0, 2, 4, 6], growth_final_norm[i][j], color='black', alpha=0.05, linewidth=0.05)
                four_drugs_growth.append(growth_final_norm[i][j][3])
    
    
    plt.subplot(2,4,2)
    for i in range(len(growth_final_norm[0])):
        plt.plot([0, 2, 4, 6], growth_final_norm[0][i], color='red', alpha=0.5, linewidth=2.0)
        plt.plot([0, 2, 4, 6], growth_final_norm[156][i], color='red', alpha=0.5, linewidth=2.0)
        plt.plot([0, 2, 4, 6], growth_final_norm[312][i], color='red', alpha=0.5, linewidth=2.0)
        plt.plot([0, 2, 4, 6], growth_final_norm[468][i], color='red', alpha=0.5, linewidth=2.0)
        plt.plot([0, 2, 4, 6], growth_final_norm[624][i], color='red', alpha=0.5, linewidth=2.0)
        
    plt.subplot(2,4,3)
    for i in range(len(growth_final_norm[0])):
        plt.plot([0, 2, 4, 6], growth_final_norm[0][i], color='red', alpha=0.5, linewidth=2.0)
        plt.plot([0, 2, 4, 6], growth_final_norm[156][i], color='red', alpha=0.5, linewidth=2.0)
        plt.plot([0, 2, 4, 6], growth_final_norm[312][i], color='red', alpha=0.5, linewidth=2.0)
        plt.plot([0, 2, 4, 6], growth_final_norm[468][i], color='red', alpha=0.5, linewidth=2.0)
        plt.plot([0, 2, 4, 6], growth_final_norm[624][i], color='red', alpha=0.5, linewidth=2.0)
        
    plt.subplot(2,4,4)
    for i in range(len(growth_final_norm[0])):
        plt.plot([0, 2, 4, 6], growth_final_norm[0][i], color='red', alpha=0.5, linewidth=2.0)
        plt.plot([0, 2, 4, 6], growth_final_norm[156][i], color='red', alpha=0.5, linewidth=2.0)
        plt.plot([0, 2, 4, 6], growth_final_norm[312][i], color='red', alpha=0.5, linewidth=2.0)
        plt.plot([0, 2, 4, 6], growth_final_norm[468][i], color='red', alpha=0.5, linewidth=2.0)
        plt.plot([0, 2, 4, 6], growth_final_norm[624][i], color='red', alpha=0.5, linewidth=2.0)
    
  
    two_drugs_totRes = []
    three_drugs_totRes = []
    four_drugs_totRes = []
    
    for i in range(len(growth_final_norm)):
        for j in range(len(growth_final_norm[i])):
            plt.subplot(2,4,unq_drugs[i]+4)
            if unq_drugs[i] == 1:
                plt.plot([0, 2, 4, 6, 8], totRes_final[i][j], color='red', alpha=0.5, linewidth=2.0)
            if unq_drugs[i] == 2:
                plt.plot([0, 2, 4, 6, 8], totRes_final[i][j], color='black', alpha=0.30, linewidth=0.30)
                two_drugs_totRes.append(totRes_final[i][j][4])
            if unq_drugs[i] == 3:
                plt.plot([0, 2, 4, 6, 8], totRes_final[i][j], color='black', alpha=0.10, linewidth=0.10)
                three_drugs_totRes.append(totRes_final[i][j][4])
            if unq_drugs[i] == 4:
                plt.plot([0, 2, 4, 6, 8], totRes_final[i][j], color='black', alpha=0.05, linewidth=0.05)
                four_drugs_totRes.append(totRes_final[i][j][4])

    plt.subplot(2,4,6)
    for i in range(len(growth_final_norm[0])):
        plt.plot([0, 2, 4, 6, 8], totRes_final[0][i], color='red', alpha=0.5, linewidth=2.0)
        plt.plot([0, 2, 4, 6, 8], totRes_final[156][i], color='red', alpha=0.5, linewidth=2.0)
        plt.plot([0, 2, 4, 6, 8], totRes_final[312][i], color='red', alpha=0.5, linewidth=2.0)
        plt.plot([0, 2, 4, 6, 8], totRes_final[468][i], color='red', alpha=0.5, linewidth=2.0)
        plt.plot([0, 2, 4, 6, 8], totRes_final[624][i], color='red', alpha=0.5, linewidth=2.0)
        
    plt.subplot(2,4,7)
    for i in range(len(growth_final_norm[0])):
        plt.plot([0, 2, 4, 6, 8], totRes_final[0][i], color='red', alpha=0.5, linewidth=2.0)
        plt.plot([0, 2, 4, 6, 8], totRes_final[156][i], color='red', alpha=0.5, linewidth=2.0)
        plt.plot([0, 2, 4, 6, 8], totRes_final[312][i], color='red', alpha=0.5, linewidth=2.0)
        plt.plot([0, 2, 4, 6, 8], totRes_final[468][i], color='red', alpha=0.5, linewidth=2.0)
        plt.plot([0, 2, 4, 6, 8], totRes_final[624][i], color='red', alpha=0.5, linewidth=2.0)
        
    plt.subplot(2,4,8)
    for i in range(len(growth_final_norm[0])):
        plt.plot([0, 2, 4, 6, 8], totRes_final[0][i], color='red', alpha=0.5, linewidth=2.0)
        plt.plot([0, 2, 4, 6, 8], totRes_final[156][i], color='red', alpha=0.5, linewidth=2.0)
        plt.plot([0, 2, 4, 6, 8], totRes_final[312][i], color='red', alpha=0.5, linewidth=2.0)
        plt.plot([0, 2, 4, 6, 8], totRes_final[468][i], color='red', alpha=0.5, linewidth=2.0)
        plt.plot([0, 2, 4, 6, 8], totRes_final[624][i], color='red', alpha=0.5, linewidth=2.0)
    
    for i in range(8): 
        plt.subplot(2,4,i+1)
        plt.ylim([-2, 30])

    maxVal = -100
    minVal = 100
    for i in range(len(growth_final_norm)):
        for j in range(len(growth_final_norm[i])):
            if growth_final_norm[i][j][3] < minVal:
                minVal = growth_final_norm[i][j][3]
            if growth_final_norm[i][j][3] > maxVal:
                maxVal = growth_final_norm[i][j][3]

    print(minVal)
    print(maxVal)

    minVal = 100
    maxVal = -100
    for i in range(len(growth_final_norm)):
        for j in range(len(growth_final_norm[i])):
            if totRes_final[i][j][4] < minVal:
                minVal = totRes_final[i][j][4]
            if totRes_final[i][j][4] > maxVal:
                maxVal = totRes_final[i][j][4]
    
    print(minVal)
    print(maxVal)
    
    
    
    plt.figure(figsize=(16.0, 16.0)) # in inches!
    
    plt.subplot(2,3,1)
    plt.hist(two_drugs_growth, bins=40, density=True, color='darkblue')
    plt.xlim(xmin=-2, xmax = 30)
    plt.xlabel('Resistance to dosing drug')
    plt.ylabel('Frequency')
    plt.title('Two drug outcomes')
    
    plt.subplot(2,3,2)
    plt.hist(three_drugs_growth, bins=40, density=True, color='darkblue')
    plt.xlim(xmin=-2, xmax = 30)
    plt.xlabel('Resistance to dosing drug')
    plt.ylabel('Frequency')
    plt.title('Three drug outcomes')
    
    plt.subplot(2,3,3)
    plt.hist(four_drugs_growth, bins=40, density=True, color='darkblue')
    plt.xlim(xmin=-2, xmax = 30)
    plt.xlabel('Resistance to dosing drug')
    plt.ylabel('Frequency')
    plt.title('Four drug outcomes')
    
    plt.subplot(2,3,4)
    plt.hist(two_drugs_totRes, bins=40, density=True, color='darkblue')
    plt.xlim(xmin=-2, xmax = 30)
    plt.xlabel('Total resistance acquired')
    plt.ylabel('Frequency')
    plt.title('Two drug outcomes')
    
    plt.subplot(2,3,5)
    plt.hist(three_drugs_totRes, bins=40, density=True, color='darkblue')
    plt.xlim(xmin=-2, xmax = 30)
    plt.xlabel('Total resistance acquired')
    plt.ylabel('Frequency')
    plt.title('Three drug outcomes')
    
    plt.subplot(2,3,6)
    plt.hist(four_drugs_totRes, bins=40, density=True, color='darkblue')
    plt.xlim(xmin=-2, xmax = 30)
    plt.xlabel('Total resistance acquired')
    plt.ylabel('Frequency')
    plt.title('Four drug outcomes')
    
    
    # Working on this section for only 2 drug cases
    stuff = []
    for i in range(5):
        stuff.append([])
        for j in range(5):
            stuff[i].append([])
            stuff[i][j].append([])
            if i != j and i < j:
                stuff[i][j][0].append(((0,0,0,0), 1000))
                stuff[i][j][0].append(((0,0,0,0), -1000))
            
    stuff2 = []
    for i in range(5):
        stuff2.append([])
        for j in range(5):
            stuff2[i].append([])
            stuff2[i][j].append([])
            if i != j and i < j:
                stuff2[i][j][0].append(((0,0,0,0), 1000))
                stuff2[i][j][0].append(((0,0,0,0), -1000))
    
    
    
    
    notherCount = 0
    
    for dose_seq in list(itertools.product(range(5), repeat=4)):
        dose_seq_set = set(dose_seq)
        if len(dose_seq_set) == 2:
            netValG = 0
            netTotRes = 0

            for j in range(len(growth_final_norm[notherCount])):
                netValG += growth_final_norm[notherCount][j][3]
                netTotRes += totRes_final[notherCount][j][4]
                    
            netValG = netValG / len(growth_final_norm[notherCount])
            netTotRes = netTotRes / len(growth_final_norm[notherCount])

            if netValG < stuff[list(dose_seq_set)[0]][list(dose_seq_set)[1]][0][0][1]: 
                stuff[list(dose_seq_set)[0]][list(dose_seq_set)[1]][0][0] = (dose_seq, netValG)
                
            if netValG > stuff[list(dose_seq_set)[0]][list(dose_seq_set)[1]][0][1][1]: 
                stuff[list(dose_seq_set)[0]][list(dose_seq_set)[1]][0][1] = (dose_seq, netValG)
                
            if netTotRes < stuff2[list(dose_seq_set)[0]][list(dose_seq_set)[1]][0][0][1]: 
                stuff2[list(dose_seq_set)[0]][list(dose_seq_set)[1]][0][0] = (dose_seq, netTotRes)
                
            if netTotRes > stuff2[list(dose_seq_set)[0]][list(dose_seq_set)[1]][0][1][1]: 
                stuff2[list(dose_seq_set)[0]][list(dose_seq_set)[1]][0][1] = (dose_seq, netTotRes)

                
            
        notherCount += 1
          
    
    mins_netGrow = [6.682, 3.06625, 8.65325, 3.327, 2.80575, 8.7665, 3.56225, 3.624, 1.25875, 2.908]
    maxs_netGrow = [8.1385, 9.95875, 15.838, 8.47075, 10.0920, 16.86275, 11.893, 14.08625, 3.4285, 12.28025]
    
    mins_netTot = [13.62575, 10.5570, 13.23075, 10.8155, 10.44875, 13.47175, 10.22125, 10.05375, 7.1525, 9.82625]
    maxs_netTot = [15.00125, 12.408, 14.94175, 11.4225, 12.202, 14.43575, 11.216, 11.84275, 8.62325, 11.3835]
    
    drug_pairs = ["0-1", "0-2", "0-3", "0-4", "1-2", "1-3", "1-4", "2-3", "2-4", "3-4"]
    
    width = 0.4
    
    indices = np.arange(10)
    
    plt.figure()
    plt.subplot(1,2,1)
    plt.bar(indices, mins_netGrow, width, color='darkblue', label='Best')
    plt.bar(indices + width, maxs_netGrow, width, color='red', label='Worst')
    plt.xticks(indices + width / 2, drug_pairs)
    plt.legend(loc='best')
    
    plt.subplot(1,2,2)
    plt.bar(indices, mins_netTot, width, color='darkblue', label='Best')
    plt.bar(indices + width, maxs_netTot, width, color='red', label='Worst')
    plt.xticks(indices + width / 2, drug_pairs)
    plt.legend(loc='best')

    
    plt.show()

if __name__ == '__main__':
    main()
        
        