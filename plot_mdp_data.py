import pickle as pkl
import numpy as np
import copy
import random
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from scipy.integrate import simpson
from matplotlib.colors import TwoSlopeNorm

import itertools
import csv

def main():

    with open('opt_policy(0-2)td.pickle', 'rb') as handle:
        opt_policy1_td = pkl.load(handle)
    with open('opt_policy(2-4)td.pickle', 'rb') as handle:
        opt_policy2_td = pkl.load(handle)
    with open('opt_policy(4-6)td.pickle', 'rb') as handle:
        opt_policy3_td = pkl.load(handle)
    with open('opt_policy(6-8)td.pickle', 'rb') as handle:
        opt_policy4_td = pkl.load(handle)

    with open('opt_policy(0-2).pickle', 'rb') as handle:
        opt_policyd2 = pkl.load(handle)
    with open('opt_policy(0-4).pickle', 'rb') as handle:
        opt_policyd4 = pkl.load(handle)
    with open('opt_policy(0-6).pickle', 'rb') as handle:
        opt_policyd6 = pkl.load(handle)
    with open('opt_policy(0-8).pickle', 'rb') as handle:
        opt_policyd8 = pkl.load(handle)


    min_resistance = -5
    max_resistance = 25


    initial_state = 4771525


    evoSteps=20
    opt_rewards_TD = np.zeros((1000,evoSteps*4))
    opt_rewards_d2 = np.zeros((1000,evoSteps*4))
    opt_rewards_d4 = np.zeros((1000,evoSteps*4))
    opt_rewards_d6 = np.zeros((1000,evoSteps*4))
    opt_rewards_d8 = np.zeros((1000,evoSteps*4))

    with open('opt_rewards_TD.pickle', 'rb') as handle:
        opt_rewards_TD = pkl.load(handle)

    with open('opt_rewards_d2.pickle', 'rb') as handle:
        opt_rewards_d2 = pkl.load(handle)
    
    with open('opt_rewards_d4.pickle', 'rb') as handle:
        opt_rewards_d4 = pkl.load(handle)
    
    with open('opt_rewards_d6.pickle', 'rb') as handle:
        opt_rewards_d6 = pkl.load(handle)

    with open('opt_rewards_d8.pickle', 'rb') as handle:
        opt_rewards_d8 = pkl.load(handle)




    opt_rewards_d2_mean = running_mean(np.insert(np.mean(opt_rewards_d2,axis=0),0,0),4)
    opt_rewards_d4_mean = running_mean(np.insert(np.mean(opt_rewards_d4,axis=0),0,0),4)
    opt_rewards_d6_mean = running_mean(np.insert(np.mean(opt_rewards_d6,axis=0),0,0),4)
    opt_rewards_d8_mean = running_mean(np.insert(np.mean(opt_rewards_d8,axis=0),0,0),4)
    opt_rewards_TD_mean = running_mean(np.insert(np.mean(opt_rewards_TD,axis=0),0,0),4)


    xval_running = np.linspace(4,evoSteps*4,evoSteps*4-3)
    x = np.linspace(0,evoSteps*4,evoSteps*4+1)
    plt.figure(figsize=(8,6))
    plt.subplot(1,1,1)
    #plt.plot(xval_running, running_mean(np.mean(opt_rewards_TD, axis=0),4),'r', linestyle='solid', label="day2 opt")
    plt.axvline(x=evoSteps, color='k', linestyle='--', linewidth=3)
    plt.axvline(x=evoSteps*2, color='k', linestyle='--', linewidth=3)
    plt.axvline(x=evoSteps*3, color='k', linestyle='--', linewidth=3)
    plt.plot(xval_running, opt_rewards_d2_mean, color='#00ffe3', label="Day2")
    plt.plot(xval_running, opt_rewards_d4_mean, color='#03cffc', label="Day4")
    plt.plot(xval_running, opt_rewards_d6_mean, color='#0689f9', label="Day6")
    plt.plot(xval_running, opt_rewards_d8_mean, color='#050ffa', label="Day8")
    plt.plot(xval_running, opt_rewards_TD_mean, color='#000000',  label="TD opt")

    # Calculate the standard error of the mean
    sem_d2 = running_mean(np.insert(np.std(opt_rewards_d2,axis=0),0,0),4) / np.sqrt(opt_rewards_d2.shape[0]/10)
    sem_d4 = running_mean(np.insert(np.std(opt_rewards_d4,axis=0),0,0),4) / np.sqrt(opt_rewards_d4.shape[0]/10)
    sem_d6 = running_mean(np.insert(np.std(opt_rewards_d6,axis=0),0,0),4) / np.sqrt(opt_rewards_d6.shape[0]/10)
    sem_d8 = running_mean(np.insert(np.std(opt_rewards_d8,axis=0),0,0),4) / np.sqrt(opt_rewards_d8.shape[0]/10)
    sem_TD = running_mean(np.insert(np.std(opt_rewards_TD,axis=0),0,0),4) / np.sqrt(opt_rewards_TD.shape[0]/10)

    # Plot the standard error of the mean as shaded areas
    plt.fill_between(xval_running, opt_rewards_d2_mean - sem_d2,
                    opt_rewards_d2_mean + sem_d2, color='#00ffe3', alpha=0.3)

    plt.fill_between(xval_running, opt_rewards_d4_mean - sem_d4,
                    opt_rewards_d4_mean + sem_d4, color='#03cffc', alpha=0.3)

    plt.fill_between(xval_running, opt_rewards_d6_mean - sem_d6,
                    opt_rewards_d6_mean + sem_d6, color='#0689f9', alpha=0.3)

    plt.fill_between(xval_running, opt_rewards_d8_mean - sem_d8,
                    opt_rewards_d8_mean + sem_d8, color='#050ffa', alpha=0.3)

    plt.fill_between(xval_running, opt_rewards_TD_mean - sem_TD,
                    opt_rewards_TD_mean + sem_TD, color='#000000', alpha=0.3)
    
    plt.ylim([-6,26])
    plt.legend(loc='upper right', fontsize=16)
    plt.ylabel('Resistance to Applied Drug', fontsize=20)
    plt.xlabel('Time (Evolutionary Steps)', fontsize=20)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)

    plt.tight_layout()
    plt.savefig('reward_over_time.png', dpi=300)

    plt.figure(figsize=(7,6))
    area_2 = simpson(opt_rewards_d2_mean, dx=1)
    area_4 = simpson(opt_rewards_d4_mean, dx=1)
    area_6 = simpson(opt_rewards_d6_mean, dx=1)
    area_8 = simpson(opt_rewards_d8_mean, dx=1)
    area_TD = simpson(opt_rewards_TD_mean, dx=1)

    # Sample data for the 5 bars
    data = [area_2, area_4, area_6, area_8, area_TD]

    # Specify the colors for the 5 bars
    colors = ['#00ffe3', '#03cffc', '#0689f9', '#050ffa', '#000000']

    # Position of the bars on the x-axis
    x_pos = np.arange(len(data))

    # Create the bar plot
    plt.bar(x_pos, data, color=colors, edgecolor='black')

    # Customize the plot
    plt.xlabel('Treatment schedule', fontsize=20)
    plt.ylabel('Total resistance to drug', fontsize=20)
    plt.xticks(x_pos, ['Day2', 'Day4', 'Day6', 'Day8', 'TD opt'], fontsize=16)
    plt.yticks(fontsize=16)

    plt.tight_layout()
    plt.savefig('area_under.png', dpi=300)

    #
    # Compare optimal solutions.
    #

    plt.figure(figsize=(15,4))

    frequencies_d2 = np.bincount(opt_policyd2, minlength=5)/len(opt_policyd2)
    frequencies_d4 = np.bincount(opt_policyd4, minlength=5)/len(opt_policyd2)
    frequencies_d6 = np.bincount(opt_policyd6, minlength=5)/len(opt_policyd2)
    frequencies_d8 = np.bincount(opt_policyd8, minlength=5)/len(opt_policyd2)

    frequencies_td2 = np.bincount(opt_policy1_td, minlength=5)/len(opt_policyd2)
    frequencies_td4 = np.bincount(opt_policy2_td, minlength=5)/len(opt_policyd2)
    frequencies_td6 = np.bincount(opt_policy3_td, minlength=5)/len(opt_policyd2)
    frequencies_td8 = np.bincount(opt_policy4_td, minlength=5)/len(opt_policyd2)

    barWidth = 0.50
    x_labels = np.arange(5)*1.25
    x2_labels = x_labels + barWidth
    drugs=["CIP", "DAP", "DOX", "CRO", "LZD"]

    plt.subplot(1,4,1)

    plt.bar(x_labels, frequencies_d2, width=barWidth, color='#00ffe3', edgecolor='black', label="Day2")
    plt.bar(x2_labels, frequencies_td2, width=barWidth, color='#000000', edgecolor='black', label="TD2")
    plt.xticks([r*1.25 + barWidth/2 for r in range(len(drugs))], drugs, fontsize=12)
    plt.yticks(fontsize=12)
    plt.legend(loc='upper left', fontsize=12)
    plt.ylabel('Use Frequency', fontsize=16)
    plt.xlabel('Antibiotic', fontsize=16)

    plt.subplot(1,4,2)

    plt.bar(x_labels, frequencies_d4, width=barWidth, color='#03cffc', edgecolor='black', label="Day4")
    plt.bar(x2_labels, frequencies_td4, width=barWidth, color='#000000', edgecolor='black', label="TD4")
    plt.xticks([r*1.25  + barWidth/2 for r in range(len(drugs))], drugs, fontsize=12)
    plt.yticks(fontsize=12)
    plt.legend(loc='best', fontsize=12)
    plt.ylabel('Use Frequency', fontsize=16)
    plt.xlabel('Antibiotic', fontsize=16)

    plt.subplot(1,4,3)

    plt.bar(x_labels, frequencies_d6, width=barWidth, color='#0689f9', edgecolor='black', label="Day6")
    plt.bar(x2_labels, frequencies_td6, width=barWidth, color='#000000', edgecolor='black', label="TD6")
    plt.xticks([r*1.25  + barWidth/2 for r in range(len(drugs))], drugs, fontsize=12)
    plt.yticks(fontsize=12)
    plt.legend(loc='best', fontsize=12)
    plt.ylabel('Use Frequency', fontsize=16)
    plt.xlabel('Antibiotic', fontsize=16)

    plt.subplot(1,4,4)

    plt.bar(x_labels, frequencies_d8, width=barWidth, color='#050ffa', edgecolor='black', label="Day8")
    plt.bar(x2_labels, frequencies_td8, width=barWidth, color='#000000', edgecolor='black', label="TD8")
    plt.xticks([r*1.25  + barWidth/2 for r in range(len(drugs))], drugs, fontsize=12)
    plt.yticks(fontsize=12)
    plt.ylabel('Use Frequency', fontsize=16)
    plt.xlabel('Antibiotic', fontsize=16)


    plt.legend(loc='best', fontsize=12)

    plt.tight_layout()
    plt.savefig('drug_distributions.png', dpi=300)


    plt.figure(figsize=(18,4))

    with open('state_TD.pickle', 'rb') as handle:
        TD_state = pkl.load(handle)

    with open('state_d2.pickle', 'rb') as handle:
        d2_state = pkl.load(handle)
    
    with open('state_d4.pickle', 'br') as handle:
        d4_state = pkl.load(handle)
    
    with open('state_d6.pickle', 'rb') as handle:
        d6_state = pkl.load(handle)

    with open('state_d8.pickle', 'rb') as handle:
        d8_state = pkl.load(handle)

    TD_state = TD_state.T
    d2_state = d2_state.T
    d4_state = d4_state.T
    d6_state = d6_state.T
    d8_state = d8_state.T

    day2 = np.zeros((5,5))
    day4 = np.zeros((5,5))
    day6 = np.zeros((5,5))
    day8 = np.zeros((5,5))

    day2[:,0] = d2_state[:,0]
    day2[:,1] = d4_state[:,0]
    day2[:,2] = d6_state[:,0]
    day2[:,3] = d8_state[:,0]
    day2[:,4] = TD_state[:,0]

    day4[:,0] = d2_state[:,1]
    day4[:,1] = d4_state[:,1]
    day4[:,2] = d6_state[:,1]
    day4[:,3] = d8_state[:,1]
    day4[:,4] = TD_state[:,1]

    day6[:,0] = d2_state[:,2]
    day6[:,1] = d4_state[:,2]
    day6[:,2] = d6_state[:,2]
    day6[:,3] = d8_state[:,2]
    day6[:,4] = TD_state[:,2]

    day8[:,0] = d2_state[:,3]
    day8[:,1] = d4_state[:,3]
    day8[:,2] = d6_state[:,3]
    day8[:,3] = d8_state[:,3]
    day8[:,4] = TD_state[:,3]

    norm = TwoSlopeNorm(vmin=-5, vcenter=0, vmax=25)

    plt.subplot(1,4,1)
    plt.imshow(day2, cmap='bwr', norm=norm)
    plt.yticks(np.arange(5), drugs, fontsize=12)
    plt.xticks(np.arange(5), ['Day2', 'Day4', 'Day6', 'Day8', 'TD opt'], fontsize=12)
    #plt.title('Day 2', fontsize=16)
    plt.hlines(y=np.arange(0, 5)+0.5, xmin=np.full(5, 0)-0.5, xmax=np.full(5, 5)-0.5, color="black")
    plt.vlines(x=np.arange(0, 5)+0.5, ymin=np.full(5, 0)-0.5, ymax=np.full(5, 5)-0.5, color="black")
    plt.xlabel('Treatment Schedule Applied', fontsize=16)
    plt.ylabel('Resistance to Drug', fontsize=16)
    #plt.colorbar()

    plt.subplot(1,4,2)
    plt.imshow(day4, cmap='bwr', norm=norm)
    plt.yticks(np.arange(5), drugs, fontsize=12)
    plt.xticks(np.arange(5), ['Day2', 'Day4', 'Day6', 'Day8', 'TD opt'], fontsize=12)
    plt.hlines(y=np.arange(0, 5)+0.5, xmin=np.full(5, 0)-0.5, xmax=np.full(5, 5)-0.5, color="black")
    plt.vlines(x=np.arange(0, 5)+0.5, ymin=np.full(5, 0)-0.5, ymax=np.full(5, 5)-0.5, color="black")
    plt.xlabel('Treatment Schedule Applied', fontsize=16)
    plt.ylabel('Resistance to Drug', fontsize=16)
    #plt.title('Day 4', fontsize=16)
    #plt.colorbar()

    plt.subplot(1,4,3)
    plt.imshow(day6, cmap='bwr', norm=norm)
    plt.yticks(np.arange(5), drugs, fontsize=12)
    plt.xticks(np.arange(5), ['Day2', 'Day4', 'Day6', 'Day8', 'TD opt'], fontsize=12)
    plt.hlines(y=np.arange(0, 5)+0.5, xmin=np.full(5, 0)-0.5, xmax=np.full(5, 5)-0.5, color="black")
    plt.vlines(x=np.arange(0, 5)+0.5, ymin=np.full(5, 0)-0.5, ymax=np.full(5, 5)-0.5, color="black")
    plt.xlabel('Treatment Schedule Applied', fontsize=16)
    plt.ylabel('Resistance to Drug', fontsize=16)
    #plt.title('Day 6', fontsize=16)
    #plt.colorbar()

    plt.subplot(1,4,4)
    plt.imshow(day8, cmap='bwr', norm=norm)
    plt.yticks(np.arange(5), drugs, fontsize=12)
    plt.xticks(np.arange(5), ['Day2', 'Day4', 'Day6', 'Day8', 'TD opt'], fontsize=12)
    plt.hlines(y=np.arange(0, 5)+0.5, xmin=np.full(5, 0)-0.5, xmax=np.full(5, 5)-0.5, color="black")
    plt.vlines(x=np.arange(0, 5)+0.5, ymin=np.full(5, 0)-0.5, ymax=np.full(5, 5)-0.5, color="black")
    plt.xlabel('Treatment Schedule Applied', fontsize=16)
    plt.ylabel('Resistance to Drug', fontsize=16)
    #plt.title('Day 8', fontsize=16)
    #plt.colorbar()

    data2 = [
        day2,
        day4,
        day6,
        day8
    ]

    # Save to CSV
    with open('opt_rewards_means.csv', 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerows(data2)

    plt.tight_layout()
    plt.savefig('resistance_profiles.png', dpi=300)

    plt.show()


def simulate_policy(env_statePrime, policy, state_vectors, evoSteps, initial_state):

    currentState = copy.deepcopy(initial_state)
    reward = np.zeros((evoSteps,1))
    state = np.zeros((evoSteps,1))

    for tStep in range(evoSteps):

        nextAction = policy[currentState]

        state[tStep] = currentState
        reward[tStep] = state_vectors[currentState,nextAction]

        nextState = random.choice(env_statePrime[currentState,nextAction,:])
        currentState = copy.deepcopy(nextState)

    return np.transpose(reward), currentState
    #return np.transpose(state), np.transpose(reward), currentState


def simulate_random_policy(env_statePrime, state_vectors, evoSteps, action_count, initial_state):

    currentState = copy.deepcopy(initial_state)
    reward = np.zeros((evoSteps,1))
    state = np.zeros((evoSteps,1))

    for tStep in range(evoSteps):
        nextAction = random.randrange(action_count)

        state[tStep] = currentState
        reward[tStep] = state_vectors[currentState,nextAction]

        nextState = random.choice(env_statePrime[currentState,nextAction,:])
        currentState = copy.deepcopy(nextState)

    return np.transpose(state), np.transpose(reward)



def one_step_lookahead(env_statePrime, state_vectors, states, V, gamma):
    """
    Helper function to calculate the value for all actions in a given state.

    Args:
        state: The state to consider (int)
        V: The value to use as an estimator, Vector of length # of states.

    Returns:
        A vector of length # of actions containing the expected value of each action.
    """

    A = np.zeros((len(state_vectors),5,1))
    state_vectors_test = np.zeros((len(state_vectors),5,1))
    state_vectors_test[:,:,0] = state_vectors
    
    for mut in range(4):
        A += 0.25 * (state_vectors_test + gamma*V[env_statePrime[:,:,mut]])

    return A
         

def calcOptimalPolicy(env_statePrime, state_vectors, gamma):

    V = np.zeros((len(state_vectors),1))
    theta = 0.001
    MaxIterations = 7500
    count = 0
    print(gamma)
    #delta = 1 #stopping condition

    while True:
        delta = 0

        v = copy.deepcopy(V)
        A = one_step_lookahead(env_statePrime, state_vectors, len(state_vectors), V, gamma)
        V = np.min(A, axis=1)

        delta = max(delta, np.max(np.abs(v-V)))

        count += 1
        if (count % 100 == 0):
            print(count)
            print(delta)

        if (delta < theta) or (count > MaxIterations):
            break

    policy = np.zeros((len(state_vectors),1))

    A = one_step_lookahead(env_statePrime, state_vectors, len(state_vectors), V, gamma)
    best_action = np.argmin(A, axis=1)
    policy = best_action

    return policy, V


def get_state_vectors(start,end,size):
    state_vector = np.arange(start, end + 1)
    state_vectors = np.array(list(itertools.product(state_vector, repeat=size)))
    return state_vectors

def get_next_resistance_state(a, b, r_min, r_max):
    return tuple([min(r_max, max(r_min, x + y)) for x, y in zip(a, b)])


def running_mean(x, N):
    cumsum = np.cumsum(x)
    return (cumsum[N:] - cumsum[:-N]) / float(N)

def discretize_data(data):
    for i in range(5):
        for j in range(20):
            if data[i,j] < -2:
                data[i,j] = -2
            elif data[i,j] < -0.25:
                data[i,j] = -1
            elif data[i,j] < 0.25:
                data[i,j] = 0
            elif data[i,j] < 2:
                data[i,j] = 1
            elif data[i,j] > 2:
                data[i,j] = 2
    return data

if __name__ == '__main__':
    main()
    #cProfile.run('main()')