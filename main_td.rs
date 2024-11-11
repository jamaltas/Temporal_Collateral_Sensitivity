use std::time::Instant;
use itertools::{Itertools, MultiProduct};
use std::fs::File;
use std::error::Error;
use csv::Reader;
use std::collections::HashMap;
use ndarray::{Array, Array1, Array2, Array3, Axis};
use rand::{thread_rng, Rng};
use rayon::prelude::*;


pub fn main() {
    // time the program
    let now = Instant::now();

    // set up the state space
    let min_resistance = -5;
    let max_resistance = 25;
    let drugs: usize = 5;
    let replicates: usize = 4;
    let resistance_levels: Vec<i8> = (min_resistance..=max_resistance).collect();
    let mut state_vectors: Vec<Vec<i8>> = Vec::new();
    for k2 in resistance_levels.iter().copied().product_repeat(drugs) {
        state_vectors.push(k2);
    }

    let mut state_to_state_dictionary: HashMap<Vec<i8>, usize> = HashMap::new();

    for i in 0..state_vectors.len() {
        state_to_state_dictionary.insert(state_vectors[i].clone(), i);
    }

    // read in and discretize the biological data
    let cs_data_6to8 = read_csv("day8_rust.csv").unwrap();
    let discretized_cs_data_6to8 = discretize_data(cs_data_6to8);

    // create the transition matrix
    let mut transition_matrix_6to8: Array3<usize> = Array::zeros((state_vectors.len(), drugs, replicates));
    for s in 0..state_vectors.len(){
        for d in 0..drugs {
            for r in 0..replicates {
                let current_state = state_vectors[s].clone();
                let next_state = get_next_resistance_state(&current_state, &discretized_cs_data_6to8[d*4 + r], min_resistance, max_resistance);
                transition_matrix_6to8[[s,d,r]] = state_to_state_dictionary[&next_state];
            }
        }
    }

    // calculate the MDP optimal for timepoints 6-8
    let gamma = 0.99;
    let optimal_policy_6to8: Vec<usize> = calc_optimal_policy(drugs, replicates, &transition_matrix_6to8, &state_vectors, gamma);


    //
    // Done calculating the optimal for time period 6-8. Next step is to simulate the outcome from each starting state. 
    // This allows us to identify the most optimal states to enter the 6-8 time period from. 
    // We can use this information as a reward for the 4-6 time period.
    // As a result the optimal 4-6 treatment is not the one that minimizes the reward over 4-6, but instead the one that enters the final period from the best state.
    // 
    
    // Using the calculated optimal, calculate the reward for using that optimal from each initial state.
    let time_steps = 20;
    let future_rewards_6to8 = calc_future_rewards(&transition_matrix_6to8, &state_vectors, &optimal_policy_6to8, replicates, time_steps);

    // read in and discretize the biological data
    let cs_data_4to6 = read_csv("day6_rust.csv").unwrap();
    let discretized_cs_data_4to6 = discretize_data(cs_data_4to6);

    // create the transition matrix
    let mut transition_matrix_4to6: Array3<usize> = Array::zeros((state_vectors.len(), drugs, replicates));
    for s in 0..state_vectors.len(){
        for d in 0..drugs {
            for r in 0..replicates {
                let current_state = state_vectors[s].clone();
                let next_state = get_next_resistance_state(&current_state, &discretized_cs_data_4to6[d*4 + r], min_resistance, max_resistance);
                transition_matrix_4to6[[s,d,r]] = state_to_state_dictionary[&next_state];
            }
        }
    }
    
    // calculate the MDP optimal for timepoints 4-6 by using the data from timepoints 6-8.
    let optimal_policy_4to6: Vec<usize> = calc_optimal_policy_td(drugs, replicates, &transition_matrix_4to6, &state_vectors, &future_rewards_6to8, gamma);

    //
    // Done calculating the optimal for time period 4-6. Next step is to simulate the outcome from each starting state. 
    // Here we simulate the 4-6 time period, and with the final state we put that into the future rewards vector to determine
    // what the final reward at the end of timepoint 8 will be.
    // 

    // Using the calculated optimal for 4-6, calculate the reward for using that optimal for each initial state
    // This is different from the previous step because the optimal is entirely based on entering the final spot in the best position.
    let future_rewards_4to6 = calc_future_rewards_projection(&transition_matrix_4to6, &state_vectors, &optimal_policy_4to6, &future_rewards_6to8, replicates, time_steps);

    // read in and discretize the biological data
    let cs_data_2to4 = read_csv("day4_rust.csv").unwrap();
    let discretized_cs_data_2to4 = discretize_data(cs_data_2to4);

    // create the transition matrix
    let mut transition_matrix_2to4: Array3<usize> = Array::zeros((state_vectors.len(), drugs, replicates));
    for s in 0..state_vectors.len(){
        for d in 0..drugs {
            for r in 0..replicates {
                let current_state = state_vectors[s].clone();
                let next_state = get_next_resistance_state(&current_state, &discretized_cs_data_2to4[d*4 + r], min_resistance, max_resistance);
                transition_matrix_2to4[[s,d,r]] = state_to_state_dictionary[&next_state];
            }
        }
    }

    // calculate the MDP optimal for timepoints 4-6 by using the data from timepoints 6-8.
    let optimal_policy_2to4: Vec<usize> = calc_optimal_policy_td(drugs, replicates, &transition_matrix_2to4, &state_vectors, &future_rewards_4to6, gamma);

    //
    // Done calculating the optimal for time period 2-4. Next step is to simulate the outcome from each starting state. 
    // Here we simulate the 2-4 time period, and with the final state we put that into the future rewards vector to determine
    // what the final reward at the end of timepoint 8 will be.
    // 

    // Using the calculated optimal for 4-6, calculate the reward for using that optimal for each initial state
    // This is different from the previous step because the optimal is entirely based on entering the final spot in the best position.
    let future_rewards_2to4 = calc_future_rewards_projection(&transition_matrix_2to4, &state_vectors, &optimal_policy_2to4, &future_rewards_4to6, replicates, time_steps);

    // read in and discretize the biological data
    let cs_data_0to2 = read_csv("day2_rust.csv").unwrap();
    let discretized_cs_data_0to2 = discretize_data(cs_data_0to2);

    // create the transition matrix
    let mut transition_matrix_0to2: Array3<usize> = Array::zeros((state_vectors.len(), drugs, replicates));
    for s in 0..state_vectors.len(){
        for d in 0..drugs {
            for r in 0..replicates {
                let current_state = state_vectors[s].clone();
                let next_state = get_next_resistance_state(&current_state, &discretized_cs_data_0to2[d*4 + r], min_resistance, max_resistance);
                transition_matrix_0to2[[s,d,r]] = state_to_state_dictionary[&next_state];
            }
        }
    }

    // calculate the MDP optimal for timepoints 0-2 by using the data from timepoints 2-4.
    let optimal_policy_0to2: Vec<usize> = calc_optimal_policy_td(drugs, replicates, &transition_matrix_0to2, &state_vectors, &future_rewards_2to4, gamma);


    let mut file = std::fs::File::create("opt_policy(0-2).pickle").unwrap();
    serde_pickle::to_writer(&mut file, &optimal_policy_0to2, serde_pickle::SerOptions::new()).unwrap();

    let mut file = std::fs::File::create("opt_policy(2-4).pickle").unwrap();
    serde_pickle::to_writer(&mut file, &optimal_policy_2to4, serde_pickle::SerOptions::new()).unwrap();

    let mut file = std::fs::File::create("opt_policy(4-6).pickle").unwrap();
    serde_pickle::to_writer(&mut file, &optimal_policy_4to6, serde_pickle::SerOptions::new()).unwrap();

    let mut file = std::fs::File::create("opt_policy(6-8).pickle").unwrap();
    serde_pickle::to_writer(&mut file, &optimal_policy_6to8, serde_pickle::SerOptions::new()).unwrap();

    let mut file2 = std::fs::File::create("transition_matrix1.pickle").unwrap();
    serde_pickle::to_writer(&mut file2, &transition_matrix_0to2, serde_pickle::SerOptions::new()).unwrap();

    let mut file2 = std::fs::File::create("transition_matrix2.pickle").unwrap();
    serde_pickle::to_writer(&mut file2, &transition_matrix_2to4, serde_pickle::SerOptions::new()).unwrap();

    let mut file2 = std::fs::File::create("transition_matrix3.pickle").unwrap();
    serde_pickle::to_writer(&mut file2, &transition_matrix_4to6, serde_pickle::SerOptions::new()).unwrap();

    let mut file2 = std::fs::File::create("transition_matrix4.pickle").unwrap();
    serde_pickle::to_writer(&mut file2, &transition_matrix_6to8, serde_pickle::SerOptions::new()).unwrap();

    println!("{}", now.elapsed().as_secs());
}



// Code for iterrools to work how I want
pub trait ProductRepeat: Iterator + Clone
    where Self::Item: Clone {
    fn product_repeat(self, repeat: usize) -> MultiProduct<Self> {
        std::iter::repeat(self)
        .take(repeat)
        .multi_cartesian_product()
    }
}

impl<T: Iterator + Clone> ProductRepeat for T
    where T::Item: Clone {}


fn read_csv(file_path: &str) -> Result<Vec<Vec<f64>>, Box<dyn Error>> {
    let file = File::open(file_path)?;
    let mut csv_reader = Reader::from_reader(file);

    let mut cs_data: Vec<Vec<f64>> = Vec::new();

    for result in csv_reader.records() {
        let record = result?;
        let row: Vec<f64> = record.iter()
            .map(|field| field.parse::<f64>().unwrap_or(0.0))
            .collect();
        cs_data.push(row);
    }

    Ok(cs_data)
}

fn discretize_data(cs_data: Vec<Vec<f64>>) -> Vec<Vec<i8>> {
    let mut disc_cs_data: Vec<Vec<i8>> = Vec::new();

    for i in 0..cs_data.len() {
        let mut row: Vec<i8> = Vec::new();
        for j in 0..cs_data[i].len() {
            let value = cs_data[i][j];
            let discretized_value = if value < -2.0 { 
                -2 
            } else if value < -0.25 {
                -1
            } else if value < 0.25 {
                0
            } else if value < 2.0 {
                1
            } else {
                2
            };
            row.push(discretized_value);
        }
        disc_cs_data.push(row);
    }

    disc_cs_data
}

fn get_next_resistance_state(a: &[i8], b: &[i8], r_min: i8, r_max: i8) -> Vec<i8> {
    let next_state: Vec<i8> = a.iter()
        .zip(b.iter())
        .map(|(&x, &y)| x + y)
        .map(|sum| sum.min(r_max).max(r_min))
        .collect();

    next_state
}

fn calc_optimal_policy(drugs: usize, replicates: usize, transition_matrix: &Array3<usize>, state_vectors: &Vec<Vec<i8>>, gamma: f64) -> Vec<usize> {

    // Turn Vec<Vec<i8>> into 2d Array for vectorized math.
    let flattened_states = state_vectors.into_iter().flatten().map(|&x| f64::from(x)).collect();
    let array_states: Array2<f64> = Array::from_shape_vec((state_vectors.len(), state_vectors[0].len()), flattened_states).unwrap();

    let theta: f64 = 0.001;
    let max_iterations = 7500;
    let mut v_new: Array1<f64> = Array1::from(vec![0.0; state_vectors.len()]);

    println!("Computing the optimal policy for gamma: {}", gamma);
    let mut count = 0;
    loop {
        let mut delta: f64 = 0.0;

        let v_old = v_new.clone();
        let reward_mat = one_step_lookahead(drugs, replicates, &transition_matrix, &array_states, &v_new, gamma);
        v_new = reward_mat.map_axis(Axis(1), |view| *view.into_iter().min_by(|a, b| a.partial_cmp(b).unwrap()).unwrap());

        let difference = v_old.clone()-v_new.clone();
        let abs_values: Vec<f64> = difference.map(|&x| x.abs()).to_vec();
        let max_val = abs_values.iter().copied().max_by(|a, b| a.partial_cmp(b).unwrap()).unwrap();

        if delta < max_val {
            delta = max_val;
        } else {
            delta = delta;
        } 

        count += 1;
        if count % 100 == 0 {
            println!("{}", count);
            println!("{}", delta);
        }
        if delta < theta || count > max_iterations {
            break;
        }

    }

    let reward_mat = one_step_lookahead(drugs, replicates, &transition_matrix, &array_states, &v_new, gamma);
    let policy: Vec<usize> = reward_mat.map_axis(Axis(1), |view| view.iter().position(|&x| x == *view.iter().min_by(|a, b| a.partial_cmp(b).unwrap()).unwrap()).unwrap()).to_vec();
    policy

}

fn one_step_lookahead(drugs: usize, replicates: usize, transition_matrix: &Array3<usize>, array_states: &Array2<f64>, v_new: &Array1<f64>, gamma: f64) -> Array2<f64> {

    // Get the size of the 2D array
    let shape = array_states.shape();
    
    let mut reward_mat: Array2<f64> = Array::zeros((shape[0], drugs));

    for r in 0..replicates {
        
        let replicate_matrix: Array2<usize> = transition_matrix.index_axis(Axis(2), r).into_owned();

        let reward_replicate = replicate_matrix
            .iter()
            .map(|&idx| v_new[idx])
            .collect::<Array1<f64>>()
            .into_shape(replicate_matrix.dim())
            .unwrap();

        reward_mat = reward_mat + ((1.0 as f64/(replicates as f64)) * (array_states + gamma * reward_replicate));
        
    }

    reward_mat
}

fn calc_optimal_policy_td(drugs: usize, replicates: usize, transition_matrix: &Array3<usize>, state_vectors: &Vec<Vec<i8>>, future_rewards: &Vec<f64>, gamma: f64) -> Vec<usize> {

    // Turn Vec<Vec<i8>> into 2d Array for vectorized math.
    let flattened_states = state_vectors.into_iter().flatten().map(|&x| f64::from(x)).collect();
    let array_states: Array2<f64> = Array::from_shape_vec((state_vectors.len(), state_vectors[0].len()), flattened_states).unwrap();

    // Get a 2dArray of state_reward for vectorized math
    let mut future_rewards_matrix = Array::zeros((future_rewards.len(), drugs));
    for s in 0..future_rewards.len() {
        for d in 0..drugs {
            future_rewards_matrix[[s,d]] = future_rewards[s];
        }
    }

    let theta: f64 = 0.001;
    let max_iterations = 7500;
    let mut v_new: Array1<f64> = Array1::from(vec![0.0; future_rewards.len()]);

    println!("Computing the optimal policy for gamma: {}", gamma);
    let mut count = 0;
    loop {
        let mut delta: f64 = 0.0;

        let v_old = v_new.clone();
        let reward_mat = one_step_lookahead_td(drugs, replicates, &transition_matrix, &v_new, gamma, &future_rewards_matrix, &array_states);
        v_new = reward_mat.map_axis(Axis(1), |view| *view.into_iter().min_by(|a, b| a.partial_cmp(b).unwrap()).unwrap());

        let difference = v_old.clone()-v_new.clone();
        let abs_values: Vec<f64> = difference.map(|&x| x.abs()).to_vec();
        let max_val = abs_values.iter().copied().max_by(|a, b| a.partial_cmp(b).unwrap()).unwrap();

        if delta < max_val {
            delta = max_val;
        } else {
            delta = delta;
        } 

        count += 1;
        if count % 100 == 0 {
            println!("{}", count);
            println!("{}", delta);
        }
        if delta < theta || count > max_iterations {
            break;
        }

    }

    let reward_mat = one_step_lookahead_td(drugs, replicates, &transition_matrix, &v_new, gamma, &future_rewards_matrix, &array_states);
    let policy: Vec<usize> = reward_mat.map_axis(Axis(1), |view| view.iter().position(|&x| x == *view.iter().min_by(|a, b| a.partial_cmp(b).unwrap()).unwrap()).unwrap()).to_vec();
    policy

}

fn one_step_lookahead_td(drugs: usize, replicates: usize, transition_matrix: &Array3<usize>, v_new: &Array1<f64>, gamma: f64, future_rewards_matrix: &Array2<f64>, array_states: &Array2<f64>) -> Array2<f64> {

    // Get the size of the 2D array
    let shape = future_rewards_matrix.shape();
    
    let mut reward_mat: Array2<f64> = Array::zeros((shape[0], drugs));

    for r in 0..replicates {
        
        let replicate_matrix: Array2<usize> = transition_matrix.index_axis(Axis(2), r).into_owned();

        let reward_replicate = replicate_matrix
            .iter()
            .map(|&idx| v_new[idx])
            .collect::<Array1<f64>>()
            .into_shape(replicate_matrix.dim())
            .unwrap();

        reward_mat = reward_mat + ((1.0 as f64/(replicates as f64)) * (array_states + future_rewards_matrix + gamma * reward_replicate));
        
    }

    reward_mat
}


fn simulate_policy(transition_matrix: &Array3<usize>, optimal_policy: &Vec<usize>, state_vectors: &Vec<Vec<i8>>, time_steps: usize, initial_state: usize, replicates: usize) -> (Vec<i8>, usize) {
    
    let mut current_state = initial_state.clone();

    let mut reward: Vec<i8> = Vec::with_capacity(time_steps as usize);
    let mut state: Vec<usize> = Vec::with_capacity(time_steps as usize);

    for _ in 0..time_steps {

        let next_action = optimal_policy[current_state];

        state.push(current_state);
        reward.push(state_vectors[current_state][next_action]);

        let random_mutant = thread_rng().gen_range(0..replicates);
        let next_state = transition_matrix[[current_state,next_action,random_mutant]];
        current_state = next_state;
    }

    (reward, current_state)
}

fn calc_future_rewards(transition_matrix: &Array3<usize>, state_vectors: &Vec<Vec<i8>>, optimal_policy: &Vec<usize>, replicates: usize, time_steps: usize) -> Vec<f64> {

    // number of times to simulate from an initital state to calculate reward
    let traj_to_sim = 100;
    let mut reward: Vec<f64> = vec![0.0; state_vectors.len()];
    reward.par_iter_mut().enumerate().for_each(|(initial_state, temp_reward)| {
        for _ in 0..traj_to_sim {
            let (optimal_rewards, _curr_state) = simulate_policy(&transition_matrix, &optimal_policy, &state_vectors, time_steps, initial_state, replicates);
            *temp_reward += *optimal_rewards.last().unwrap() as f64;
        }
        *temp_reward /= traj_to_sim as f64;
    });
    
    reward
}

fn calc_future_rewards_projection(transition_matrix: &Array3<usize>, state_vectors: &Vec<Vec<i8>>, optimal_policy: &Vec<usize>, future_rewards: &Vec<f64>, replicates: usize, time_steps: usize) -> Vec<f64> {

    // number of times to simulate from an initital state to calculate reward
    let traj_to_sim = 100;
    let mut reward: Vec<f64> = vec![0.0; state_vectors.len()];
    reward.par_iter_mut().enumerate().for_each(|(initial_state, temp_reward)| {
        for _ in 0..traj_to_sim {
            let (_optimal_rewards, curr_state) = simulate_policy(&transition_matrix, &optimal_policy, &state_vectors, time_steps, initial_state, replicates);
            *temp_reward += future_rewards[curr_state];
        }
        *temp_reward /= traj_to_sim as f64;
    });
    
    reward
}


fn running_mean(x: &[f32], N: usize) -> Vec<f32> {
    let cumsum: Vec<f32> = x
        .iter()
        .scan(0.0, |state, &value| {
            *state += value;
            Some(*state)
        })
        .collect();

    let mut result = Vec::with_capacity(x.len() - N + 1);

    for i in N..x.len() {
        let mean = (cumsum[i] - cumsum[i - N]) / N as f32;
        result.push(mean);
    }

    result
}