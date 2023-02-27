# import custom libraries
import itertools
from convergence_ABM import *
from graph_generators import *
from modularity_util import *

# distribution parameters
param1, param2 = 0, 1

# initial values
s = 1000
V0 = np.random.normal(param1, param2, s) 
features = [AgentFeatures(V0[i], 1, 1) for i in range(s)]

# initialise models
config_param_c = {# "graph" : G,
               "features" : features,
               "selection" : random_selection, # select neighbours randomly
               "interaction" : comm_mean, # calculate the average at each interaction
               "next_move" : next_expo, # random expovariate
               "max_time" : 80, # total simulation time
               "log_interval" : 5, # logging interval
               "event_logger" : False, # log all events
               "time_logger" : True # log values at intervals
               }

results = {}
rounds = 100    

try:
    for i in range(1000):
        # set parameters
        k = np.random.randint(2,101)
        m = np.random.randint(2,16)
        beta = np.random.uniform(0.45,0.95)
        alpha = np.random.uniform(0.1,0.9)
        gamma = np.random.uniform(0.1,0.9)
        # generate graph
        G, node_comm = get_connected_community(s,k,m,beta,alpha,gamma)
        # update graph in configuration
        config_param_c["graph"] = G 
        # for community selection G
        config_param_c["selection"] = community_selection
        _, _, err_m, m_m = simulate_bulk(config_param_c, rounds) 
        # for random selection G
        config_param_c["selection"] = random_selection
        _, _, err_c, m_c = simulate_bulk(config_param_c, rounds) 
        # record results
        res = {"k": k, "m": m, "beta": beta, "alpha": alpha, "gamma": gamma, 
        "err_m_f": -np.log10(err_m[-1]), "m_m_f": m_m, 
        "err_c_f": -np.log10(err_c[-1]), "m_c_f": m_c,        
        } 
        # save results
        stats_net = get_struct_stats(G, node_comm)
        res.update(stats_net)
        results[i] = res
        # inform user
        print(str(i) + " done")

finally:
    # put in dataframe
    df = pd.DataFrame(results).T

    df.to_csv("results.csv")
