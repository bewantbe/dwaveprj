# factoring program for QPU

import dwavebinarycsp as dbc

# Use a dimod test sampler that gives the BQM value for all values of its variables
from dimod import ExactSolver
exact_sampler = ExactSolver()

# Set an integer to factor
P = 11*13 # 25777 # 42

n_factor_bits = 4 # 3

# A binary representation of P ("{:06b}" formats for 6-bit binary)
bP = ("{:0%db}"%(2*n_factor_bits)).format(P)

csp = dbc.factories.multiplication_circuit(n_factor_bits)
# sorted(csp.variables.keys())

# Convert the CSP into BQM bqm
bqm = dbc.stitch(csp, min_classical_gap=0.1)
# Print a sample coefficient (one of the programable inputs to a D-Wave system)
print("# of qbits:", len(bqm.linear.values()))
print("# of couples:", len(bqm.quadratic.values()))

# To see helper functions, select Jupyter File Explorer View from the Online Learning page
#from helpers import draw
#draw.circuit_from(bqm)


# Our multiplication_circuit() creates these variables
p_vars = ['p'+str(i) for i in range(2*n_factor_bits)]

# Convert P from decimal to binary
fixed_variables = dict(zip(reversed(p_vars), bP))
fixed_variables = {var: int(x) for(var, x) in fixed_variables.items()}

# Fix product variables
for var, value in fixed_variables.items():
    bqm.fix_variable(var, value)
    
# Confirm that a P variable has been removed from the BQM, for example, "p0"
print("Variable p0 in BQM: ", 'p0' in bqm)
print("Variable a0 in BQM: ", 'a0' in bqm)

# Solve the BQM
#solution = exact_sampler.sample(bqm)
#list(solution.data())

# Setting Up a Solver
from dwave.system.samplers import DWaveSampler

#from helpers.solvers import default_solver
#my_solver, my_token = default_solver()
# Use a D-Wave system as the sampler
#sampler = DWaveSampler(solver=my_solver, token=my_token)
sampler = DWaveSampler()
_, target_edgelist, target_adjacency = sampler.structure

import dimod
#from helpers.embedding import embeddings

# Set a pre-calculated minor-embeding
#embedding = embeddings[sampler.solver.id]
#bqm_embedded = dimod.embed_bqm(bqm, embedding, target_adjacency, 3.0)

# Minor-Embedding
import minorminer

# Find an embedding
embedding = minorminer.find_embedding(bqm.quadratic, target_edgelist)
if bqm and not embedding:
    raise ValueError("no embedding found")

# Apply the embedding to the factoring problem to map it to the QPU
bqm_embedded = dimod.embed_bqm(bqm, embedding, target_adjacency, 3.0)


# Request num_reads samples
kwargs = {}
kwargs['num_reads'] = 1000
if 'answer_mode' in sampler.parameters:
    kwargs['answer_mode'] = 'histogram'
response = sampler.sample(bqm_embedded, **kwargs)

# Map back to the BQM's graph (nodes labeled "a0", "b0" etc,)
response = dimod.unembed_response(response, embedding, source_bqm=bqm)
print("\nThe solution in problem variables: \n",next(response.data(fields=['sample'])))

## Viewing the Solution
from collections import OrderedDict

#from helpers.convert import to_base_ten

def to_base_ten(qs):
    def getbits(kyref):
        li = [ky for ky in qs.keys() if ky[0]==kyref and ky[1]<'9']
        v = sum([(qs[ky] * 2**ii) for ii, ky in enumerate(sorted(li))])
        return v
    a = getbits('a')
    b = getbits('b')
    return (a,b)

# Function for converting the response to a dict of integer values
def response_to_dict(response):
    results_dict = OrderedDict()
    for sample, energy in response.data(['sample', 'energy']):
        # Convert A and B from binary to decimal
        a, b = to_base_ten(sample)
        # Aggregate results by unique A and B values (ignoring internal circuit variables)
        if (a, b) not in results_dict:
            results_dict[(a, b)] = energy
            
    return results_dict

# Function for converting the response to a dict of integer values
def response_to_dict_eng(response):
    results_dict = OrderedDict()
    for sample, energy, nocc in response.data(['sample', 'energy', 'num_occurrences']):
        # Convert A and B from binary to decimal
        a, b = to_base_ten(sample)
        # Aggregate results by unique A and B values (ignoring internal circuit variables)
        if (a, b, energy) not in results_dict:
            results_dict[(a, b, energy)] = nocc
            
    return results_dict

def response_to_engprob(response):
    results_dict = OrderedDict()
    for energy, nocc in response.data(['energy', 'num_occurrences']):
        if energy not in results_dict:
            results_dict[energy] = nocc
        else:
            results_dict[energy] += nocc
            
    return results_dict

# Convert the dimod.Response object to a dict and display it
results = response_to_dict(response)
draw.energy_of(results)

results = response_to_engprob(response)
draw.energy_of(results)

