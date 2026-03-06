import numpy as np
import matplotlib.pyplot as plt
import qiskit
import math
#import aer 
#from qiskit.providers.aer.library import save_statevector
import matplotlib.pyplot as plt
from qiskit_aer import AerSimulator
from qiskit import  QuantumCircuit
from qiskit.compiler import transpile
from qiskit.transpiler import generate_preset_pass_manager
# Generate a pass manager without providing a backend
from qiskit.primitives import StatevectorEstimator
#from qiskit.primitives import StatevectorSampler as Sampler
from qiskit_ibm_runtime import SamplerV2 as Sampler
from qiskit_ibm_runtime import EstimatorV2 as Estimator
from qiskit_ibm_runtime.fake_provider import FakeManilaV2
#from qiskit_ibm_runtime import SamplerV2 as QiskitRuntimeService
from qiskit_ibm_runtime import Session, SamplerV2 as Sampler, QiskitRuntimeService
from qiskit.quantum_info import SparsePauliOp
from qiskit.circuit.library import UnitaryGate
from qiskit.quantum_info import Statevector
from qiskit.circuit.library import grover_operator, MCMTGate, ZGate
from qiskit.circuit.library import MCMTVChain
from qiskit.visualization import plot_distribution
import binembedP
 
fa_triplets = """
0 1 2
0 2 2
0 3 1
1 1 2
1 2 3
1 3 1
"""
#fa_triplets = """
#0 1 1
#0 2 2
#1 1 2
#1 2 2
#"""


#Qubits for embedding the FA
#Computation can be obtained from binembed.py for automation
embedding_ancillae = 1 

#Design the FA matrix
matrix, meta = binembedP.build_reversible_fa_binary(fa_triplets, embedding_ancillae)
unitary_fa_qubits= meta['total_bits']
input_qubits = meta['n_i']
state_qubits = meta['n_s']
unitary = matrix.todense()
gate = UnitaryGate(unitary)
 
#print(gate)
#Input seuqnece length
input_seq_length =3 

#Search requirements
search_addons = 1

def insertQgate(circuit, gate, meta):
#Construct the quantum circuit to simulate the FA for the input sequence
    for j in range(len(meta['states'])):
        for i in range(input_seq_length):
            qubits = list()
            qubits.append(i)
            for e in range(meta['n_a']):
                #qubits.append(input_seq_length+input_seq_length*meta['n_a']*j+e+i)
                qubits.insert(0,input_seq_length+input_seq_length*meta['n_a']*j+e+i)
            for s in range(meta['n_s']):
                #qubits.append(input_seq_length+len(meta['states'])*embedding_ancillae*input_seq_length+j*meta['n_s']+s)
                qubits.insert(0,input_seq_length+len(meta['states'])*embedding_ancillae*input_seq_length+j*meta['n_s']+s)
            print("State: ",j,qubits)
            circuit.unitary(gate, qubits)
    #        circuit.append(gate, qubits)


def grover_oracle(marked_states, qubits, targets):
    """Build a Grover oracle for multiple marked states

    Here we assume all input marked states have the same number of bits

    Parameters:
        marked_states (str or list): Marked states of oracle

    Returns:
        QuantumCircuit: Quantum circuit representing Grover oracle
    """
    if not isinstance(marked_states, list):
        marked_states = [marked_states]
    # Compute the number of qubits in circuit
#    num_qubits = len(marked_states[0])
    num_qubits = qubits
    qc = QuantumCircuit(num_qubits)
    # Mark each target state in the input list
    for target in marked_states:
        # Flip target bit-string to match Qiskit bit-ordering
        rev_target = target[::-1]
        # Find the indices of all the '0' elements in bit-string
        zero_inds = [
            ind+targets[0]
            for ind in range(num_qubits)
            if rev_target.startswith("0", ind)
        ]
        # Add a multi-controlled Z-gate with pre- and post-applied X-gates (open-controls)
        # where the target bit-string has a '0' entry
        if zero_inds:
            qc.x(zero_inds)
        qc.compose(MCMTGate(ZGate(), len(targets) - 1, 1), targets, inplace=True)
        #qc.compose(MCMTGate(ZGate(), num_qubits - 1, 1), targets, inplace=True)
        if zero_inds:
            qc.x(zero_inds)
    return qc

print(input_seq_length)
print(len(meta['states'])*input_seq_length*embedding_ancillae)
print(len(meta['states'])*meta['n_s'])
total_qubits = input_seq_length+len(meta['states'])*input_seq_length*embedding_ancillae+len(meta['states'])*meta['n_s']+search_addons
circuit = QuantumCircuit(total_qubits,2)
print(f'Circuit on {total_qubits} qubits created')



#prepare the inputs
for i in range(input_seq_length):
    circuit.h(i)
#Test
#circuit.x(0)
#circuit.x(1)
#circuit.append(gateic, [0,1])

#prepare the states
states = meta['state_map']
for i in range(len(meta['states'])):
    state = meta['state_map'][str(i+1)]
    for j in range(len(state)):
        if (state[j] == '1'):
            #circuit.x(input_seq_length+(i*meta['n_s'])+j)
            circuit.x(input_seq_length+len(meta['states'])*embedding_ancillae*input_seq_length+j*meta['n_s']+j)
 
#Construct the quantum circuit to simulate the FA for the input sequence
insertQgate(circuit, gate, meta)


#Create the Quantum Grover Oracle
oracle = QuantumCircuit(total_qubits)
#print(meta['state_map'])
ccontrols = []
for j in range(len(meta['states'])*meta['n_s']):
    ccontrols.append(total_qubits - 1 - len(meta['states'])*meta['n_s'] +j)
print(ccontrols)
states = meta['state_map']
marked = []
for i in range(len(meta['states'])):
#for i in range(1):
    oracle.h(total_qubits-1)
    negs = ((" ".join(states[f'{i+1}']*3)).replace(" ",""))
    marked.append(negs)
    for j in range(len(negs)):
        if (negs[j] == '0'):
            oracle.x(total_qubits - 1 - len(meta['states'])*meta['n_s'] +j)
    oracle.compose(MCMTGate(ZGate(), len(ccontrols)-1, 1),ccontrols, inplace=True)#, inplace=True)
#    oracle.h(total_qubits-1)
#    oracle.mcx(ccontrols,total_qubits-1)
#    oracle.h(total_qubits-1)
    for j in range(len(negs)):
        if (negs[j] == '0'):
            oracle.x(total_qubits - 1 - len(meta['states'])*meta['n_s'] +j)
    #oracle.h(total_qubits-1)

oracle1 = grover_oracle(marked, total_qubits, ccontrols)
#oracle1.draw(output='mpl')

oracle.draw(output='mpl')
plt.show()
#ccontrols.append(total_qubits-1)
grover_op = grover_operator(oracle1, reflection_qubits=ccontrols, insert_barriers=True)
grover_op.draw(output='mpl')
plt.show()

optimal_num_iterations = math.floor(
    math.pi / (4 * math.asin(math.sqrt(len(meta['states']) / 2**len(meta['states']))))
    #math.pi / (4 * math.asin(math.sqrt(len(meta['states']) / 2**len(ccontrols))))
    #/ (4 * math.asin(math.sqrt(len(meta['states']) / 2**grover_op.num_qubits)))
)
optimal_num_iterations=100
circuit.compose(grover_op.power(optimal_num_iterations), inplace=True)

circuit = transpile(circuit, optimization_level=3)#, basis_gates=["u3", "cx"])
#print(dict(circuit.count_ops()))
print(grover_op.decompose().draw(output="text"))
#print(circuit.draw(output='text'))

#print(circuit.draw(output='text'))

# Measure all qubits
#circuit.measure(9,0)
#circuit.measure(10,1)
#circuit.measure_all()
#circuit.draw(output="text", style="iqp")





#Generate the observable string for the circuit
#iden = 'I'*(total_qubits-len(meta['states'])*meta['n_s']-1)+'Z'*(len(meta['states'])*meta['n_s'])+ 'I'
observables = []
for j in range(len(negs)):
    observables.append(SparsePauliOp('I'*(total_qubits-len(meta['states'])*meta['n_s']-1+j)+'Z'+'I'*((len(meta['states'])*meta['n_s'])-j)))
print(observables) 
#observable = SparsePauliOp(observable,[1,1,1,1,1,1])#, ["IIIZZ"])
#observable = SparsePauliOp(["IIIXX"],[0.5])#, ["IIIZZ"])
#observable = SparsePauliOp(["IIIII", "XXXXX", "YYYYY", "ZZZZZ"], coeffs=[1, 1, -1, 1])
#observables = [
#    [SparsePauliOp(["XX", "IY"], [0.5, 0.5])],
#    [SparsePauliOp("XX")],
#    [SparsePauliOp("IY")],
#]
# value(s) for the circuit parameter(s)
#parameter_values = [[0], [np.pi / 6], [np.pi / 2]]

#service = QiskitRuntimeService()

# Run the sampler job locally using AerSimulator.
# Session syntax is supported but ignored because local mode doesn't support sessions.
aer_sim = AerSimulator()
pm = generate_preset_pass_manager(backend=aer_sim, optimization_level=1)
isa_circuit = pm.run(circuit)
with Session(backend=aer_sim) as session:
    sampler = Sampler()
    sampler.options.default_shots = 10_000
    result = sampler.run([isa_circuit]).result()
    dist = result[0].data#.meas.get_counts()
    psi = Statevector(isa_circuit)
    print(psi)


 
# Run the sampler job locally using FakeManilaV2
#fake_manila = FakeManilaV2()
#pm = generate_preset_pass_manager(backend=fake_manila,optimization_level=1)
#isa_circuit = pm.run(circuit)
layout = isa_circuit.layout
observables = [
    [observable.apply_layout(layout) for observable in observable_set]
    for observable_set in observables
]

 

# You can use a fixed seed to get fixed results.
#options = {"simulator": {"seed_simulator": 42}}
#sampler = Sampler(mode=fake_manila, options=options)

#Run on hardware
#sampler.options.default_shots = 10_000
#result = sampler.run([isa_circuit]).result()
#dist = result[0].data.meas.get_counts()

#Use the statevector estimator
#estimator = StatevectorEstimator()
#result = Estimator.run([(circuit, observables)]).result()
#print(result[0].data)
#dist = result[0].data.meas.get_counts()
#result = job.result()

#print(result[0].data.evs)
#result = job.result()

# Error-bar information is also available, but the error is 0
# for this StatevectorEstimator.
#result.data.stds
 
# Pull out the array-based expectation value estimate data from the
# result and plot a trace for each observable.
#for idx, pauli in enumerate(observables):
#    plt.plot(result.data.evs[idx], label=pauli)
#plt.legend()

sv = Statevector(circuit)
probas = sv.probabilities()
shots = sv.sample_counts(100)
plot_distribution(shots)
#print(sv)
#print(shots)

s_idx = np.argsort(sv.probabilities())
#print(s_idx)
l = len(s_idx)
form = '{0:0'+str(total_qubits)+'b}'
print(form.format(s_idx[l-1]), probas[s_idx[l-1]])
print(form.format(s_idx[l-2]), probas[s_idx[l-2]])
print(form.format(s_idx[l-3]), probas[s_idx[l-3]])
print(form.format(s_idx[l-4]), probas[s_idx[l-4]])
print(form.format(s_idx[l-1]), probas[s_idx[l-5]])
print(form.format(s_idx[l-2]), probas[s_idx[l-6]])
print(form.format(s_idx[l-3]), probas[s_idx[l-7]])
print(form.format(s_idx[l-4]), probas[s_idx[l-8]])
print(form.format(s_idx[l-1]), probas[s_idx[l-9]])
print(form.format(s_idx[l-2]), probas[s_idx[l-10]])
print(form.format(s_idx[l-3]), probas[s_idx[l-11]])
print(form.format(s_idx[l-4]), probas[s_idx[l-12]])
##print(Statevector(circuit).probabilities(qargs=[0]) )
#print(Statevector(circuit).probabilities(qargs=[1]) )
#print(Statevector(circuit).probabilities(qargs=[2]) )
#print(Statevector(circuit).probabilities(qargs=[3]) )
#print(Statevector(circuit).probabilities(qargs=[4]) )
#print(Statevector(circuit).probabilities(qargs=[5]) )
#print(Statevector(circuit).probabilities(qargs=[6]) )
#print(Statevector(circuit).probabilities(qargs=[7]) )
#print(Statevector(circuit).probabilities(qargs=[8]) )

print(
    f"The result of the submitted job had {len(result)} PUB and has a value:\n {result}\n"
)
print(
    f"The associated PubResult of this job has the following data bins:\n {result[0].data}\n"
)
print(f"And this DataBin has attributes: {result[0].data.keys()}")
print(
    "Recall that this shape is due to our array of parameter binding sets having shape (100, 2) -- where 2 is the\n\
         number of parameters in the circuit -- combined with our array of observables having shape (3, 1). \n"
)
#print(
#    f"The expectation values measured from this PUB are: \n{result[0].data.evs}"
#)

#print(f" > Result class: {type(result)}")
#print(result)
#statevector = result.get_statevector()
