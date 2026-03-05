from qiskit import QuantumCircuit
from qiskit.circuit import Parameter
from qiskit.quantum_info import SparsePauliOp
import numpy as np
# Generate a pass manager without providing a backend
from qiskit.transpiler import generate_preset_pass_manager
from qiskit.primitives import StatevectorEstimator
from qiskit.circuit.library import UnitaryGate
from qiskit.circuit.library import grover_operator, MCMTGate, ZGate
import binembed
from oracle import grover_oracle as go
import math
 
 
fa_triplets = """
0 A B
0 B C
0 C B
1 A C
1 B A
1 C A
"""

marked_states = ["000000", "010101", "101010"]

goracle = go(marked_states)
goracle.draw(output="mpl", style="iqp")
grover_op = grover_operator(goracle, reflection_qubits=[0, 1, 2], insert_barriers=True)
grover_op = grover_operator(goracle)
grover_op.decompose().draw(output="mpl", style="iqp")

optimal_num_iterations = math.floor(
    math.pi
    / (4 * math.asin(math.sqrt(len(marked_states) / 2**grover_op.num_qubits)))
)


qc = QuantumCircuit(grover_op.num_qubits)
# Create even superposition of all basis states
qc.h(range(grover_op.num_qubits))
# Apply Grover operator the optimal number of times
qc.compose(grover_op.power(optimal_num_iterations), inplace=True)
# Measure all qubits
qc.measure_all()
qc.draw(output="mpl", style="iqp")


matrix, meta = binembed.build_binary_reversible_matrix(fa_triplets)
total_bits= meta['total_bits']
unitary = matrix.todense()
gate = UnitaryGate(unitary)

circuit = QuantumCircuit(total_bits)
#circuit = QuantumCircuit(5)
#circuit.h(0)
circuit.append(gate, [0, 1, 2, 3, 4])
#circuit.draw("mpl", style="iqp")



# circuit for which you want to obtain the expected value
#circuit = QuantumCircuit(5)
circuit.ry(Parameter("theta"), 0)
#circuit.h(0)
#circuit.cx(0, 1)
circuit.draw("mpl", style="iqp")

 
# observable(s) whose expected values you want to compute
 
observable = SparsePauliOp(["IIIII", "XXXXX", "YYYYY", "ZZZZZ"], coeffs=[1, 1, -1, 1])
 
# value(s) for the circuit parameter(s)
parameter_values = [[0]]
#parameter_values = [[0], [np.pi / 6], [np.pi / 2]]

pm = generate_preset_pass_manager(optimization_level=1)
isa_circuit = pm.run(circuit)
isa_observable = observable.apply_layout(isa_circuit.layout)


estimator = StatevectorEstimator()
job = estimator.run([(circuit, observable, parameter_values)])
result = job.result()
print(f" > Result class: {type(result)}")

print(f" > Expectation value: {result[0].data.evs}")
print(f" > Metadata: {result[0].metadata}")
