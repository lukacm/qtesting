import numpy as np
from scipy.sparse import csr_matrix

def build_reversible_fa_matrix(triplets_str):
    """
    Builds a reversible permutation matrix for an FA.
    Ordering: |input, ancilla, state>
    Logic: Output_Row = Matrix * Input_Column
    """
    # 1. Parse the triplets
    lines = [line.split() for line in triplets_str.strip().split('\n') if line.strip()]
    
    inputs_list = sorted(list(set(l[0] for l in lines)))
    states_list = sorted(list(set(l[1] for l in lines) | set(l[2] for l in lines)))
    
    i_map = {val: idx for idx, val in enumerate(inputs_list)}
    s_map = {val: idx for idx, val in enumerate(states_list)}
    
    num_i = len(inputs_list)
    num_s = len(states_list)
    num_a = num_s  # Ancilla dimension matches state dimension
    
    # Transition function delta(input, state) -> next_state
    delta = {}
    for inp, curr, nxt in lines:
        delta[(i_map[inp], s_map[curr])] = s_map[nxt]

    # 2. Build the Matrix
    # Total dimension D = num_i * num_a * num_s
    dim = num_i * num_a * num_s
    rows = []
    cols = []
    data = []

    # Iterate through all possible basis states in the vector
    for i in range(num_i):
        for a in range(num_a):
            for s in range(num_s):
                # Column Index (Original State)
                # index = i*(Na*Ns) + a*(Ns) + s
                col_idx = (i * num_a * num_s) + (a * num_s) + s
                
                # Compute Next State logic
                # If transition not defined, we treat next_state as 0 (or stay put)
                nxt_idx = delta.get((i, s), 0)
                
                # Reversible step: a_new = (a_old + nxt_idx) mod num_s
                a_new = (a + nxt_idx) % num_a
                
                # Row Index (Target State)
                row_idx = (i * num_a * num_s) + (a_new * num_s) + s
                
                rows.append(row_idx)
                cols.append(col_idx)
                data.append(1.0)

    # Construct Sparse Matrix
    # v_out = M * v_in
    matrix = csr_matrix((data, (rows, cols)), shape=(dim, dim))
    
    return matrix, inputs_list, states_list

def get_basis_vector(i_val, a_val, s_val, inputs, states):
    """Helper to create a column vector representing |i, a, s>"""
    num_i = len(inputs)
    num_s = len(states)
    num_a = num_s
    
    i_idx = inputs.index(i_val)
    a_idx = states.index(a_val)
    s_idx = states.index(s_val)
    
    idx = (i_idx * num_a * num_s) + (a_idx * num_s) + s_idx
    v = np.zeros(num_i * num_a * num_s)
    v[idx] = 1
    return v.reshape(-1, 1)

def decode_vector(v, inputs, states):
    """Helper to turn a column vector back into readable labels"""
    idx = np.argmax(v)
    num_s = len(states)
    num_a = num_s
    
    # Reverse the index formula
    s_idx = idx % num_s
    a_idx = (idx // num_s) % num_a
    i_idx = idx // (num_a * num_s)
    
    return f"|In: {inputs[i_idx]}, Ancilla: {states[a_idx]}, State: {states[s_idx]}>"

## --- Example ---
## Non-reversible FA: both inputs lead to state B eventually
#fa_input = """
#0 A B
#0 B B
#0 C C
#1 A A
#1 B B
#1 C A
#"""
#
#matrix, inputs, states = build_reversible_fa_matrix(fa_input)
#
#print(f"Matrix Dimension: {matrix.shape[0]}x{matrix.shape[0]}")
#
## Test a transition: Input '0', Current State 'A', Ancilla 'A' (index 0)
## We expect the 'next state' to be added to the ancilla
#v_in = get_basis_vector('0', 'A', 'A', inputs, states)
#v_out = matrix @ v_in
#
#print("\nApplying Matrix to Column Vector:")
#print(f"Input:  {decode_vector(v_in, inputs, states)}")
#print(f"Output: {decode_vector(v_out, inputs, states)}")
#
## Check reversibility
#is_reversible = np.allclose((matrix.T @ matrix).toarray(), np.eye(matrix.shape[0]))
#print(f"\nIs the matrix reversible? {is_reversible}")
#print(matrix.todense())
