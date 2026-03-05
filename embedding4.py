import numpy as np
from scipy.sparse import csr_matrix

def build_reversible_fa_matrix(triplet_text):
    """
    Builds a reversible permutation matrix from FA transitions.
    Order of registers: |Input, Ancilla, State>
    """
    # 1. Parse Input
    lines = [line.split() for line in triplet_text.strip().split('\n') if line.strip()]
    
    inputs_found = sorted(list(set(l[0] for l in lines)))
    states_found = sorted(list(set(l[1] for l in lines) | set(l[2] for l in lines)))
    
    input_to_idx = {val: i for i, val in enumerate(inputs_found)}
    state_to_idx = {val: i for i, val in enumerate(states_found)}
    
    num_i = len(inputs_found)
    num_s = len(states_found)
    num_a = num_s  # Ancilla size must match state size to allow modular addition
    
    # 2. Build the transition map delta(i, s) -> s'
    # Default to s' = 0 if a transition is not specified
    delta = {}
    for inp, curr, nxt in lines:
        delta[(input_to_idx[inp], state_to_idx[curr])] = state_to_idx[nxt]

    # 3. Construct the Permutation Matrix
    # Total dimension D = |I| * |A| * |S|
    dim = num_i * num_a * num_s
    rows = []
    cols = []
    data = []

    for i in range(num_i):
        for a in range(num_a):
            for s in range(num_s):
                # Calculate the row index for |i, a, s>
                # Formula: i*(Na*Ns) + a*(Ns) + s
                row_idx = (i * num_a * num_s) + (a * num_s) + s
                
                # Get next state from FA
                next_state_idx = delta.get((i, s), 0) 
                
                # Reversible logic: Ancilla register tracks the transition result
                # a_new = (a_old + delta(i, s)) mod num_states
                a_new = (a + next_state_idx) % num_a
                
                # Calculate the column index for |i, a_new, s>
                col_idx = (i * num_a * num_s) + (a_new * num_s) + s
                
                rows.append(row_idx)
                cols.append(col_idx)
                data.append(1.0)

    # Create Sparse Matrix (more efficient for large state spaces)
    matrix = csr_matrix((data, (rows, cols)), shape=(dim, dim))
    
    return matrix, inputs_found, states_found

def is_permutation_matrix(matrix):
    """Checks if the matrix is a valid reversible permutation matrix."""
    # Check if every row and every column sums to 1
    row_sums = np.array(matrix.sum(axis=1)).flatten()
    col_sums = np.array(matrix.sum(axis=0)).flatten()
    return np.all(row_sums == 1) and np.all(col_sums == 1)

# --- Example Usage ---
# Format: Input CurrentState NextState
fa_triplets = """
0 S0 S1
0 S1 S2
0 S2 S1
1 S0 S0
1 S1 S0
1 S2 S1
"""

matrix, inputs, states = build_reversible_fa_matrix(fa_triplets)

print(f"FA parsed: {len(inputs)} inputs, {len(states)} states.")
print(f"Matrix Dimension: {matrix.shape[0]}x{matrix.shape[0]}")
print(f"Register Order: |Input, Ancilla, State>")
print(f"Is Reversible: {is_permutation_matrix(matrix)}")

# Print the mapping for the first few entries
num_s = len(states)
print("\nSample Mappings (|i, a, s> -> |i, a', s>):")
for i_idx in range(len(inputs)):
    for s_idx in range(num_s):
        # We look at ancilla = 0 to see where the transition is 'stored'
        row = (i_idx * num_s * num_s) + (0 * num_s) + s_idx
        col = matrix.indices[matrix.indptr[row]] # Find column where '1' is
        
        # Decode col index back to (i, a_new, s)
        new_a = (col // num_s) % num_s
        print(f"|Input:{inputs[i_idx]}, Ancilla:0, State:{states[s_idx]}>  -->  "
              f"|Input:{inputs[i_idx]}, Ancilla:{states[new_a]}, State:{states[s_idx]}>")
print(matrix.todense())
