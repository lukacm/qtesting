import numpy as np
from scipy.sparse import csr_matrix
import math

def build_reversible_fa_binary(triplets_str, num_ancilla_bits):
    """
    Builds a reversible binary permutation matrix.
    Order: |Input Bits, Ancilla Bits, State Bits>
    Matrix maps column vectors (inputs) to row vectors (outputs).
    """
    # 1. Parse Input
    lines = [line.split() for line in triplets_str.strip().split('\n') if line.strip()]
    unique_inputs = sorted(list(set(l[0] for l in lines)))
    unique_states = sorted(list(set(l[1] for l in lines) | set(l[2] for l in lines)))
    
    # Calculate required bits for Input and State
    n_i = math.ceil(math.log2(len(unique_inputs))) if len(unique_inputs) > 1 else 1
    n_s = math.ceil(math.log2(len(unique_states))) if len(unique_states) > 1 else 1
    n_a = num_ancilla_bits

    # Create Binary Encodings
    input_to_idx = {val: i for i, val in enumerate(unique_inputs)}
    state_to_idx = {val: i for i, val in enumerate(unique_states)}

    # Helper to format as binary string
    def to_bin_str(idx, bits): return f"{idx:0{bits}b}"

    input_encodings = {val: to_bin_str(i, n_i) for val, i in input_to_idx.items()}
    state_encodings = {val: to_bin_str(i, n_s) for val, i in state_to_idx.items()}
    
    # Mapping labels to integers
    i_map = {val: i for i, val in enumerate(unique_inputs)}
    s_map = {val: i for i, val in enumerate(unique_states)}
    
    # Transition function delta(input, state) -> next_state
    delta = {}
    for inp, curr, nxt in lines:
        delta[(i_map[inp], s_map[curr])] = s_map[nxt]

    # Total bits N = n_i + n_a + n_s
    total_bits = n_i + n_a + n_s
    dim = 2**total_bits
    
    rows = []
    cols = []
    data = []

    # Ancilla mask to ensure XOR stays within the parameter bit-width
    a_mask = (1 << n_a) - 1
    
    skip = 2**n_a
    ros = list('1' for va in range(0,2**total_bits))
    # Iterate through all possible basis states
    for i_val in range(2**n_i):
        for a_val in range(2**n_a):
            for s_val in range(len(unique_states)):
            #for s_val in range(2**n_s):
                # Logic: Find the next state
                # If the transition isn't defined, we use 0 (standard for FA padding)
                next_state_idx = delta.get((i_val, s_val), 0)
#                print(next_state_idx)
                if (a_val == 0):
                    next_ros_idx = (i_val << (n_a + n_s)) | (0 << n_s) | next_state_idx
                    if (ros[next_ros_idx] == '1'):
                       row_idx = next_ros_idx
                       ros[next_ros_idx] = '0'
                    else:
                       for i in range(1,2**n_a):
                          if (ros[(i_val << (n_a + n_s)) | (i << n_s) | next_state_idx] == '1'):
                             row_idx = (i_val << (n_a + n_s)) | (i << n_s) | next_state_idx
                             ros[(i_val << (n_a + n_s)) | (i << n_s) | next_state_idx] = '0'
                             break
                else:  
#                     print(ros)
                     for i in range(0,2**total_bits):
                       if (ros[i] == '1'):
                          ros[i] = '0'
                          row_idx = i
                          break
                col_idx = (i_val << (n_a + n_s)) | (a_val << n_s) | s_val
                 
                rows.append(row_idx)
                cols.append(col_idx)
                data.append(1.0)
    #Fill in the rows and cols completely unused
    for i_val in range(2**n_i):
        for a_val in range(2**n_a):
            for s_val in range(len(unique_states), 2**n_s):
                col_idx = (i_val << (n_a + n_s)) | (a_val << n_s) | s_val
                for i in range(0,2**total_bits):
                    if (ros[i] == '1'):
                       ros[i] = '0'
                       row_idx = i
                       break
  
                rows.append(row_idx)
                cols.append(col_idx)
                data.append(1.0)
    # Construct the sparse matrix M where v_out = M * v_in
    matrix = csr_matrix((data, (rows, cols)), shape=(dim, dim))
    
    meta = {'n_i': n_i, 'n_a': n_a, 'n_s': n_s, 'inputs': unique_inputs, 'states': unique_states, 'total_bits' : total_bits, 'state_map': state_encodings}
    return matrix, meta

def create_vector(i_lab, a_bits, s_lab, meta):
    """Helper to create a column vector representing the binary state |i, a, s>"""
    i_idx = meta['inputs'].index(i_lab)
    s_idx = meta['states'].index(s_lab)
    # a_bits is passed as an integer (e.g., 0)
    
    idx = (i_idx << (meta['n_a'] + meta['n_s'])) | (a_bits << meta['n_s']) | s_idx
    v = np.zeros(2**(meta['n_i'] + meta['n_a'] + meta['n_s']))
    v[idx] = 1.0
    return v.reshape(-1, 1)

def decode_vector(v, meta):
    """Decodes column vector index back to labels"""
    idx = np.argmax(v)
    
    s_mask = (1 << meta['n_s']) - 1
    a_mask = (1 << meta['n_a']) - 1
    
    s_idx = idx & s_mask
    a_val = (idx >> meta['n_s']) & a_mask
    i_idx = (idx >> (meta['n_a'] + meta['n_s']))
    
    i_lab = meta['inputs'][i_idx] if i_idx < len(meta['inputs']) else f"bin({i_idx})"
    s_lab = meta['states'][s_idx] if s_idx < len(meta['states']) else f"bin({s_idx})"
    
    return f"|Input: {i_lab}, Ancilla (val): {a_val}, State: {s_lab}>"

# --- Example Execution ---
fa_data = """
0 1 2
0 2 3
0 3 2
1 1 3
1 2 1
1 3 1
"""

old = """
0 A B
0 B C
0 C B
1 A C
1 B A
1 C A
"""
#
## PARAMETER: Number of ancilla bits
N_ANCILLA = 1 
#
matrix, metadata = build_reversible_fa_binary(old, N_ANCILLA)
unitary = matrix.todense()
#print(unitary)
#print(f"Total Bits: {metadata['n_i'] + metadata['n_a'] + metadata['n_s']}")
#print(f"Matrix Shape: {matrix.shape}")
## Test a transition: Input '0', Current State 'A', Ancilla starts at 0
#v_in = create_vector('0', 0, 'A', metadata)
#v_out = matrix @ v_in
#
#print("\nResult of v_out = M * v_in:")
#print(f"v_in:  {decode_vector(v_in, metadata)}")
#print(f"v_out: {decode_vector(v_out, metadata)}")
#
## Verify Reversibility: M.T @ M should be Identity
#is_rev = np.allclose((matrix.T @ matrix).toarray(), np.eye(matrix.shape[0]))
#print(f"\nMatrix is Reversible: {is_rev}")
