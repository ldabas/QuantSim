# # C ++ Tutorial

## quantum state

### Generation of quantum state
The following code generates a quantum state of <code> n </code> qubit.
The generated quantum state is initialized to \ f $ | 0 \ rangle ^ {\ otimes n} \ f $.
```cpp
#include <cppsim / state.hpp>

int main () {
// Generate 5-qubit state
unsigned int n = 5;
QuantumState state (n);
// Initialize to | 00000>
state.set_zero_state ();
return 0;
}
```
If there is not enough memory, the program will terminate.


### Initialization of quantum state
The generated quantum state can be initialized to a computational basis or to a random state.
```cpp
#include <cppsim / state.hpp>

int main () {
unsigned int n = 5;
QuantumState state (n);
state.set_zero_state ();
// Initialize to | 00101>
state.set_computational_basis (0b00101);
// Generate random initial state
state.set_Haar_random_state ();
// Generate random initial state by specifying seed
state.set_Haar_random_state (0);
return 0;
}
```


### Copy and load quantum state data
You can duplicate a quantum state or load data from another quantum state.
```cpp
#include <cppsim / state.hpp>

int main () {
unsigned int n = 5;
QuantumState state (n);
state.set_computational_basis (0b00101);

// copy and create new quantum state
auto second_state = state.copy ();

// Create a new quantum state and copy the existing state vector
QuantumState third_state (n);
third_state.load (& state);
return 0;
}
```


### Classic Register Operations
Quantum states have classical registers and can be read and written.
```cpp
#include <cppsim / state.hpp>

int main () {
unsigned int n = 5;
QuantumState state (n);
state.set_zero_state ();

// write register
int register_position = 3;
int register_value = 1;
state.set_classical_bit (register_position, register_value);

// read register
int obtained_value;
obtained_value = state.get_classical_bit (register_position);
return 0;
}
```

### Calculation on quantum state
The following processing is possible as a calculation that does not change the quantum state.
Calculations that change the quantum state are always performed via quantum gates and quantum circuits.
```cpp
#include <cppsim / state.hpp>

int main () {
unsigned int n = 5;
QuantumState state (n);
state.set_zero_state ();

// Calculate norm
double norm = state.get_squared_norm ();
// Calculate entropy when measured on Z basis
double entropy = state.get_entropy ();

// calculate index-th qubit with Z basis to get 0
unsigned int index = 3;
double zero_prob = state.get_zero_probability (index);

// Calculate marginal probabilities (The following is an example of the probability that 0,3-th qubit is measured as 0 and 1,2-th qubit is measured as 1)
std :: vector <unsigned int> value_list = {0,1,1,0,2};
double marginal_prob = state.get_marginal_probability (value_list);
return 0;
}
```

### Dot product of quantum states
The inner product can be calculated with the <code> inner_product </ code> function.
```cpp
#include <cppsim / state.hpp>

int main () {
unsigned int n = 5;
QuantumState state_ket (n);
state_ket.set_zero_state ();

QuantumState state_bra (n);
state_bra.set_Haar_random_state ();

std :: complex <double> value = state :: inner_product (& state_ket, & state_bra);
return 0;
}
```


### Acquisition of quantum state data
Get an array of length \ f $ 2 ^ n \ f $ representing the quantum states.
Note that creating a quantum state on a GPU, especially for large \ f $ n \ f $, is a very heavy operation.
```cpp
#include <cppsim / state.hpp>

int main () {
unsigned int n = 5;
QuantumState state (n);
state.set_zero_state ();

// For GNU C ++, get double _Complex array
// Get std :: complex <double> array for MSVC
const CTYPE * raw_data_c = state.data_c ();

// get std :: complex <double> array
const CPPCTYPE * raw_data_cpp = state.data_cpp ();
}
```

If you want to set the quantum state directly to the specified array, we recommend that you create the corresponding quantum gate and use it as an action of the quantum gate.


## Quantum gate

### Creation and action of quantum gate
The quantum gate implemented by default is created through the function of gate_factory, and the quantum state pointer is acted on as an argument. The quantum gate created by gate_factory is not released automatically, so the user must release it.

```cpp
#define _USE_MATH_DEFINES
#include <cmath>
#include <cppsim / state.hpp>
#include <cppsim / gate_factory.hpp>

int main () {
unsigned int n = 5;
QuantumState state (n);
state.set_zero_state ();

// X gate action
unsigned int index = 3;
auto x_gate = gate :: X (index);
x_gate-> update_quantum_state (& state);

// PI / 2 rotation at Y
double angle = M_PI / 2.0;
auto ry_gate = gate :: RY (index, angle);
ry_gate-> update_quantum_state (& state);

delete x_gate;
delete ry_gate;
return 0;
}
```

Gates defined in <code> gate </ code> namespace are as follows.
-single-qubit Pauli operation: Identity, X, Y, Z
-single-qubit Clifford operation: H, S, Sdag, T, Tdag, sqrtX, sqrtXdag, sqrtY, sqrtYdag
-two-qubit Clifford operation: CNOT, CZ, SWAP
-single-qubit Pauli rotation: RX, RY, RZ
-General Pauli operation: Pauli, PauliRotation
-IBMQ basis-gate: U1, U2, U3
-General gate: DenseMatrix
-Measurement: Measurement
-Noise: BitFlipNoise, DephasingNoise, IndepenedentXZNoise, DepolarizingNoise


### Synthesis of quantum gate
Quantum gates can be synthesized to create new quantum gates.
You must release the synthesized gate yourself.
```cpp
#define _USE_MATH_DEFINES
#include <cmath>
#include <cppsim / state.hpp>
#include <cppsim / gate_factory.hpp>
#include <cppsim / gate_merge.hpp>
#include <cppsim / gate_matrix.hpp>

int main () {
unsigned int n = 5;
QuantumState state (n);
state.set_zero_state ();

unsigned int index = 3;
auto x_gate = gate :: X (index);

double angle = M_PI / 2.0;
auto ry_gate = gate :: RY (index, angle);

// Create a gate that acts in the order of X, RY
auto x_and_ry_gate = gate :: merge (x_gate, ry_gate);

x_and_ry_gate-> update_quantum_state (& state);

delete x_gate;
delete ry_gate;
delete x_and_ry_gate;
return 0;
}
```

### Sum of quantum gate matrices
You can take the sum of the gate elements of a quantum gate.
(Do not use the sum when there is a control-qubit because the current operation is undefined.)
```cpp
#define _USE_MATH_DEFINES
#include <cmath>
#include <cppsim / state.hpp>
#include <cppsim / gate_factory.hpp>
#include <cppsim / gate_merge.hpp>
#include <cppsim / gate_matrix.hpp>

int main () {
auto gate00 = gate :: merge (gate :: P0 (0), gate ::P0 (1));
auto gate11 = gate :: merge (gate :: P1 (0), gate :: P1 (1));
// | 00> <00 | + | 11> <11 |
auto proj_00_or_11 = gate :: add (gate00, gate11);
std :: cout << proj_00_or_11 << std :: endl;

auto gate_ii_zz = gate :: add (gate :: Identity (0), gate :: merge (gate :: Z (0), gate :: Z (1)));
auto gate_ii_xx = gate :: add (gate :: Identity (0), gate :: merge (gate :: X (0), gate :: X (1)));
auto proj_00_plus_11 = gate :: merge (gate_ii_zz, gate_ii_xx);
// ((| 00> + | 11>) (<00 | + <11 |)) / 2 = (II + ZZ) (II + XX) / 4
proj_00_plus_11-> multiply_scalar (0.25);
std :: cout << proj_00_plus_11 << std :: endl;
return 0;
}
```

### Special quantum gates and general quantum gates
The basic quantum gate in cppsim is divided into the following two.
-Special gate: A function that has a dedicated speed-up function for the function of the gate.
-General gate: A gate that holds a matrix and acts by multiplying the matrix.

The former is faster than the latter because a dedicated function is created, but operations that change the action of the quantum gate, such as increasing the number of control qubits, cannot be performed later.
If you want to make these changes, you have to convert the special gates to general gates.

This can be achieved with <code> gate :: convert_to_matrix_gate </ code>.
Here is an example:

```cpp
#include <cppsim / state.hpp>
#include <cppsim / gate_factory.hpp>
#include <cppsim / gate_merge.hpp>
#include <cppsim / gate_matrix.hpp>

int main () {
unsigned int n = 5;
QuantumState state (n);
state.set_zero_state ();

unsigned int index = 3;
auto x_gate = gate :: X (index);

// Add control qubit so that it works only when 1st-qubit is 0
auto x_mat_gate = gate :: to_matrix_gate (x_gate);
unsigned int control_index = 1;
unsigned int control_with_value = 0;
x_mat_gate-> add_control_qubit (control_index, control_with_value);

x_mat_gate-> update_quantum_state (& state);

delete x_gate;
delete x_mat_gate;
return 0;
}
```

See the API documentation for a list of dedicated quantum gates.


### Get gate matrix of quantum gate
You can get the gate matrix of the generated quantum gate. Control qubits are not included in the gate matrix. In particular, be careful of gates that do not have a gate matrix (for example, an n-qubit Pauli rotation gate), which requires a very large amount of memory and time to acquire.
```cpp
#include <iostream>
#include <cppsim / state.hpp>
#include <cppsim / gate_factory.hpp>
#include <cppsim / gate_merge.hpp>

int main () {
unsigned int index = 3;
auto x_gate = gate :: X (index);

// Get matrix elements
// ComplexMatrix is ​​a complex matrix type made RowMajor with Eigen :: MatrixXcd
ComplexMatrix matrix;
x_gate-> set_matrix (matrix);
std :: cout << matrix << std :: endl;
return 0;
}
```


### Obtaining quantum gate information

Debug information of quantum gate can be displayed by pouring into <code> ostream </ code>. Note that if the gate matrix of the quantum gate is very large, it takes a long time. Quantum gates with dedicated functions do not display their own gate matrix.

```cpp
#include <iostream>
#include <cppsim / state.hpp>
#include <cppsim / gate_factory.hpp>
#include <cppsim / gate_merge.hpp>

int main () {

unsigned int index = 3;
auto x_gate = gate :: X (index);

std :: cout << x_gate << std :: endl;

delete x_gate;
return 0;
}
```


### Realization of general quantum gate
cppsim realizes various maps of quantum information in the following form.

#### Unitary operation
Implemented as a quantum gate.

#### Projection operator, Claus operator, etc.
Implemented as a quantum gate. In general, the norm of the quantum state is not preserved after the action. It can be generated by <code> DenseMatrix </ code> function.
```cpp
#define _USE_MATH_DEFINES
#include <cmath>
#include <cppsim / state.hpp>
#include <cppsim / gate_factory.hpp>
#include <cppsim / gate_merge.hpp>
#include <cppsim / gate_matrix.hpp>
#include <cppsim / gate_general.hpp>

int main () {
ComplexMatrix one_qubit_matrix (2, 2);
one_qubit_matrix << 0, 1, 1, 0;
auto one_qubit_gate = gate :: DenseMatrix (0, one_qubit_matrix);
std :: cout << one_qubit_gate << std :: endl;

ComplexMatrix two_qubit_matrix (4,4);
two_qubit_matrix <<
1, 0, 0, 0,
0, 1, 0, 0,
0, 0, 0, 1,
0, 0, 1, 0;
auto two_qubit_gate = gate :: DenseMatrix ({0,1}, two_qubit_matrix);
std :: cout << two_qubit_gate << std :: endl;
return 0;
}
```



#### Stochastic unitary operation
Use <code> Probabilistic </ code> function to create multiple unitary operations and probability distributions.
```cpp
#define _USE_MATH_DEFINES
#include <cmath>
#include <cppsim / state.hpp>
#include <cppsim / gate_factory.hpp>
#include <cppsim / gate_merge.hpp>
#include <cppsim / gate_matrix.hpp>
#include <cppsim / gate_general.hpp>

int main () {
unsigned int n = 5;
QuantumState state (n);
state.set_zero_state ();

unsigned int index = 3;
auto x_gate = gate :: X (index);
auto z_gate = gate :: Z (index);

auto probabilistic_xz_gate = gate :: Probabilistic ({0.1,0.2}, {x_gate, z_gate});
auto depolarizing_gate = gate :: DepolarizingNoise (index, 0.3);

depolarizing_gate-> update_quantum_state (& state);
probabilistic_xz_gate-> update_quantum_state (& state);
return 0;
}
```


#### CPTP-map
It is created by giving the <code> CPTP </ code> function as a list of Claus operators satisfying completeness.
```cpp
#define _USE_MATH_DEFINES
#include <cmath>
#include <cppsim / state.hpp>
#include <cppsim / gate_factory.hpp>
#include <cppsim / gate_merge.hpp>
#include <cppsim / gate_matrix.hpp>
#include <cppsim / gate_general.hpp>

int main () {
unsigned int n = 5;
QuantumState state (n);
state.set_zero_state ();

unsigned int index = 3;
auto p0 = gate :: P0 (index);
auto p1_fix = gate :: merge (gate :: P1 (index), gate :: X (index));

auto correction = gate :: CPTP ({p0, p1_fix});
auto noise = gate :: BitFlipNoise (index, 0.1);

noise-> update_quantum_state (& state);
correction-> update_quantum_state (& state);
return 0;
}
```


#### POVM
Since it is the same as Instrument in numerical calculation, it is realized as Instrument.

#### Instrument
Instrument is an operation to get the subscript of the Claus operator that acts randomly, in addition to the general CPTP-map operation. For example, measurement on the Z basis is equivalent to acting on a CPTP-map consisting of <code> P0 </ code> and <code> P1 </ code> and knowing which one worked.
In cppsim, the <code> Instrument </ code> function is realized by specifying the information of the CPTP-map and the address of the classic register in which the subscript of the operated Claus operator is written.

```cpp
#define _USE_MATH_DEFINES
#include <cmath>
#include <cppsim / state.hpp>
#include <cppsim / gate_factory.hpp>
#include <cppsim / gate_merge.hpp>
#include <cppsim / gate_matrix.hpp>
#include <cppsim / gate_general.hpp>

int main () {
auto gate00 = gate :: merge (gate :: P0 (0), gate :: P0 (1));
auto gate01 = gate :: merge (gate :: P0 (0), gate :: P1 (1));
auto gate10 = gate :: merge (gate :: P1 (0), gate :: P0 (1));
auto gate11 = gate :: merge (gate :: P1 (0), gate :: P1 (1));

std :: vector <QuantumGateBase *> gate_list = {gate00, gate01, gate10, gate11};
unsigned int classical_pos = 0;
auto gate = gate :: Instrument (gate_list, classical_pos);

QuantumState state (2);
state.set_Haar_random_state ();

std :: cout << state << std :: endl;
gate-> update_quantum_state (& state);
unsigned int result = state.get_classical_value (classical_pos);
std :: cout << state << std :: endl;
std :: cout << result << std :: endl;
return 0;
}
```


#### Adaptive
Perform or not perform operations depending on the value written to the classic register. cppsim implements this by specifying a function that takes a register of type <code> [unsigned int] </ code> as an argument and returns a <code> bool </ code> type.

```cpp
#define _USE_MATH_DEFINES
#include <cmath>
#include <cppsim / state.hpp>
#include <cppsim / gate_factory.hpp>
#include <cppsim / gate_merge.hpp>
#include <cppsim / gate_matrix.hpp>
#include <cppsim / gate_general.hpp>

int main () {
unsigned int n = 5;
QuantumState state (n);
state.set_zero_state ();

unsigned int index = 3;
auto h = gate :: H (index);
h-> update_quantum_state (& state);

auto meas = gate :: Measurement (index, 0);
meas-> update_quantum_state (& state);

auto condition = [] (const std :: vector <UINT> reg) {
return reg [0] == 1;
};
auto correction = gate :: Adaptive (gate :: X (index), condition);
correction-> update_quantum_state (& state);
return 0;
}
```

#### CP-map
If Kraus-rank is 1, treat it as a single Claus operator as described above. In other cases, after adjusting the Claus operator to be TP, adjust the <code> Identity </ code> operator multiplied by a constant with the <code> multiply_scalar </ code> function to adjust please.



## Quantum circuit

### Quantum circuit configuration
A quantum circuit is represented as a set of quantum gates.
For example, a quantum circuit can be configured as follows.

```cpp
#include <cppsim / state.hpp>
#include <cppsim / gate_factory.hpp>
#include <cppsim / circuit.hpp>

int main () {
unsigned int n = 5;
QuantumState state (n);
state.set_zero_state ();

// define quantum circuit
QuantumCircuit circuit (n);

// Add a gate to the quantum circuit
for (int i = 0; i <n; ++ i) {
circuit.add_H_gate (i);
}

// You can add your own defined gate
for (int i = 0; i <n; ++ i) {
circuit.add_gate (gate :: H (i));
}

// act on quantum circuits
circuit.update_quantum_state (& state);
return 0;
}
```

The quantum circuit added by <code> add_gate </ code> is released when the quantum circuit is released. Therefore, the assigned gate cannot be reused.
If you want to reuse the gate given as an argument, use the <code> add_gate_copy </ code> function. However, in this case you need to release the gate yourself.

### Quantum circuit optimization

By combining quantum gates into a single quantum gate, the number of quantum gates can be reduced and the time required for numerical calculations can be reduced. (Of course, when the number of target qubits increases, or when a quantum gate with a dedicated function is synthesized into a quantum gate without a dedicated function, the total calculation time may decrease. It depends.)

The following code uses the <code> optimize </ code> function to repeat the greedy synthesis of the quantum gate of the quantum circuit until the target qubit becomes three.

```cpp
#include <cppsim / state.hpp>
#include <cppsim / gate_factory.hpp>
#include <cppsim / circuit.hpp>
#include <cppsim / circuit_optimizer.hpp>

int main () {
unsigned int n = 5;
unsigned int depth = 10;
QuantumCircuit circuit (n);
for (int d = 0; d <depth; ++ d) {
for (int i = 0; i <n; ++ i) {
circuit.add_gate (gate :: H (i));
}
}

// Quantum circuit optimization
QuantumCircuitOptimizer opt;
unsigned int max_block_size = 3;
opt.optimize (& circuit, max_block_size);
return 0;
}
```

### Information debugging of quantum circuits
Like a quantum gate, a quantum circuit can display debug information by flowing it into <code> ostream </ code>.

```cpp
#include <cppsim / state.hpp>
#include <cppsim / gate_factory.hpp>
#include <cppsim / circuit.hpp>

int main () {
unsigned int n = 5;
unsigned int depth = 10;
QuantumCircuit circuit (n);
for (int d = 0; d <depth; ++ d) {
for (int i = 0; i <n; ++ i) {
circuit.add_gate (gate :: H (i));
}
}

// Output quantum circuit information
std :: cout << circuit << std :: endl;
return 0;
}
```

## Observable

### Generate observables
Observables are represented as a set of Pauli operators.
The Pauli operator can be defined as follows:
```cpp
#include <cppsim / observable.hpp>
#include <string>

int main () {
unsigned int n = 5;
double coef = 2.0;
std :: string Pauli_string = "X 0 X 1 Y 2 Z 4";
Observable observable (n);
observable.add_operator (coef, Pauli_string.c_str ());
return 0;
}
```

### Link with OpenFermion
In addition, observables can be generated from files in the following formats generated using OpenFermion. At this time, the observable is the minimum size necessary to compose it. For example, it is possible to read observables obtained using openfermion as shown below and generate observables.
```python
from openfermion.ops import FermionOperator
from openfermion.transforms import bravyi_kitaev

h_00 = h_11 = -1.252477
h_22 = h_33 = -0.475934
h_0110 = h_1001 = 0.674493
h_2332 = h_3323 = 0.697397
h_0220 = h_0330 = h_1221 = h_1331 = h_2002 = h_3003 = h_2112 = h_3113 = 0.663472
h_0202 = h_1313 = h_2130 = h_2310 = h_0312 = h_0132 = 0.181287

fermion_operator = FermionOperator ('0 ^ 0', h_00)
fermion_operator + = FermionOperator ('1 ^ 1', h_11)
fermion_operator + = FermionOperator ('2 ^ 2', h_22)
fermion_operator + = FermionOperator ('3 ^ 3', h_33)

fermion_operator + = FermionOperator ('0 ^ 1 ^ 1 0', h_0110)
fermion_operator + = FermionOperator ('2 ^ 3 ^ 3 2', h_2332)
fermion_operator + = FermionOperator ('0 ^ 3 ^ 3 0', h_0330)
fermion_operator + = FermionOperator ('1 ^ 2 ^ 2 1', h_1221)

fermion_operator + = FermionOperator ('0 ^ 2 ^ 2 0', h_0220-h_0202)
fermion_operator + = FermionOperator ('1 ^ 3 ^ 3 1', h_1331-h_1313)

fermion_operator + = FermionOperator ('0 ^ 1 ^ 3 2', h_0132)
fermion_operator + = FermionOperator ('2 ^ 3 ^ 1 0', h_0132)

fermion_operator + = FermionOperator ('0 ^ 3 ^ 1 2', h_0312)
fermion_operator + = FermionOperator ('2 ^ 1 ^ 3 0', h_0312)

## Bravyi-Kitaev transformation
bk_operator = bravyi_kitaev (fermion_operator)

## output
fp = open ("H2.txt", 'w')
fp.write (str (bk_operator))
fp.close ()
```
At this time, the <code> H2.txt </ code> file generated by the above python code has the following format.
```txt
(-0.8126100000000005 + 0j) [] +
(0.04532175 + 0j) [X0 Z1 X2] +
(0.04532175 + 0j) [X0 Z1 X2 Z3] +
(0.04532175 + 0j) [Y0 Z1 Y2] +
(0.04532175 + 0j) [Y0 Z1 Y2 Z3] +
(0.17120100000000002 + 0j) [Z0] +
(0.17120100000000002 + 0j) [Z0 Z1] +
(0.165868 + 0j) [Z0 Z1 Z2] +
(0.165868 + 0j) [Z0 Z1 Z2 Z3] +
(0.12054625 + 0j) [Z0 Z2] +
(0.12054625 + 0j) [Z0 Z2 Z3] +
(0.16862325 + 0j) [Z1] +
(-0.22279649999999998 + 0j) [Z1 Z2 Z3] +
(0.17434925 + 0j) [Z1 Z3] +
(-0.22279649999999998 + 0j) [Z2]
```
To create an observable from such a file, you can create an observable through a function as follows:

```cpp
#include <cppsim / observable.hpp>
#include <string>

int main () {
unsigned int n = 5;
std :: string filename = "H2.txt";
Observable * observable = observable :: create_observable_from_openfermion_file (filename);
delete observable;
return 0;
}
```


### Observable evaluation
You can evaluate the expected value of the observable against the state.
```cpp
#include <cppsim / observable.hpp>
#include <cppsim / state.hpp>
#include <string>

int main () {
unsigned int n = 5;
double coef = 2.0;
std :: string Pauli_string = "X 0 X 1 Y 2 Z 4";
Observable observable (n);
observable.add_operator (coef, Pauli_string.c_str ());
</ S> </ s> </ s>
QuantumState state (n);
observable.get_expectation_value (& state);
return 0;
}
```

### Observable rotation
Rotate observable \ f $ H \ f $ \ f $ e ^ {i \ theta H} \ f $ by Trotter expansion. <code> num_repeats </ code> defaults to the following code, but can be specified by the user.
```cpp
#include <cppsim / circuit.hpp>
#include <cppsim / state.hpp>
#include <cppsim / observable.hpp>

int main () {
UINT n;
UINT num_repeats;
double theta = 0.1;
Observable * observable = observable :: create_observable_from_openfermion_file ("../ test / cppsim / H2.txt");

n = observable-> get_qubit_count ();
QuantumState state (n);
state.set_computational_basis (0);

QuantumCircuit circuit (n);
num_repeats = (UINT) std :: ceil (theta * (double) n * 100.);
circuit.add_observable_rotation_gate (* observable, theta, num_repeats);
circuit.update_quantum_state (& state);

auto result = observable-> get_expectation_value (& state);
std :: cout << result << std :: endl;
delete observable;
return 0;
}
```


## Variational quantum circuit
Defining a quantum circuit as a ParametricQuantumCircuit class allows you to use some functions that are useful for optimizing quantum circuits using variational methods, in addition to the usual functions of the QuantumCircuit class.

### Variational quantum circuit usage example

Quantum gates with one rotation angle (X-rot, Y-rot, Z-rot, multi_qubit_pauli_rotation) can be added to quantum circuits as parametric quantum gates. For quantum gates added as parametric gates, the number of parametric gates can be extracted after the quantum circuit is configured, and the rotation angle can be changed later.

```cpp
#include <cppsim / state.hpp>
#include <vqcsim / parametric_circuit.hpp>
#include <cppsim / utility.hpp>

int main () {
const UINT n = 3;
const UINT depth = 10;

// create n-qubit parametric circuit
ParametricQuantumCircuit * circuit = new ParametricQuantumCircuit (n);
Random random;
for (UINT d = 0; d <depth; ++ d) {
// add parametric X, Y, Z gate with random initial rotation angle
for (UINT i = 0; i <n; ++ i) {
circuit-> add_parametric_RX_gate (i, random.uniform ());
circuit-> add_parametric_RY_gate (i, random.uniform ());
circuit-> add_parametric_RZ_gate (i, random.uniform ());
}
// add neighboring two-qubit ZZ rotation
for (UINT i = d% 2; i + 1 <n; i + = 2) {
circuit-> add_parametric_multi_Pauli_rotation_gate ({i, i + 1}, {3,3}, random.uniform ());
}
}

// get parameter count
UINT param_count = circuit-> get_parameter_count ();

// get current parameter, and set shifted parameter
for (UINT p = 0; p <param_count; ++ p) {
double current_angle = circuit-> get_parameter (p);
circuit-> set_parameter (p, current_angle + random.uniform ());
}

// create quantum state and update
QuantumState state (n);
circuit-> update_quantum_state (& state);

// output state and circuit info
std :: cout << state << std :: endl;
std :: cout << circuit << std :: endl;

// release quantum circuit
delete circuit;
}
```

