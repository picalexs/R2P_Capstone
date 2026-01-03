#!/usr/bin/env python3
import numpy as np
from json.encoder import JSONEncoder
from json.decoder import JSONDecoder
from functools import partial
import os

from qiskit_ibm_runtime import QiskitRuntimeService
from qiskit_serverless import distribute_task, get_arguments, get, save_result
from qiskit_addon_sqd.fermion import SCIResult, diagonalize_fermionic_hamiltonian, solve_sci

from functools import partial
from qiskit_addon_sqd.fermion import (
    SCIResult,
    diagonalize_fermionic_hamiltonian,
)

### Argument retrieval
args                      = get_arguments() 

data                      = args["data"]                           # Chemistry Data
energy_tol                = args["energy_tol"]                     # SQD option
occupancies_tol           = args["occupancies_tol"]                # SQD option
max_iterations            = args["max_iterations"]                 # SQD option
symmetrize_spin           = args["symmetrize_spin"]                # Eigenstate solver option
carryover_threshold       = args["carryover_threshold"]            # Eigenstate solver option
num_batches               = args["num_batches"]                    # Eigenstate solver option
samples_per_batch         = args["samples_per_batch"]              # Eigenstate solver option
max_cycle                 = args["max_cycle"]                      # Eigenstate solver option
mem                       = args["mem"]                            # Memory per Worker

# --- fan‑out target: 1 CPU + mem GB RAM per call -------------
@distribute_task(target={"cpu": 1, "mem": mem * 1024**3})
def _solve_sci_worker(ix, ci_strs, one_body_tensor, two_body_tensor, norb, nelec, spin_sq):
    print(f">>>>> HELLO! I AM WORKER {ix}! Leave a 30% tip? (y/n)")
    res = solve_sci(ci_strs,
                     one_body_tensor,
                     two_body_tensor,
                     norb=norb,
                     nelec=nelec,
                     spin_sq=spin_sq)
    
    print(f">>>>> WORKER {ix} COMPLETE")     
    return res

def distribute_solve_sci_batch(
    ci_strings: list[tuple[np.ndarray, np.ndarray]],
    one_body_tensor: np.ndarray,
    two_body_tensor: np.ndarray,
    norb: int,
    nelec: tuple[int, int],
    *,
    spin_sq: float | None = None,
    **kwargs,
) -> list[SCIResult]:
    """Diagonalize Hamiltonian in subspaces, parallelizing across vCPUs in the Serverless environment.

    Args:
        ci_strings: List of pairs (strings_a, strings_b) of arrays of spin-alpha CI
            strings and spin-beta CI strings whose Cartesian product give the basis of
            the subspace in which to perform a diagonalization.
        one_body_tensor: The one-body tensor of the Hamiltonian.
        two_body_tensor: The two-body tensor of the Hamiltonian.
        norb: The number of spatial orbitals.
        nelec: The numbers of alpha and beta electrons.
        spin_sq: Target value for the total spin squared for the ground state.
            If ``None``, no spin will be imposed.
        **kwargs: Keyword arguments to pass to `pyscf.fci.selected_ci.kernel_fixed_space <https://pyscf.org/pyscf_api_docs/pyscf.fci.html#pyscf.fci.selected_ci.kernel_fixed_space>`_

    Returns:
        The results of the diagonalizations in the subspaces given by ci_strings.
    """
    inputs = [(ix, ci_strs, one_body_tensor, two_body_tensor, norb, nelec, spin_sq) for ix, ci_strs in enumerate(ci_strings)] 
    
    # fan‑out: spawn one worker per input tuple
    print(">>>>> ENTERING WORKER FAN-OUT")
    refs = [_solve_sci_worker(*input_) for input_ in inputs] 
    print(">>>>> WAITING ON WORKERS TO FINISH TASKS")

    # fan‑in: block until every worker finishes
    results = get(refs) 
    print(">>>>> DISTRIBUTED JOBS COMPLETED")

    return results 

# A caveat of executing a Python program remotely is that the inputs to the remote program must be passed over an internet network. Similarly, the outputs 
# must be passed back to the local program via the same structure. Python objects are not always able to be passed over a network, and must be encoded in a 
# JSON serializable format. 
i_data = JSONDecoder().decode(data)

# i_data has all of the information needed from the local program to pick up where the computation left off
# after its submission to the remote environment. 
[job_id, hcore, eri, num_orbitals, nuclear_repulsion_energy, num_elec_a, num_elec_b] = i_data

# Re-convert data back into numpy formay, after serialization 
hcore                     = np.array(hcore)
eri                       = np.array(eri)
nuclear_repulsion_energy  = np.float64(nuclear_repulsion_energy)

# Instantiate Runtume Service to retrieve the bitstrings from the QPU job. We provided these credntials upon Serverless setup. 
service = QiskitRuntimeService( 
    channel=os.environ.get("QISKIT_IBM_CHANNEL"),
    token=os.environ.get("QISKIT_IBM_TOKEN"),
    instance=os.environ.get("QISKIT_IBM_INSTANCE")
)

# retrieving the QPU job data from the Serverless side 
job = service.job(job_id)
primitive_result = job.result()
pub_result = primitive_result[0]
bit_array = pub_result.data.meas # Getting the bitstrings

# Pass options to the built-in eigensolver
sci_solver = partial(distribute_solve_sci_batch, spin_sq=0.0, max_cycle=max_cycle)

# List to capture intermediate results
result_history = []

def callback(results: list[SCIResult]):
    result_history.append(results)
    iteration = len(result_history)
    print(f">>>>> SQD ITERATION {iteration}")
    for i, result in enumerate(results):
        print(f">>>>> SUBSAMPLE {i}")
        print(f">>>>> \tENERGY: {result.energy + nuclear_repulsion_energy}")
        print(
                f">>>>> \tSUBSPACE DIMENSION: {np.prod(result.sci_state.amplitudes.shape)}"
        )

result = diagonalize_fermionic_hamiltonian(
    hcore,
    eri,
    bit_array,
    samples_per_batch=samples_per_batch,
    norb=num_orbitals,
    nelec=(num_elec_a, num_elec_b),
    num_batches=num_batches,
    energy_tol=energy_tol,
    occupancies_tol=occupancies_tol,
    max_iterations=max_iterations,
    sci_solver=sci_solver,
    symmetrize_spin=symmetrize_spin,
    carryover_threshold=carryover_threshold,
    callback=callback,
    seed=12345,
    )

print(">>>>> EXACT DIAGONALIZATION COMPLETE. CLEANING UP, SERIALIZING DATA.")
# Numpy arrays are not JSON serializable, convert them to List objects before using the JSONEncoder
o_data = JSONEncoder().encode([result.energy+nuclear_repulsion_energy,
                               result.energy,
                               result.rdm1.tolist(),
                               result.rdm2.tolist(),
                               [x.tolist() for x in result.orbital_occupancies],
                               [result.sci_state.nelec,
                               result.sci_state.norb,
                               [x.tolist() for x in result.sci_state.orbital_occupancies()],
                               [x.tolist() for x in result.sci_state.rdm()]]])

# JSON-safe package
save_result({"outputs": o_data})  # single JSON blob returned to client

