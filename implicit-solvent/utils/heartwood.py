#!/usr/bin/env python3

# Heartwood Algorithm subroutines
import os
import numpy as np
import pyscf
from pyscf import ao2mo, cc
from functools import reduce
import ffsim
from qiskit import QuantumCircuit, QuantumRegister

from json.encoder import JSONEncoder
from json.decoder import JSONDecoder
import time

from functools import partial
from qiskit_addon_sqd.fermion import (
    SCIResult,
    diagonalize_fermionic_hamiltonian,
    solve_sci_batch,
)

from utils.helpers import feedback_serverless

def update_rdm(casci_object, dmas):
    """
      Inputs:
        mc   -> CASCI object
        dmas -> Spin-summed 1-particle reduced density matrix

      This function returns the CASCI/SQD one-body density matrix in the full basis of atomic orbitals, written as the sum (last line) of two terms:
         - a contribution from the core orbitals, np.dot(mocore, mocore.conj().T) * 2, (core = inactive and doubly-occupied)
         - a contribution from the active-space orbitals and electrons (dmas) rotated from the active-space to the AO basis (the reduce operation)

      Outputs:
        rho_approximation: The CASCI/SQD one-body density matrix in the full basis of atomic orbitals

    """
    mo_coeff = casci_object.mo_coeff
    ncore = casci_object.ncore
    ncas = casci_object.ncas
    mocore = mo_coeff[:,:ncore]
    mocas = mo_coeff[:,ncore:ncore+ncas]
    dm1 = np.dot(mocore, mocore.conj().T) * 2

    rho_approximation = dm1 + reduce(np.dot, (mocas, dmas, mocas.conj().T))
    return rho_approximation

def run_active_space_calculation(h1e_cas, h2e_cas, norb, ne_act, orbs, fermilevel, ecore):
    # ----- perform an HF and a CCSD calculation in the active space
    from pyscf import tools
    from datetime import datetime

    now = datetime.now().strftime("%H:%M:%S")
    print(">>>>> ACTIVE SPACE CALCULATIONS ")
    tools.fcidump.from_integrals(f'as_fcidump_{now}.txt', h1e_cas, h2e_cas, norb, ne_act, ms=0, nuc=ecore) # Forcefully represents the active space in the correct structure
    mf_as = tools.fcidump.to_scf(f'as_fcidump_{now}.txt')
    os.remove(f'as_fcidump_{now}.txt')
    mf_as.kernel()
    print(">>>>> RUNNING CCSD")

    mf_cc = cc.CCSD(mf_as)
    mf_cc.kernel()
    orbts = mf_as.mo_coeff
    t1,t2 = mf_cc.t1, mf_cc.t2
    print(">>>>> UPDATED t1, t2 PARAMETERS")

    # ----- update the HF orbitals
    active = list(range(fermilevel-ne_act//2,fermilevel-ne_act//2+norb))
    orbs[:,active] = np.dot(orbs[:,active],orbts)
    return orbs, t1, t2

def get_lucj(norb, num_elec_a, num_elec_b, t1, t2, n_reps=1):
    print(f">>>>> CONSTRUCTING LUCJ CIRCUIT")

    alpha_alpha_indices = [(p, p + 1) for p in range(norb - 1)]
    alpha_beta_indices = [(p, p) for p in range(0, norb, 4)]

    ucj_op = ffsim.UCJOpSpinBalanced.from_t_amplitudes(
        t1=t1, # <---- Update t1 each loop
        t2=t2, # <---- Update t2 each loop
        n_reps=n_reps,
        interaction_pairs=(alpha_alpha_indices, alpha_beta_indices),
        )
    nelec = (num_elec_a, num_elec_b)

    # create an empty quantum circuit
    qubits = QuantumRegister(2 * norb, name="q")
    circuit = QuantumCircuit(qubits)

    # prepare Hartree-Fock state as the reference state and append it to the quantum circuit
    circuit.append(ffsim.qiskit.PrepareHartreeFockJW(norb, nelec), qubits)

    # apply the UCJ operator to the reference state
    circuit.append(ffsim.qiskit.UCJOpSpinBalancedJW(ucj_op), qubits)
    circuit.measure_all()

    return circuit

# Classical Diagonalization Engine sent to HPC
def classically_diagonalize(bit_array=None,                   # Bit string array (only needed if locally processing data)
                            nuclear_repulsion_energy=None,    # Electronic energy from the core orbitals
                            hcore=None,                       # 1-electron hamiltonian integrals
                            eri=None,                         # 2-electron hamiltonian integrals
                            num_orbitals=None,                # Number of spatial orbitals
                            nelec=None,                       # Number of electrons
                            num_elec_a=None,                  # Alpha orbitals
                            num_elec_b=None,                  # Beta orbitals
                            job_id=None,                      # QPU bitstring Job ID
                            client=None,                      # Diagonalization engine worker
                            energy_tol=1e-4,                  # SQD option
                            occupancies_tol=1e-3,             # SQD option
                            max_iterations=12,                # SQD option
                            num_batches = 8,                  # Eigenstate solver option
                            samples_per_batch=300,            # Eigenstate solver option
                            symmetrize_spin=False,            # Eigenstate solver option
                            carryover_threshold=1e-5,         # Eigenstate solver option
                            max_cycle=200,                    # Eigenstate solver option
                            local=True,                       # Remote vs Local Diagonalization
                            mem=16.                           # Memory per Serverless Worker (Gb)
                            ):           

      print(">>>>> STARTING DIAGONALIZATION ENGINE ")
      # Pass options to the built-in eigensolver. If you just want to use the defaults,
      # you can omit this step, in which case you would not specify the sci_solver argument
      # in the call to diagonalize_fermionic_hamiltonian below.
      if local:
         sci_solver = partial(solve_sci_batch, spin_sq=0.0, max_cycle=max_cycle)

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
            nelec=(nelec//2, nelec//2),
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

         result = (result.energy, result.rdm1, result.rdm2)

      else:
         # Serverless Logic
         print(f">>>>> SENDING QISKIT RUNTIME JOB {job_id} TO QISKIT SERVERLESS")

         data = [job_id, hcore.tolist(), eri.tolist(), int(num_orbitals), float(nuclear_repulsion_energy), int(num_elec_a), int(num_elec_b)]

         # Encode the execution dependencies with the JSONEncoder
         data_e = JSONEncoder().encode(data)

         # Send to Serverless
         worker = client.load("diagonalization_engine") 
         serverless_job = worker.run(data=data_e,                        
                        energy_tol=energy_tol,                      # SQD option
                        occupancies_tol=occupancies_tol,            # SQD option
                        max_iterations=max_iterations,              # SQD option
                        symmetrize_spin=symmetrize_spin,            # Eigenstate solver option
                        carryover_threshold=carryover_threshold,    # Eigenstate solver option
                        num_batches=num_batches,                    # Eigenstate solver option
                        samples_per_batch=samples_per_batch,        # Eigenstate solver option
                        max_cycle=max_cycle,                        # Eigenstate solver option
                        mem=mem,                                    # Memory per Worker (Gb)      
                        )

         # Wait for the job to execute
         _ = feedback_serverless(serverless_job)

         o_data = JSONDecoder().decode(serverless_job.result()["outputs"])
         result = (o_data[1], np.array(o_data[2]), np.array(o_data[3]))

         print(f">>>>>>>>>> Active Space Energy: {o_data[1]}")
         print(f">>>>>>>>>> rdm1: {o_data[2]}")
         print(f">>>>>>>>>> rdm2: {o_data[3]}")

      return result


if __name__ == "__main__":
    raise SystemExit("This file is a utility module, not meant to be run directly.")
