#!/usr/bin/env python3
import json
import time
from pathlib import Path

def init_results_file():
    molecules = ["Methylamine"]#, "Ethanol"]
    entries = []
    for m in molecules:
        entries.append({
            "molecule": m,
            "hartree_fock_E": float("nan"),
            "casci_E": float("nan"),
            "sqd_E": float("nan"),
        })
    filename = f"./results/results.json"
    with open(filename, "w") as f:
        json.dump(entries, f, indent=2, allow_nan=True)

def report_result(molecule_name: str,
                  hartree_fock_E: float,
                  casci_E: float,
                  sqd_E: list[tuple[int, float, str]]) -> None:
    """
    Call this function at the end of every demonstration to append your results to the results.yaml file. Results will be appended, not overwritten. So if you think you can improve your results
    then it is safe to do so, without the fear that a previous good run will be ruined.
    """
    # Validation checks
    assert molecule_name == "Methylamine" or molecule_name == "Ethanol"
    
    # Load existing file
    path = Path("./results/results.json")
    with open(path, "r") as f:
        data = json.load(f)
   
    # Update the matching molecule
    updated = False
    for entry in data:
        if entry["molecule"] == molecule_name:
            entry["hartree_fock_E"] = hartree_fock_E
            entry["casci_E"] = casci_E
            entry["sqd_E"] = sqd_E
            updated = True
            break

    if not updated:
        raise ValueError(f"Molecule {molecule_name} not found in {path}")

    # Write back
    with open(path, "w") as f:
        json.dump(data, f, indent=2, allow_nan=True)

def feedback_serverless(serverless_job):
    import time
    # Wait for the job to execute
    print(f">>>>> Serverless status: {serverless_job.job_id}")
    timer = 0
    while (timer < 10000):
        if serverless_job.status() == "QUEUED" \
            or serverless_job.status() == "INITIALIZING" \
                or serverless_job.status() == "RUNNING":
            print(f">>>>> [{timer}s] Serverless job {serverless_job.job_id}: \
                {serverless_job.status()}")
            time.sleep(10)
            timer += 10

        elif serverless_job.status() == "ERROR":
            print(f">>>>> Serverless job {serverless_job.job_id}: {serverless_job.status()}")
            print(">>>>> Logs:")
            print(serverless_job.logs())
            break

        elif serverless_job.status() == "DONE":
            print(f">>>>> Serverless job {serverless_job.job_id}: {serverless_job.status()}")
            break

        else:
            break

    return

if __name__ == "__main__":
    raise SystemExit("This file is a utility module, not meant to be run directly.")