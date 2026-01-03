from qiskit_ibm_runtime import QiskitRuntimeService
def save_account(token, instance):
    
    service1 = QiskitRuntimeService.save_account(
        token=token,
        instance=instance,
        name='qdc-2025',
        overwrite=True)
    service2 = QiskitRuntimeService.save_account(
        token=token,
        instance=instance,
        name='r2p-2025',
        overwrite=True)