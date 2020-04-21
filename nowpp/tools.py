def start_client(nci_machine, threads_per_worker=4, memory_limit=None):
    import os
    from dask.distributed import LocalCluster, Client
    if nci_machine == 'vdi':
        local_dir = "/local/e14/gs9353/dask-workers/"
        cluster = LocalCluster(threads_per_worker=threads_per_worker,
                               processes=True, local_directory=local_dir)
        client = Client(cluster)
    elif nci_machine == 'gadi':
        local_dir = os.path.join(os.environ['PBS_JOBFS'],
                                 'dask-worker-space')
        n_workers = int(os.environ['PBS_NCPUS']) // threads_per_worker
        if memory_limit is None:
            memory_limit = f'{3.9 * threads_per_worker}gb'
        client = Client(n_workers=n_workers,
                        threads_per_worker=threads_per_worker,
                        processes=True,
                        memory_limit=memory_limit,
                        local_directory=local_dir)
    elif nci_machine == 'rajin':
        raise ValueError('Raijin has been decommited')
    else:
        raise ValueError('No such machine')
    return client
