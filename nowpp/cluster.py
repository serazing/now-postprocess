import os
from dask.distributed import LocalCluster, Client
from dask_jobqueue import PBSCluster

PROJECT = os.environ['PROJECT']
USER = os.environ['USER']


def vdi_cluster(**kwargs):
    """
    Start a local cluster on VDI using dask

    Parameters
    ----------
    **kwargs : optional
        Additional parameters to pass to dask.distributed.LocalCluster

    Returns
    -------
    cluster : dask.distributed.LocalCluster
        The dask local cluster
    client : dask.distributed.Client
        The dask client connected to the local cluster
    """
    local_dir = "/local/%s/%s/dask-workers/" % (PROJECT, USER)
    cluster = LocalCluster(local_dir=local_dir, **kwargs)
    client = Client(cluster)
    return cluster, client


def raijin_cluster(cores=14, memory='64GB', processes=1, queue='expressbw',
                   walltime='01:00:00', **kwargs):
    """
    Start a local cluster on VDI using dask

    Parameters
    ----------
    cores : int, optional
        The number of cores to use on Raijin, must not be larger than the
        number of cores available on one node
    memory : str, optional
        The memory asked on the node
    processes : int or str, optional
        Number of processes
    queue : str, optional
        The name of the queue. Default is to run on the Broadwell nodes with
        the express queue
    walltime : str, optional
        The time limit for the job. Default is 1 h.
    **kwargs : optional
        Additional parameters to pass to dask_jobqueue.PBS_cluster

    Returns
    -------
    cluster :dask_jobqueue.PBS_cluster
        The dask PBS cluster
    client : dask.distributed.Client
        The dask client connected to the local cluster
    """
    cluster = PBSCluster(cores=cores,
                         memory=memory,
                         processes=processes,
                         ip=get_interface_ip('ib0'),
                         dashboard_address=get_interface_ip('vlan192'),
                         **kwargs)
    client = Client(cluster)
    cluster.job_header = (('#!/usr/bin/env bash\n'
                           '#PBS -P %s\n'
                           '#PBS -N dask-worker\n'
                           '#PBS -q %s\n'
                           '#PBS -l ncpus=%s\n'
                           '#PBS -l mem=%s\n'
                           '#PBS -l walltime=%s\n')
                          % (PROJECT, queue, cores, memory, walltime)
                          + 'JOB_ID=${PBS_JOBID%.*}'
                          )
    return cluster, client


def get_interface_ip(ifname):
    import socket
    import fcntl
    import struct
    s = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
    stpack = struct.pack('256s', bytes(ifname[:15], 'utf-8'))
    return socket.inet_ntoa(fcntl.ioctl(s.fileno(), 0x8915, stpack)[20:24])