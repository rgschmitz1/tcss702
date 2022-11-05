from minio import Minio
import json
import os
import subprocess
from .Inspector import Inspector
from .sebs import SeBS

# Create Minio client
def minio_client():
    with open("/var/openfaas/secrets/minio-access-key") as f:
        access_key = f.read()
    with open("/var/openfaas/secrets/minio-secret-key") as f:
        secret_key = f.read()

    return Minio(os.environ['minio_hostname'],
        access_key=access_key,
        secret_key=secret_key,
        secure=False)

def handle(event, context):
    # Start inspector
    inspector = Inspector()
    inspector.inspectAll()

    # Parse input json data
    body = json.loads(event.body)

    #payload = body['payload']
    #bucket = body['bucket']

    # Get minio client object to download/upload data
    mc = minio_client()
    # Initalize SeBS object
    sebs = SeBS(mc)

    # Download sample
    #mc.fget_object(bucket, inputfile, f'/tmp/{inputfile}')

    # Uncompress sample
    #uncompress_sample = subprocess.run(f'tar -xf /tmp/{inputfile} -C /tmp', shell=True)
    #print(uncompress_sample, flush=True)

    # Put aligned normal and tumor bam in minio
    #mc.fput_object(bucket, aligned_bam, f'/tmp/{aligned_bam}')

    # Collect inspector deltas
    inspector.inspectAllDeltas()

    iret = inspector.finish()
    ret = {
        "status": 200,
        "body": iret
    }
    return ret
