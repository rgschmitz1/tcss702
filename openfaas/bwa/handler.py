from minio import Minio
import json
import os
from .SAAF import Inspector
from .BWA import BWA

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

    # assume input file file is tar and compressed
    inputfile = body['inputfile']
    bucket = body['bucket']

    # Get minio client object to download/upload data
    mc = minio_client()

    # Initalize BWA object
    bwa = BWA()

    # Align input file using bwa
    bwa.process(mc, inputfile, bucket)

    # Collect inspector deltas
    inspector.inspectAllDeltas()

    # Include functionName
    inspector.addAttribute("functionName", 'bwa-mem')

    iret = inspector.finish()

    # Removed aligned sample, we have no use for it after workload is complete
    mc.remove_object(bucket, bwa.get_aligned_sample())

    # Construct json return from function
    ret = {
        "status": 200,
        "body": iret
    }
    return ret