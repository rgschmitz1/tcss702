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
    thread_cnt = body['thread_cnt']

    # Get minio client object to download/upload data
    mc = minio_client()

    # Initalize BWA object
    bwa = BWA()

    # Align input file using bwa
    ret = bwa.process(mc, inputfile, bucket, thread_cnt)

    # Collect inspector deltas
    inspector.inspectAllDeltas()

    # Include functionName
    inspector.addAttribute("functionName", 'bwa-mem')

    iret = inspector.finish()

    if ret['status'] == 0:
        # Removed aligned sample, we have no use for it after workload is complete
        mc.remove_object(bucket, ret['body'])
        ret = {
            "body": {
                "status": 200,
                "saaf": iret
            }
        }
    else:
        ret = {
            "body": {
                "status": ret['status'],
                "message": ret['body']
            }
        }

    return ret
