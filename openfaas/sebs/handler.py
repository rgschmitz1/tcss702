from minio import Minio
import json
import os
from .SAAF import Inspector
from .SeBS import SeBS

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

    # Initalize SeBS object
    sebs = SeBS()

    # Use dictionary to switch between functions
    fn = {
        "dna_visualization": sebs.dna_visualization,
        "graph_bfs": sebs.graph_bfs,
        "graph_mst": sebs.graph_mst,
        "graph_pagerank": sebs.graph_pagerank
    }

    fn_name = body['fn']
    if fn_name == "dna_visualization":
        # Get minio client object to download/upload data
        mc = minio_client()
        key = body['key']
        bucket = body['bucket']
        fn[fn_name](mc, key, bucket, bucket)
    else:
        size = body['size']
        fn[fn_name](size)

    # Collect inspector deltas
    inspector.inspectAllDeltas()

    # Include function name
    inspector.addAttribute("functionName", f'sebs-{fn_name}')

    iret = inspector.finish()

    # If DNA visualization was executed, remove output object
    if fn_name == "dna_visualization":
        mc.remove_object(bucket, sebs.get_dna_visualization_output())

    ret = {
        "body": {
            "status": 200,
            "saaf": iret
        }
    }
    return ret
