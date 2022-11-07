from minio import Minio
import json
import os
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

    # Get minio client object to download/upload data
    mc = minio_client()

    # Initalize SeBS object
    sebs = SeBS(mc)

    # Use dictionary to switch between functions
    fn = {"dna_visualization": sebs.dna_visualization,
          "graph_bfs": sebs.graph_bfs,
          "graph_mst": sebs.graph_mst,
          "graph_pagerank", sebs.graph_pagerank}

    fn_name = body['fn']
    if fn_name == "dna_visualization":
        bucket = body['bucket']
        key = body['key']
        fn_ret = fn[fn_name](key, bucket, bucket)
    else:
        size = body['size']
        fn_ret = fn[fn_name](size)
    print(fn_ret, flush=True)

    # Collect inspector deltas
    inspector.inspectAllDeltas()

    iret = inspector.finish()
    ret = {
        "status": 200,
        "body": iret
    }
    return ret
