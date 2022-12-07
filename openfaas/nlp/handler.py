from minio import Minio
import json
import os
from .SAAF import Inspector
from .topic_model import topic_model

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
    mc = minio_client()

    tm = topic_model(mc)

    # Collect data
    inspector = Inspector()
    inspector.inspectAll()

    body = json.loads(event.body)
    fn_name = body['fn']

    fn = {"preprocess": tm.preprocess,
          "train": tm.train,
          "query": tm.query}

    fn[fn_name]()

    inspector.inspectAllDeltas()

    # Include functionName
    inspector.addAttribute("functionName", f'nlp-{fn_name}')

    iret = inspector.finish()

    ret = {
        "body": {
            "status": 200,
            "saaf": iret
        }
    }
    return ret
