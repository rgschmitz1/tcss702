from minio import Minio
import json
import os
from .Inspector import Inspector

def handle(event, context):
    with open("/var/openfaas/secrets/minio-access-key") as f:
        access_key = f.read()
    with open("/var/openfaas/secrets/minio-secret-key") as f:
        secret_key = f.read()

    mc = Minio(os.environ['minio_hostname'],
               access_key=access_key,
               secret_key=secret_key,
               secure=False)

    # Collect data
    inspector = Inspector()
    inspector.inspectAll()
    # Add custom message and finish the function
#    if "startWallClock" in event:
#        inspector.addAttribute("startWallClock", event['startWallClock'])

    body = json.loads(event.body)
#    print(body['fn'], flush=True)

#    fn = {"p": tm.preprocess,
#          "t": tm.train,
#          "q": tm.query}

#    fn[body['fn']]()

    inspector.inspectAllDeltas()
    # Include functionName
#    inspector.addAttribute("functionName", fn[body['fn']].__name__)

    iret = inspector.finish()
    ret = {
        "status": 200,
        "body": iret
    }
    return ret
