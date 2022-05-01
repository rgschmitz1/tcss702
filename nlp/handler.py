from minio import Minio
import json
import os
from .Inspector import Inspector
from .topic_model import topic_model

#def handle(event):
def handle(event, context):
    mc = Minio(os.environ['minio_hostname'],
               access_key=os.environ['minio_access_key'],
               secret_key=os.environ['minio_secret_key'],
               secure=False)

    tm = topic_model(mc)

    # Collect data
    inspector = Inspector()
    inspector.inspectAll()
    # Add custom message and finish the function
#    if "startWallClock" in event:
#        inspector.addAttribute("startWallClock", event['startWallClock'])

    body = json.loads(event.body)
    print(body['fn'], flush=True)

    fn = {"p": tm.preprocess,
          "t": tm.train,
          "q": tm.query}

    fn[body['fn']]()

    inspector.inspectAllDeltas()
    # Include functionName
    inspector.addAttribute("functionName", body['fn'])

    iret = inspector.finish()
    ret = {
        "status": 200,
        "body": iret
    }
    return ret
