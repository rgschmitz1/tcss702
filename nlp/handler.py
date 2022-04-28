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

    # Trying to figure out event
    body = json.loads(event.body)
    print(body['fn'], flush=True)
#    print(event.headers, flush=True)
#    print(event.method, flush=True)
#    print(event.query, flush=True)
#    print(event.path, flush=True)

    # Just execute the first function and see what happens
    fn = {1: tm.lambda_function_1,
          2: tm.lambda_function_2,
          3: tm.lambda_function_3}

    fn[body['fn']]()

    inspector.inspectAllDeltas()
    # Include lambdaName and functionName
#    inspector.addAttribute("lambdaName", inspector.getAttribute('functionName'))
#    inspector.addAttribute("functionName", event['function_name'])

    iret = inspector.finish()
    ret = {
        "status": 200,
        "body": iret
    }
    print(ret, flush=True)
    return ret
