from minio import Minio
import os
from .Inspector import Inspector
from .topic_model import topic_model

def handle(ret):
#def handle(event, context):

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
#    topic_model.run[event['function_name']]()
    #print(event)
    tmret = tm.lambda_function_1()
#    print(event.headers.get('Content-Type'))

    inspector.inspectAllDeltas()
    # Include lambdaName and functionName
#    inspector.addAttribute("lambdaName", inspector.getAttribute('functionName'))
#    inspector.addAttribute("functionName", event['function_name'])

    iret = inspector.finish()
    ret = {
        "status": 200,
        "body": tmret
    }
    return ret
