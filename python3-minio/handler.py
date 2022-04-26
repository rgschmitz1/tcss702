from minio import Minio
#import requests
import json
import os

def handle(req):
    """handle a request to the function
    Args:
        req (str): request body
    """
    mc = Minio(os.environ['minio_hostname'],
               access_key=os.environ['minio_access_key'],
               secret_key=os.environ['minio_secret_key'],
               secure=False)
    err = None
    try:
        mc.make_bucket('abucket')
        buckets = mc.list_buckets()
        for bucket in buckets:
            print(bucket.name, bucket.creation_date)
    except Exception as e:
        err = e.message

    resp = { 'success': False, 'error': err } if err else { 'success': True }
    
    print(json.dumps(resp))
    return resp
