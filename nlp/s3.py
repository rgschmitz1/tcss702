import logging
import boto3
from botocore.client import Config
from botocore.exceptions import ClientError
import os
"""
Utility module for S3 buckets

@author: Bob Schmitz
"""

hostname=os.getenv['minio_hostname']
_s3_client = boto3.client('s3',
    endpoint_url=f'http://{hostname}',
    aws_access_key_id=os.getenv['minio_access_key'],
    aws_secret_access_key=os.getenv['minio_secret_key'],
    config=Config(signature_version='s3v4'),
    region_name='us-east-1'
)

def s3_delete(bucket, object_name):
    """Download a file to an S3 bucket

    :param bucket: Bucket to download to
    :param object_name: S3 object name.
    :return: True if file is deleted, else False
    """
    # Delete the file
    try:
        response = _s3_client.delete_object(bucket, object_name)
    except ClientError as e:
        logging.error(e)
        return False
    return True


def s3_download(file_name, bucket, object_name=None):
    """Download a file to an S3 bucket

    :param file_name: File to download
    :param bucket: Bucket to download to
    :param object_name: S3 object name. If not specified then file_name is used
    :return: True if file was downloaded, else False
    """
    # If S3 object_name was not specified, use file_name
    if object_name is None:
        object_name = os.path.basename(file_name)

    # Download the file
    try:
        response = _s3_client.download_file(bucket, object_name, file_name)
    except ClientError as e:
        logging.error(e)
        return False
    return True


def s3_upload_file(file_name, bucket, object_name=None):
    """Upload a file to an S3 bucket

    :param file_name: File to upload
    :param bucket: Bucket to upload to
    :param object_name: S3 object name. If not specified then file_name is used
    :return: True if file was uploaded, else False
    """
    # If S3 object_name was not specified, use file_name
    if object_name is None:
        object_name = os.path.basename(file_name)

    # Upload the file
    try:
        response = _s3_client.upload_file(file_name, bucket, object_name)
    except ClientError as e:
        logging.error(e)
        return False
    return True


def s3_upload_object(buffer_object, bucket, object_name):
    """Upload a file to an S3 bucket

    :param buffer_object: buffer object to upload
    :param bucket: Bucket to upload to
    :param object_name: S3 object name.
    :return: True if file was uploaded, else False
    """
    try:
        _s3_client.put_object(Body=buffer_object, Bucket=bucket, Key=object_name)
    except ClientError as e:
        logging.error(e)
        return False
    return True
