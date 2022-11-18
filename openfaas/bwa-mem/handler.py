from minio import Minio
import json
import os
import subprocess
from .SAAF import Inspector

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

    # Get minio client object to download/upload data
    mc = minio_client()

    # Download sample
    mc.fget_object(bucket, inputfile, f'/tmp/{inputfile}')

    # Uncompress sample
    uncompress_sample = subprocess.run(f'tar -xf /tmp/{inputfile} -C /tmp', shell=True)
    print(uncompress_sample, flush=True)

    # Extract BWA index
    uncompress_bwa_idx = subprocess.run('tar -xf /home/app/function/genome/GRCh38.d1.vd1_BWA.tar.zst -C /tmp', shell=True)
    print(uncompress_bwa_idx, flush=True)

    # Create an empty reference file in tmp, bwa will reference this for name only
    with open('/tmp/GRCh38.d1.vd1.fa', 'w') as f:
        pass

    # Align sample and convert to bam
    inputfile = inputfile.split('.')[0]
    fastq = [f'{inputfile}_1.fq', f'{inputfile}_2.fq']
    aligned_bam = f'{inputfile}_realign.bam'
    align_sample = subprocess.run(f'bwa mem -t 8 -T 0 /tmp/GRCh38.d1.vd1.fa /tmp/{fastq[0]} /tmp/{fastq[1]} | samtools sort -o /tmp/{aligned_bam}', shell=True)
    print(align_sample, flush=True)

    # Put aligned normal and tumor bam in minio
    mc.fput_object(bucket, aligned_bam, f'/tmp/{aligned_bam}')

    # Collect inspector deltas
    inspector.inspectAllDeltas()

    # Include functionName
    inspector.addAttribute("functionName", 'bwa-mem')

    iret = inspector.finish()
    ret = {
        "status": 200,
        "body": iret
    }
    return ret
