from minio import Minio
import json
import os
import subprocess
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

    # Download normal and tumor samples
    mc.fget_object('bwa-mem', 'normal.tar.xz', '/tmp/normal.tar.xz')
    mc.fget_object('bwa-mem', 'tumor.tar.xz', '/tmp/tumor.tar.xz')

    # Uncompress normal and tumor samples
    uncompress_normal_sample = subprocess.run('tar -xf /tmp/normal.tar.xz -C /tmp', shell=True)
    print(uncompress_normal_sample)
    uncompress_tumor_sample = subprocess.run('tar -xf /tmp/tumor.tar.xz -C /tmp', shell=True)
    print(uncompress_tumor_sample)

    # Extract reference sequence and BWA index
    uncompress_ref_seq = subprocess.run('tar -xf /home/app/function/genome/GRCh38.d1.vd1.fa.tar.gz -C /tmp', shell=True)
    print(uncompress_ref_seq)
    uncompress_bwa_idx = subprocess.run('tar -xf /home/app/function/genome/GRCh83.d1.vd1_BWA.tar.gz -C /tmp', shell=True)
    print(uncompress_bwa_idx)

    # Align normal
    align_normal = subprocess.run('bwa mem -t 8 -T 0 /tmp/GRCh38.d1.vd1.fa /tmp/normal_1.fq /tmp/normal_2.fq | samtools sort -o /tmp/normal_realign.bam', shell=True)
    print(align_normal)
    # Align tumor
    align_tumor = subprocess.run('bwa mem -t 8 -T 0 /tmp/GRCh38.d1.vd1.fa /tmp/tumor_1.fq /tmp/tumor_2.fq | samtools sort -o /tmp/tumor_realign.bam', shell=True)
    print(align_tumor)

    # Put aligned normal and tumor bam in minio
    mc.fput_object('bwa-mem', 'normal_realign.bam', '/tmp/normal_realign.bam')
    mc.fput_object('bwa-mem', 'tumor_realign.bam', '/tmp/tumor_realign.bam')

    # Collect data
    inspector = Inspector()
    inspector.inspectAll()

    # Add custom message and finish the function
#    if "startWallClock" in event:
#        inspector.addAttribute("startWallClock", event['startWallClock'])

    body = json.loads(event.body)

    inspector.inspectAllDeltas()
    # Include functionName
#    inspector.addAttribute("functionName", fn[body['fn']].__name__)

    iret = inspector.finish()
    ret = {
        "status": 200,
        "body": iret
    }
    return ret
