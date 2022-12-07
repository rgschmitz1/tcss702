"""
BWA utility wrapper

@author: Bob Schmitz
@version: 2022-11-29
"""
from shutil import rmtree
from subprocess import run
from tempfile import mkdtemp
from uuid import uuid4

class BWA:
    def process(self, mc, inputfile: str, bucket: str, thread_cnt: int=1) -> None:
        # Create a temp directory
        tmp = mkdtemp()

        # Download sample
        try:
            mc.fget_object(bucket, inputfile, f'{tmp}/{inputfile}')
        except:
            self._cleanup(tmp)
            return {'status': 500, 'body': f'failed to download sample {inputfile}'}

        # Uncompress sample
        uncompress_sample = run(f'tar -xf {tmp}/{inputfile} -C {tmp}', shell=True)
        print(uncompress_sample, flush=True)
        if uncompress_sample.returncode:
            self._cleanup(tmp)
            return {'status': uncompress_sample.returncode, 'body': 'failed to uncompress sample'}

        inputfile = inputfile.split('.')[0]
        fastq = [f'{inputfile}_1.fq', f'{inputfile}_2.fq']

        aligned_bam = f'{inputfile}_realign-{uuid4()}.bam'

        # Align sample and convert to bam
        align_sample = run(f'bwa mem -t {thread_cnt} -T 0 /tmp/GRCh38.d1.vd1.fa {tmp}/{fastq[0]} {tmp}/{fastq[1]} | samtools sort -o {tmp}/{aligned_bam}', shell=True)
        print(align_sample, flush=True)
        if align_sample.returncode:
            self._cleanup(tmp)
            return {'status': align_sample.returncode, 'body': 'failed running bwa mem alignment'}

        # Put aligned bam in minio
        try:
            mc.fput_object(bucket, aligned_bam, f'{tmp}/{aligned_bam}')
        except:
            self._cleanup(tmp)
            return {'status': 500, 'body': f'failed to upload aligned sample {aligned_bam}'}

        self._cleanup(tmp)

        # Return aligned sample object name
        return {'status': 0, 'body': aligned_bam}


    def _cleanup(self, tmp) -> None:
        # Remove temp directory
        try:
            rmtree(tmp)
        except:
            print(f"encountered and error removing '{tmp}'", flush=True)
