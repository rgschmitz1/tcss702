from subprocess import run
from multiprocessing import cpu_count
from uuid import uuid4

class BWA:
    def process(self, mc, inputfile, bucket):
        # Download sample
        mc.fget_object(bucket, inputfile, f'/tmp/{inputfile}')

        # Uncompress sample
        uncompress_sample = run(f'tar -xf /tmp/{inputfile} -C /tmp', shell=True)
        print(uncompress_sample, flush=True)

        # Align sample and convert to bam
        inputfile = inputfile.split('.')[0]
        fastq = [f'{inputfile}_1.fq', f'{inputfile}_2.fq']
        aligned_bam = f'{inputfile}_realign-{uuid4()}.bam'

        threads = cpu_count()
        align_sample = run(f'bwa mem -t {threads} -T 0 /tmp/GRCh38.d1.vd1.fa /tmp/{fastq[0]} /tmp/{fastq[1]} | samtools sort -o /tmp/{aligned_bam}', shell=True)
        print(align_sample, flush=True)

        # Put aligned normal and tumor bam in minio
        mc.fput_object(bucket, aligned_bam, f'/tmp/{aligned_bam}')

        # Set aligned bam name for reference later
        self.__aligned_bam = aligned_bam

    def get_aligned_sample(self):
        return self.__aligned_bam
