# Serverless Benchmark functions
# source, https://github.com/spcl/serverless-benchmarks
from datetime import datetime, timedelta

class SeBS():
    def dna_visualization(self, mc, key, in_bucket, out_bucket):
        from io import BytesIO
        from json import dumps
        # using https://squiggle.readthedocs.io/en/latest/
        from squiggle import transform
        from uuid import uuid4

        # Download sample from bucket
        download_begin = datetime.now()
        data = mc.get_object(in_bucket, key).read()
        download_end = datetime.now()

        # Transform sample
        process_begin = datetime.now()
        result = transform(data)
        process_end = datetime.now()

        # Upload sample to bucket
        self.__key_out = f'{key}-transformed-{uuid4()}'
        with BytesIO(dumps(result).encode()) as buf:
            buf.seek(0)
            upload_begin = datetime.now()
            mc.put_object(out_bucket, self.__key_out, buf, length=-1, part_size=10*1024*1024)
            upload_end = datetime.now()

        # Compute times
        download_time = (download_end - download_begin) / timedelta(microseconds=1)
        upload_time = (upload_end - upload_begin) / timedelta(microseconds=1)
        process_time = (process_end - process_begin) / timedelta(microseconds=1)

        return {
            'measurement': {
                'download_time': download_time,
                'upload_time': upload_time,
                'process_time': process_time
            }
        }


    def get_dna_visualization_output(self):
        return self.__key_out


    def _generate_barabasi_graph(self, size):
        from igraph import Graph

        graph_generating_begin = datetime.now()
        graph = Graph.Barabasi(n=size, m=10)
        graph_generating_end = datetime.now()
        graph_generating_time = (graph_generating_end - graph_generating_begin) / timedelta(microseconds=1)
        return (graph, graph_generating_time)


    def graph_bfs(self, size):
        graph, graph_generating_time = self._generate_barabasi_graph(size)

        process_begin = datetime.now()
        graph.bfs(vid=0)
        process_end = datetime.now()
        process_time = (process_end - process_begin) / timedelta(microseconds=1)

        return {
            'measurement': {
                'graph_generating_time': graph_generating_time,
                'process_time': process_time
            }
        }


    def graph_mst(self, size):
        graph, graph_generating_time = self._generate_barabasi_graph(size)

        process_begin = datetime.now()
        graph.spanning_tree(weights=None, return_tree=False)
        process_end = datetime.now()
        process_time = (process_end - process_begin) / timedelta(microseconds=1)

        return {
            'measurement': {
                'graph_generating_time': graph_generating_time,
                'process_time': process_time
            }
        }


    def graph_pagerank(self, size):
        graph, graph_generating_time = self._generate_barabasi_graph(size)

        process_begin = datetime.now()
        graph.pagerank()
        process_end = datetime.now()
        process_time = (process_end - process_begin) / timedelta(microseconds=1)

        return {
            'measurement': {
                'graph_generating_time': graph_generating_time,
                'process_time': process_time
            }
        }
