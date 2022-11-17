# Serverless Benchmark functions
# source, https://github.com/spcl/serverless-benchmarks
import datetime

class SeBS():
    def dna_visualization(self, mc, key, in_bucket, out_bucket):
        import io, json
        # using https://squiggle.readthedocs.io/en/latest/
        from squiggle import transform

        # Download sample from bucket
        download_path = f'/tmp/{key}'
        download_begin = datetime.datetime.now()
        mc.fget_object(in_bucket, key, download_path)
        download_end = datetime.datetime.now()
        with open(download_path, "r") as f:
            data = f.read()

        # Transform sample
        process_begin = datetime.datetime.now()
        result = transform(data)
        process_end = datetime.datetime.now()

        # Upload sample to bucket
        upload_begin = datetime.datetime.now()
        buf = io.BytesIO(json.dumps(result).encode())
        buf.seek(0)
        mc.put_object(out_bucket, f'transformed_{key}', buf, length=-1, part_size=10*1024*1024)
        upload_end = datetime.datetime.now()
        buf.close()

        # Compute times
        download_time = (download_end - download_begin) / datetime.timedelta(microseconds=1)
        upload_time = (upload_end - upload_begin) / datetime.timedelta(microseconds=1)
        process_time = (process_end - process_begin) / datetime.timedelta(microseconds=1)

        return {
            'measurement': {
                'download_time': download_time,
                'upload_time': upload_time,
                'process_time': process_time
            }
        }


    def _generate_barabasi_graph(self, size):
        from igraph import Graph

        graph_generating_begin = datetime.datetime.now()
        graph = Graph.Barabasi(n=size, m=10)
        graph_generating_end = datetime.datetime.now()
        graph_generating_time = (graph_generating_end - graph_generating_begin) / datetime.timedelta(microseconds=1)
        return (graph, graph_generating_time)



    def graph_bfs(self, size):
        graph, graph_generating_time = self._generate_barabasi_graph(size)

        process_begin = datetime.datetime.now()
        graph.bfs(vid=0)
        process_end = datetime.datetime.now()
        process_time = (process_end - process_begin) / datetime.timedelta(microseconds=1)

        return {
            'measurement': {
                'graph_generating_time': graph_generating_time,
                'process_time': process_time
            }
        }


    def graph_mst(self, size):
        graph, graph_generating_time = self._generate_barabasi_graph(size)

        process_begin = datetime.datetime.now()
        graph.spanning_tree(weights=None, return_tree=False)
        process_end = datetime.datetime.now()
        process_time = (process_end - process_begin) / datetime.timedelta(microseconds=1)

        return {
            'measurement': {
                'graph_generating_time': graph_generating_time,
                'process_time': process_time
            }
        }


    def graph_pagerank(self, size):
        graph, graph_generating_time = self._generate_barabasi_graph(size)

        process_begin = datetime.datetime.now()
        graph.pagerank()
        process_end = datetime.datetime.now()
        process_time = (process_end - process_begin) / datetime.timedelta(microseconds=1)

        return {
            'measurement': {
                'graph_generating_time': graph_generating_time,
                'process_time': process_time
            }
        }
