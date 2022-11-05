# Serverless Benchmark functions
# source, https://github.com/spcl/serverless-benchmarks
import datetime

class SeBS():
    def __init__(self, mc):
        self.__mc = mc


    def dna_visualization(self, key, in_bucket, out_bucket):
        import io, json
        # using https://squiggle.readthedocs.io/en/latest/
        from squiggle import transform

        download_path = f'/tmp/{key}'

        download_begin = datetime.datetime.now()
        self.__mc.fget_object(in_bucket, key, download_path)
        download_stop = datetime.datetime.now()
        with open(download_path, "r") as f:
            data = f.read()

        process_begin = datetime.datetime.now()
        result = transform(data)
        process_end = datetime.datetime.now()

        upload_begin = datetime.datetime.now()
        buf = io.BytesIO(json.dumps(result).encode())
        buf.seek(0)
        self.__mc.put_object(out_bucket, f'transformed_{key}', buf, length=-1, part_size=10*1024*1024)
        upload_stop = datetime.datetime.now()
        buf.close()

        download_time = (download_stop - download_begin) / datetime.timedelta(microseconds=1)
        upload_time = (upload_stop - upload_begin) / datetime.timedelta(microseconds=1)
        process_time = (process_end - process_begin) / datetime.timedelta(microseconds=1)

        return {
            'result': {
                'bucket': out_bucket,
                'key': f'transformed_{key}'
            },
            'measurement': {
                'download_time': download_time,
                'upload_time': upload_time,
                'process_time': process_time
            }
        }


    def graph_bfs(self, size):
        import igraph

        graph_generating_begin = datetime.datetime.now()
        graph = igraph.Graph.Barabasi(size, 10)
        graph_generating_end = datetime.datetime.now()

        process_begin = datetime.datetime.now()
        result = graph.bfs(0)
        process_end = datetime.datetime.now()

        graph_generating_time = (graph_generating_end - graph_generating_begin) / datetime.timedelta(microseconds=1)
        process_time = (process_end - process_begin) / datetime.timedelta(microseconds=1)

        return {
            'result': result,
            'measurement': {
                'graph_generating_time': graph_generating_time,
                'process_time': process_time
            }
        }


    def graph_mst(self, size):
        import igraph

        graph_generating_begin = datetime.datetime.now()
        graph = igraph.Graph.Barabasi(size, 10)
        graph_generating_end = datetime.datetime.now()

        process_begin = datetime.datetime.now()
        result = graph.spanning_tree(None, False)[0]
        process_end = datetime.datetime.now()

        graph_generating_time = (graph_generating_end - graph_generating_begin) / datetime.timedelta(microseconds=1)
        process_time = (process_end - process_begin) / datetime.timedelta(microseconds=1)

        return {
            'result': result,
            'measurement': {
                'graph_generating_time': graph_generating_time,
                'process_time': process_time
            }
        }


    def graph_pagerank(self, size):
        import igraph

        graph_generating_begin = datetime.datetime.now()
        graph = igraph.Graph.Barabasi(size, 10)
        graph_generating_end = datetime.datetime.now()

        process_begin = datetime.datetime.now()
        result = graph.pagerank()[0]
        process_end = datetime.datetime.now()

        graph_generating_time = (graph_generating_end - graph_generating_begin) / datetime.timedelta(microseconds=1)
        process_time = (process_end - process_begin) / datetime.timedelta(microseconds=1)

        return {
            'result': result,
            'measurement': {
                'graph_generating_time': graph_generating_time,
                'process_time': process_time
            }
        }
