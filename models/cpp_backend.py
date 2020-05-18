from subprocess import Popen, PIPE
import numpy as np
import json


def dict_to_query(d):
    j = json.dumps(d, separators=(',', ':')) + '\n'
    j = str.encode(j)
    return j


class CppBackend:
    executable_path = '../cpp_backend/cmake-build-debug/cpp_backend'

    def __init__(self):
        self.p = Popen([CppBackend.executable_path], shell=True, stdout=PIPE, stdin=PIPE)

    def __del__(self):
        self.p.kill()

    def send(self, s):
        self.p.stdin.write(s)
        self.p.stdin.flush()

    def get(self, decode_json=True):
        line = self.p.stdout.readline()
        result = line.decode("utf-8")
        if decode_json: result = json.loads(result)
        # print("recieved:", result)
        return result

    def reset(self):
        # raise NotImplementedError()
        self.send(dict_to_query({"query_type" : "reset"}))

    def close(self):
        self.p.close()

