from cpp_backend import *
b = CppBackend()

class Action:
    def __init__(self, action_id, args):
        self.action_id = action_id
        self.args = args

    def __str__(self):
        return f"iid: {self.action_id} args {self.args}"

    def to_json(self):
        d = { # # # # # # #
            "query_type" : "action",
            "action_id" : self.action_id,
            "params" : self.args
        }
        return dict_to_query(d)

    def __le__(self, other):
        return (self.action_id, self.args) <= (other.action_id, other.args)

    def __ge__(self, other):
        return (self.action_id, self.args) >= (other.action_id, other.args)



class BaseEnv:
    def __init__(self):
        self.backend = None
        self.max_time = 3
        self.reset()

    def step(self, action):
        if not isinstance(action, Action):
            iid, args = action
            action = Action(iid, args)

        self.backend.send(action.to_json())

        response = self.backend.get(False)
        response = eval(response)

        observation = self._get_observation_slow()

        shape = observation['info']['shape']
        reward = shape[0] - self.last_observation['info']['shape'][0]
        self.last_observation = observation

        self.time += 1
        done = self.time >= self.max_time

        return (self.last_observation, reward, done, None)

    def observation(self): return self.last_observation

    def reset(self):
        self.time = 0
        if self.backend is None:
            self.backend = CppBackend()
        else:
            self.backend.reset()
        self.last_observation = self._get_observation_slow()
        return self.last_observation

    def _get_observation_slow(self):
        self.backend.send(dict_to_query({"query_type" : "graph"}))
        return self.backend.get()

    def render(self):
        # with open('graph.json', 'w+') as file:
        #     print(graph, file=file)
        # with open('storage.json', 'w+') as file:
        #     print(storage, file=file)
        pass

    def save(self):
        self.saved_time = self.time
        self.saved_last_obs = self.last_observation
        self.backend.send(dict_to_query({"query_type" : "save"}))

    def load(self):
        self.time = self.saved_time
        self.last_observation = self.saved_last_obs
        self.backend.send(dict_to_query({"query_type" : "load"}))


    def close(self):
        pass

    def seed(self, seed):
        pass

