#!/usr/bin/python3
import torch
import torch.nn as nn
from torch.distributions.categorical import Categorical
from torch.optim import Adam
import numpy as np
import gym
from gym.spaces import Discrete, Box

from base_env import BaseEnv

def reward_to_go(rews):
    n = len(rews)
    rtgs = np.zeros_like(rews)
    for i in reversed(range(n)):
        rtgs[i] = rews[i] + (rtgs[i+1] if i+1 < n else 0)
    return rtgs


######################
import torch
import torch.nn.functional as F
from torch_geometric.nn import GATConv, MessagePassing
from cg_conv import CGConv
from torch_geometric.data import Data

class Net(torch.nn.Module):
    def __init__(self):
        super(Net, self).__init__()
        channels = 5
        self.lin1 = torch.nn.Linear(2, channels)
        self.conv1 = CGConv(channels=channels, dim=1)
        self.conv2 = CGConv(channels=channels, dim=1)
        self.lin2 = torch.nn.Linear(channels, 2)

    def forward(self, data):
        x = data.x
        x = self.lin1(x)
        x = self.conv1(x, data.edge_index, data.edge_attr)
        x = self.conv2(x, data.edge_index, data.edge_attr)

        x = self.lin2(x)
        return x
##########################

def obs_to_data(obs):
    num_edges = len(obs['e'])
    num_vertices = len(obs['v'])

    d = Data()
    d.edge_index = torch.zeros(2, num_edges).long()
    d.edge_attr = torch.zeros(num_edges, 1).float()

    d.x = torch.zeros(num_vertices, 2).float()
    d.y = torch.zeros(num_vertices).long()

    for i, (id, tag, _) in enumerate(obs['v']):
        d.x[id][0] = tag
        d.x[id][1] = id
        # d.y[id] = tag

    for i, (v, u, tag) in enumerate(obs['e']):
        d.edge_index[0][i] = v
        d.edge_index[1][i] = u
        d.edge_attr[i][0] = tag

    return d

def train(lr=1e-2, epochs=1000, batch_size=500, render=False):
    env = BaseEnv()
    logits_net = Net()
    print('start')

    # make function to compute action distribution
    def get_policy(obs):
        data = obs_to_data(obs)
        logits0 = logits_net(data)[:,0]
        logits1 = logits_net(data)[:,1]
        return Categorical(logits=logits0.reshape(-1)), Categorical(logits=logits1.reshape(-1))

    # make action selection function (outputs int actions, sampled from policy)
    def get_action(obs):
        # iid = int(np.random.choice([0, 1, 2]))
        iid = int(np.random.choice([0]))

        p1, p2 = get_policy(obs)
        a1 = p1.sample().item()
        a2 = p2.sample().item()

        return (iid, (a1, a2))

    # make loss function whose gradient, for the right data, is policy gradient

    def compute_loss(batch_obs, batch_acts, batch_weights):
        loss = 0
        for obs, act, weights in zip(batch_obs, batch_acts, batch_weights):
            p1, p2 = get_policy(obs)
            a1 = torch.as_tensor(act[1][0], dtype=torch.float32)
            a2 = torch.as_tensor(act[1][1], dtype=torch.float32)
            x = p1.log_prob(a1) + p2.log_prob(a2)
            loss += x * weights

        loss = -loss
        loss /= len(batch_obs)
        return loss

    # make optimizer
    optimizer = Adam(logits_net.parameters(), lr=lr)

    # for trainint policy
    def train_one_epoch():
        batch_obs = []
        batch_acts = []
        batch_weights = []
        batch_rets = []
        batch_lens = []

        obs = env.reset()
        done = False
        ep_rews = []

        # render first episode of each epoch
        finished_rendering_this_epoch = False

        while True:
            if (not finished_rendering_this_epoch) and render: env.render()
            batch_obs.append(obs.copy())
            act = get_action(obs)
            obs, rew, done, _ = env.step(act)
            # print(f'step: {act}')
            batch_acts.append(act)
            ep_rews.append(rew)

            if done:
                # print('last:', obs['info']['shape'])
                ep_ret, ep_len = sum(ep_rews), len(ep_rews)
                batch_rets.append(ep_ret)
                batch_lens.append(ep_len)

                # the weight for each logprob(a_t|s_t) is reward-to-go from t
                batch_weights += list(reward_to_go(ep_rews))

                obs, done, ep_rews = env.reset(), False, []
                finished_rendering_this_epoch = True
                if len(batch_obs) > batch_size: break

        # take a single policy gradient update step
        optimizer.zero_grad()
        batch_loss = compute_loss(batch_obs, batch_acts, batch_weights)
        batch_loss.backward()
        optimizer.step()
        return batch_loss, batch_rets, batch_lens

    # training loop
    for i in range(epochs):
        batch_loss, batch_rets, batch_lens = train_one_epoch()
        print('epoch: %3d \t loss: %.3f \t return: %.3f \t median_return: %.3f \t ep_len: %.3f'%
                (i, batch_loss, np.mean(batch_rets), np.median(batch_rets),np.mean(batch_lens)))

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--env_name', '--env', type=str, default='CartPole-v0')
    parser.add_argument('--render', action='store_true')
    parser.add_argument('--lr', type=float, default=1e-2)
    args = parser.parse_args()
    # print('\nUsing reward-to-go formulation of policy gradient.\n')
    train(render=args.render, lr=args.lr)
