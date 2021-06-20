import random
import gym
from gym import spaces, logger
from gym.utils import seeding
import numpy as np
import time
from CalculateCost import CalculateCost


class AtomEnv(gym.Env):

    """

    Description:
        To find the structrue of low cost;
        3 atoms for now(num_atoms = 3);
        The tool to culculate the cost is provided by  College of Chemistry , Jilin University.

    Observation:
        Type: Discrete(3 * num_atoms)
        Num	Observation               Min             Max
        0	x of atom_1                0               1
        1	y of atom_1                0               1
        2	z of atom_1                0               1
        3	x of atom_2                0               1
        ...
        8   z of atom_3                0               1



    Actions:
        Type: Discrete(6 * num_atoms)
        Num	 Action
        0	 x of atom_1  +0.01
        1	 y of atom_1  +0.01
        2	 z of atom_1  +0.01
        3	 x of atom_2  +0.01
        ...
        17   z of atom_3  -0.01


        Note: action 0~2  belongs to atom_1;3~5 belongs to atom_2 ... They all +0.01.
              action 9~11 belongs to atom_1 ... 15~17 belongs to atom_3.They all -0.01


    Reward:
        Reward is -1 for every step taken, including the termination step

    Starting State:
        All observations are assigned a uniform random value in [0,1]

    Episode Termination:
        1.Episode length is greater than 1000.


    """



    metadata = {
        'render.modes': ['human', 'rgb_array'],
        'video.frames_per_second': 50
    }

    def __init__(self):
        self.num_atoms = 2
        self.num_coordinates = self.num_atoms * 3
        #满足条件的原子结构坐标保存在found_stru中
        self.found_stru = []
        self.done = False
        self.step_len = 1  # 步长因子，步长=step_len * 0.01

        #for an atom , it has three dimentional coordinates(x,y,z) ,and each x,y,z can "+0.01" or "-0.01"
        high = np.ones(self.num_coordinates)
        low = np.zeros(self.num_coordinates)

        self.action_space = spaces.Discrete(2 * self.num_coordinates)
        self.observation_space = spaces.Box(low, high, dtype                 =np.float32)

        self.seed()
        self.viewer = None
        self.state = None
        self.coor = None
        self.prcost = []

    def seed(self, seed=None):
        self.np_random, seed = seeding.np_random(seed)
        return [seed]

    def step(self, action):
        err_msg = "%r (%s) invalid" % (action, type(action))
        assert self.action_space.contains(action), err_msg
        done = False
        reward = -1
        # 求坐标改变前的cost，并与改变后的cost做对比，求判断reward
        prev_coordinates_for_cost = self.coor.reshape((self.num_atoms, 3))
        prev_cost = CalculateCost(prev_coordinates_for_cost)[3]

        changedCoordinate = action % self.num_coordinates
        if action < self.num_coordinates:
            self.coor[changedCoordinate] += 0.01 * self.step_len;
        else:
            self.coor[changedCoordinate] -= 0.01 * self.step_len
        while(self.coor[changedCoordinate]>1 or self.coor[changedCoordinate]<0):
            if self.coor[changedCoordinate] > 1:
                self.coor[changedCoordinate] -= 1
            elif self.coor[changedCoordinate] < 0:
                self.coor[changedCoordinate] += 1

        curr_coordinates_for_cost = self.coor.reshape((self.num_atoms, 3))
        curr_cost = CalculateCost(curr_coordinates_for_cost)[3]
        self.prcost.append(curr_cost)
        self.state = CalculateCost(curr_coordinates_for_cost)[0:3]

        # fo = open("cost_record_1.txt", "a")
        # cost_record = 'cost: ' + str(curr_cost) + '   reward:' + str(reward) + '   state:' + str(
        # self.state) + '   action: ' + str(action) + '  time: ' + str(time.asctime(time.localtime(time.time()))) + '\n'
        # fo.write(cost_record)
        # fo.close()
        self.step_len = 1
        if curr_cost < prev_cost:
            if prev_cost - curr_cost < 10:
                reward = 0.1
            elif prev_cost - curr_cost < 100:
                reward = 0.3
            else:
                print('prev_cost - curr_cost > 100!!  prev:', prev_cost, '   curr:', curr_cost)
                reward = 0.5
        elif prev_cost - curr_cost < -500:
                 print('prev_cost - curr_cost > -500!!  prev:', prev_cost, '   curr:', curr_cost)
                 reward = -5

        if curr_cost > 500 and random.random() > 0.2:
            self.step_len = 5
        elif curr_cost > 180 and random.random() > 0.2:
            self.step_len = 3

        if curr_cost < 20:
            reward = 1
        if curr_cost < 5:
            reward = 100
            done = True
            self.found_stru.append(self.state)
            fo = open("./found_stru.txt", "a+")
            goal = str(self.state) + '\t' + str(self.coor) + '\t' + str(curr_cost) + '\t' + str(action) +  '\t' + str(time.asctime(time.localtime(time.time()))) + '\r\n'
            print(goal)
            fo.write(goal)
            fo.close()


        return np.array(self.state), reward, done, {}

    def reset(self):
        #初始化state
        self.coor = self.np_random.uniform(low=0, high=1, size=(self.num_coordinates ,))
        self.coor= np.round(self.coor,4)
        c=self.coor
        c = c.reshape((self.num_atoms, 3))
        self.state = CalculateCost(c)[0:3]
        #print(self.state)
        return np.array(self.state)

    def render(self, mode='human'):
        return 0

    def close(self):
        if self.viewer:
            self.viewer.close()
            self.viewer = None
