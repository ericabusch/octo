#!/usr/bin/env python
# coding: utf-8


import numpy as np

outside_rock_prob = 0.5
distribution_outer = [outside_rock_prob, 1 - outside_rock_prob]
distribution_inner = [1/3, 1/3, 1/3]
num_food = 2
options = ['rock  ', 'blank  ', 'unique  ', 'food']
food_options = ['mussel', 'clam', 'shrimp']


def randomize_tank(length, width):
    total = length * width
    n_outer = length * 2 + width * 2 - 4
    n_inner = total - n_outer
    x=0
    outer_loop = list(np.random.choice(options[:2], size=n_outer, p=distribution_outer))
    for i in outer_loop:
        if str(i)== 'rock  ':
            x+=1

    
    inner_loop = list(np.random.choice(options[:3], size=n_inner-num_food, p=distribution_inner))
    inner_loop += list((np.random.choice(food_options, size=num_food)))
    np.random.shuffle(inner_loop)
    
    print('Outside:')
    for i in range(len(outer_loop)):
        if i < 10:
            print(f'{i}    {outer_loop[i]}')
        else:
            print(f'{i}   {outer_loop[i]}')
        
    print('\nInside')
    for i in range(len(inner_loop)):
        if i < 10:
            print(f'{i}    {inner_loop[i]}')
        else:
            print(f'{i}   {inner_loop[i]}')
    
    x,y = 0,0
    for i in outer_loop:
        if str(i) == 'rock  ':
            x+=1
    print('number of rocks in outer loop is: ' +str(x))

    x,y = 0,0
    for i in inner_loop:
        if str(i) == 'rock  ':
            x+=1
        elif str(i) == 'unique  ':
            y+=1
    print(f'number of rocks and unique in inner loop is: {x} and {y} respectively\n')
    
    print_tank(outer_loop, inner_loop, length, width)


def print_tank(outer, inner, length, width):
    tank = [['X']*length]*width
    tank[0], tank[-1] = outer[:length], outer[-length:]
    
    for x in range(1, width-1):
        tank[x][0] = 'Y'
        tank[x][-1] = 'Z'
        
    tank = np.array(tank)
    np.place(tank, tank == 'X', inner)
    np.place(tank, tank == 'Z', outer[length+(width//2):len(outer)-length])
    np.place(tank, tank == 'Y', outer[length:length+(width//2)])
    np.place(tank, tank =='rock ', 'rock')
    np.place(tank, tank == 'blank ', 'blank')
    
    print('\nyour tank looks like:\n')
    for row in tank:
        print(row)


randomize_tank(5,4)




