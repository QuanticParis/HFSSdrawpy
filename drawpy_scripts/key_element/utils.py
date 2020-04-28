from functools import wraps

def _moved(func, previous_pos, previous_ori, *args, **kwargs):
    #1 We need to keep track of the entities created during the execution of a function
    list_entities_old =args[0].modelentities_to_move()
    list_ports_old = args[0].ports_to_move()
    args[0].append_lists()

    #2 We store the current position on the chip where the piece is to be drawn.
    #2 It is modified during the execution of innner functions, which forces us to store them BEFORE the execution
    pos = args[0].current_pos
    angle = args[0].current_ori

    if len(pos)==2:
        pos.append(0)
    #3 We execute the actual function
    result = func(*args, **kwargs)

    #4 We move the entity that were created by the last function
    list_entities_new = args[0].modelentities_to_move()
    list_ports_new = args[0].ports_to_move()

#    print("to_move", [i.name for i in list_entities_new])

    #5 We move the entities_to_move with the right operation
    if len(list_entities_new)>0:
        args[0].rotate(list_entities_new, angle=angle)
        args[0].translate(list_entities_new, vector=[pos[0], pos[1], pos[2]])

    if len(list_ports_new)>0:
        args[0].rotate_port(list_ports_new, angle)
        args[0].translate_port(list_ports_new, vector=[pos[0], pos[1], pos[2]])

    #6 We empty a part of the 'to_move' dictionnaries
    if isinstance(list_entities_old[-1], list):
        a = list_entities_old.pop(-1)
        for entity in a:
            list_entities_old.append(entity)
    else:
        args[0].reset()

    if isinstance(list_ports_old[-1], list):
        a = list_ports_old.pop(-1)
        for entity in a:
            list_ports_old.append(entity)
    else:
        args[0].reset()

    #7 We reset the current_coor to the last variables saved
    args[0].set_coor(previous_pos, previous_ori)
    return result

def move(func):
    '''
    Decorator which moves the KeyElements and CustomElements (rotation+translation) with the parameters given by get_current_coor

    Input : func to move
    Output : moved is the composition of func with the rotation+translation of the ModeleEntities created during its execution
    '''
    @wraps(func)
    # At the begining of the execution, we decide that all elements created go to instances_to_move
    def moved(*args, **kwargs):
        # args[0] = body
        previous_pos = args[0].current_pos
        previous_ori = args[0].current_ori
#        print(previous_pos, previous_ori)
        return _moved(func, previous_pos, previous_ori, *args, **kwargs)
    return moved