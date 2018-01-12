from __future__ import print_function

PLEVEL = 0
def vprint(level, txt, *args):
    if PLEVEL >= level:
        print(txt.format(*args))

def parity4(data):
    '''
    Thanks to http://www.dalkescientific.com/writings/diary/archive/2016/08/15/fragment_parity_calculation.html
    '''
    if data[0] < data[1]:
        if data[2] < data[3]:
            if data[0] < data[2]:
                if data[1] < data[2]:
                    return 0 # (0, 1, 2, 3) 
                else:
                    if data[1] < data[3]:
                        return 1 # (0, 2, 1, 3) 
                    else:
                        return 0 # (0, 3, 1, 2) 
            else:
                if data[0] < data[3]:
                    if data[1] < data[3]:
                        return 0 # (1, 2, 0, 3) 
                    else:
                        return 1 # (1, 3, 0, 2) 
                else:
                    return 0 # (2, 3, 0, 1) 
        else:
            if data[0] < data[3]:
                if data[1] < data[2]:
                    if data[1] < data[3]:
                        return 1 # (0, 1, 3, 2) 
                    else:
                        return 0 # (0, 2, 3, 1) 
                else:
                    return 1 # (0, 3, 2, 1) 
            else:
                if data[0] < data[2]:
                    if data[1] < data[2]:
                        return 1 # (1, 2, 3, 0) 
                    else:
                        return 0 # (1, 3, 2, 0) 
                else:
                    return 1 # (2, 3, 1, 0) 
    else:
        if data[2] < data[3]:
            if data[0] < data[3]:
                if data[0] < data[2]:
                    return 1 # (1, 0, 2, 3) 
                else:
                    if data[1] < data[2]:
                        return 0 # (2, 0, 1, 3) 
                    else:
                        return 1 # (2, 1, 0, 3) 
            else:
                if data[1] < data[2]:
                    return 1 # (3, 0, 1, 2) 
                else:
                    if data[1] < data[3]:
                        return 0 # (3, 1, 0, 2) 
                    else:
                        return 1 # (3, 2, 0, 1) 
        else:
            if data[0] < data[2]:
                if data[0] < data[3]:
                    return 0 # (1, 0, 3, 2) 
                else:
                    if data[1] < data[3]:
                        return 1 # (2, 0, 3, 1) 
                    else:
                        return 0 # (2, 1, 3, 0) 
            else:
                if data[1] < data[2]:
                    if data[1] < data[3]:
                        return 0 # (3, 0, 2, 1) 
                    else:
                        return 1 # (3, 1, 2, 0) 
                else:
                    return 0 # (3, 2, 1, 0) 
