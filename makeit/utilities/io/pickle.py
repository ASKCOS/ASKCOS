import sys 
import six
from six.moves import cPickle as pickle

def convert_pickled_bytes_2_to_3(data):
    if isinstance(data, bytes):  return data.decode()
    if isinstance(data, dict):   return dict(map(convert_pickled_bytes_2_to_3, data.items()))
    if isinstance(data, tuple):  return tuple(map(convert_pickled_bytes_2_to_3, data))
    if isinstance(data, list):   return list(map(convert_pickled_bytes_2_to_3, data))
    return data


def load(file):
    if sys.version_info[0] < 3:
        return pickle.load(file)
    else:
        return convert_pickled_bytes_2_to_3(pickle.load(file, encoding='bytes'))

def dump(data, file, *args, **kwargs):
    '''Always use protocol 2 for backwards compatibility!'''
    pickle.dump(data, file, 2) # note: always use protocol 2 for backward compatibility