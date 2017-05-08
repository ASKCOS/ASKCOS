from __future__ import absolute_import, unicode_literals, print_function
from .celery import app
import time 

@app.task
def print_and_wait(x, y):
    print('{} + {} = ...?'.format(x, y))
    time.sleep(5)
    print('{}'.format(x+y))
    return x+y
