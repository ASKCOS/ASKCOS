import os
import time
import makeit.global_config as gc


class MyLogger:
    '''
    Create logger. Four different levels of information output. A level 3 ("FATAL") 
    log will exit the program.
    '''
    logFile = os.path.join(os.path.dirname(
        os.path.dirname(os.path.dirname(__file__))), 'log.txt')
    levels = {
        0: 'INFO',
        1: 'WARN',
        2: 'ERROR',
        3: 'FATAL'
    }

    @staticmethod
    def initialize_logFile(ROOT=os.path.dirname(os.path.dirname(os.path.dirname(__file__))), name=''):
        MyLogger.logFile = os.path.join(ROOT, '{}_log.txt'.format(name))
        if os.path.isfile(MyLogger.logFile):
            os.remove(MyLogger.logFile)
        gc.time_zero = time.time()

    @staticmethod
    def print_and_log(text, location, level=0):
        file = open(MyLogger.logFile, 'a')
        time_elapsed = time.time() - gc.time_zero

        outstr = '{:25s}: [{:04.3f}s]\t{}'.format('{}@{}'.format(MyLogger.levels[level], location)[:25],
            time_elapsed, text)
        print(outstr)
        file.write(outstr)
        file.write('\n')

        if level == 3:
            quit()
