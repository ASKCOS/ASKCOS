import os
import time
import makeit.global_config as gc

class MyLogger:
    '''
    Create logger. Four different levels of information output. A level 3 ("FATAL") log will exit the program.
    '''
    logFile = os.path.join(os.getcwd(),'log.txt')
    levels = {
        0:'INFO',
        1:'WARNING',
        2:'ERROR',
        3:'FATAL'
    }
    
    @staticmethod
    def initialize_logFile(ROOT = os.getcwd(), name = ''):
        MyLogger.logFile = os.path.join(ROOT, '{}_log.txt'.format(name))
        if os.path.isfile(MyLogger.logFile):
           os.remove(MyLogger.logFile)
        gc.time_zero = time.time() 
    
    @staticmethod
    def print_and_log(text,location,level=0):
        file = open(MyLogger.logFile,'a')
        time_elapsed = round((time.time() - gc.time_zero)*1000)/1000
        length = len(str(time_elapsed))
        pos = length
        max = 5
        
        while pos < 5:
            time_elapsed = str(time_elapsed)+'0'
            pos+=1
            
        time_pos = 20
        spaces = time_pos - len(location)
        pos = 0
        space = ''
        
        while pos < spaces:
            space+=' '
            pos+=1
            
        if spaces <= 0:
            location = location[0:len(location)+spaces-4]+'...'
            space = ' '
            
        print('{}@{}:{}[{}s]\t{}'.format(MyLogger.levels[level],location,space, time_elapsed, text))
        file.write('{}@{}:{}[{}s]\t{}\n'.format(MyLogger.levels[level],location,space, time_elapsed, text))
        
        if level == 3:
            quit()

