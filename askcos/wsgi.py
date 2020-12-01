import os
import sys

sys.path.append('/home/ubuntu/miniconda/envs/askcos/lib/python2.7/site-packages/')
sys.path.append('/home/ubuntu/ASKCOS/Make-It')
sys.path.append('/home/ubuntu/ASKCOS/ASKCOS_Website/')
sys.path.append('/home/ubuntu/ASKCOS/RDChiral/')
sys.path.append('/home/ubuntu/ASKCOS/SCScore')

os.environ['DJANGO_SETTINGS_MODULE'] = 'askcos_site.settings'

if 'PYTHONPATH' not in os.environ:
    os.environ['PYTHONPATH'] = ''
os.environ['PYTHONPATH'] = '/home/ubuntu/ASKCOS/ASKCOS_Website:/home/ubuntu/ASKCOS/Make-It:home/ubuntu/miniconda/envs/askcos/lib/python2.7/site-packages:' + os.environ['PYTHONPATH']
if 'PATH' not in os.environ:
    os.environ['PATH'] = ''
os.environ['PATH'] = r'/home/ubuntu/miniconda/envs/askcos/bin:/usr/bin:' + os.environ['PATH']
# NOTE: LD_LIBRARY_PATH set in /etc/apache2/envvars

from django.core.wsgi import get_wsgi_application
application = get_wsgi_application()

# Load libcairo
from ctypes import cdll
cdll.LoadLibrary('libcairo.so.2')
import cairo 
cairo_dir = dir(cairo)
if len(cairo_dir) < 20:
    print('######################################')
    print('## CAIRO INSTALLATION LIKELY BROKEN ##')
    print('######################################')
